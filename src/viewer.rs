use crate::buffer::DisplayBuffer;
use crate::color::ColorScheme;
use crate::reader::{open_source, RecordSource};
use crate::record::Record;
use crate::search::{MatchScratch, Pattern, SearchIndex, SearchScope};
use crate::stats::{FastqStats, StatsAccumulator, AVG_QUAL_BIN_WIDTH, AVG_QUAL_NUM_BINS};
use crate::utils::{
    calculate_record_lines, determine_min_max_phred, translate_frames,
    translate_frames_with_quality, PhredRange, ReadOrientation,
};
use crate::worker::Worker;
use anyhow::Result;
use bio::alphabets::dna;
use nix::poll::PollFlags;
use nix::poll::{poll, PollFd};
use num_format::{Locale, ToFormattedString};
use ratatui::{
    backend::TermionBackend,
    layout::*,
    style::{Color, Style},
    text::{Line, Span},
    widgets::{Bar, BarChart, BarGroup, *},
    Terminal,
};

use signal_hook::consts::SIGINT;
use signal_hook::flag;
use std::cmp::min;
use std::io::stdin;
use std::os::unix::io::AsRawFd;
use std::sync::mpsc::{self, Receiver};
use std::sync::{
    atomic::{AtomicBool, AtomicU64, Ordering},
    Arc, Mutex,
};
use std::thread;
use std::time::Duration;
use std::{fs::File, io::stdout};
use std::{io::Write, os::fd::BorrowedFd};
use termion::{
    event::Key,
    input::TermRead,
    raw::{IntoRawMode, RawTerminal},
    screen::{AlternateScreen, IntoAlternateScreen},
};

// Constants for buffer and scrolling behavior
const BUFFER_WINDOW_SIZE: usize = 1000;
const PAGE_SCROLL_AMOUNT: usize = 10;

pub struct TuiViewer {
    terminal: Terminal<TermionBackend<AlternateScreen<RawTerminal<std::io::Stdout>>>>,
    /// One entry per input file: R1, plus R2 in paired mode.
    mates: Vec<Mate>,
    current_position: u64,
    horizontal_offset: usize,
    no_wrap: bool,
    show_quality: bool,
    show_stats: bool,
    show_help: bool,
    orientation: ReadOrientation,
    show_translation: bool,
    stats_scroll: usize,
    stats_tab: usize,
    help_scroll: usize,
    phred_range: PhredRange,
    color_scheme: ColorScheme,
    search_mode: bool,
    search_input: String,
    search_query: String,
    search_pattern: Option<Pattern>,
    search_include_rc: bool,
    /// Which record fields the pattern matches: sequence, name, or both.
    search_scope: SearchScope,
    /// True when the scope was chosen automatically because the query looks
    /// like a read name rather than a nucleotide sequence.
    search_scope_auto: bool,
    /// 1-based number of the current match among all matches.
    search_match_k: u64,
}

#[derive(Clone, Copy)]
enum MatchDirection {
    AtOrAfter,
    Before,
}

/// How often background scans publish their progress.
const COUNT_PROGRESS_EVERY: u64 = 100_000;
const STATS_PROGRESS_EVERY: u64 = 100_000;
const SEARCH_PROGRESS_EVERY: u64 = 10_000;

/// Everything that belongs to one input file: R1, or R2 in paired mode.
struct Mate {
    input: Input,
    buffer: DisplayBuffer,
    /// Number of reads in the file, from a fast background count. Stdin
    /// cannot be read twice, so for stdin only loaded reads are known.
    total_reads: Arc<AtomicU64>,
    count_done: Arc<AtomicBool>,
    counter: Worker,
    stats: Arc<Mutex<FastqStats>>,
    stats_worker: Worker,
    /// Stdin statistics, accumulated on the main thread from loaded reads.
    stdin_stats: Option<StatsAccumulator>,
    search_index: Arc<Mutex<SearchIndex>>,
    search_indexed_up_to: u64,
    search_worker: Worker,
    search_scanned: Arc<AtomicU64>,
    search_done: Arc<AtomicBool>,
}

impl Mate {
    fn open(path: String, reference: Option<String>) -> Result<Self> {
        let buffer = DisplayBuffer::new(&path, reference.as_deref())?;
        Ok(Self {
            input: Input { path, reference },
            buffer,
            total_reads: Arc::new(AtomicU64::new(0)),
            count_done: Arc::new(AtomicBool::new(false)),
            counter: Worker::default(),
            stats: Arc::new(Mutex::new(FastqStats::default())),
            stats_worker: Worker::default(),
            stdin_stats: None,
            search_index: Arc::new(Mutex::new(SearchIndex::new())),
            search_indexed_up_to: 0,
            search_worker: Worker::default(),
            search_scanned: Arc::new(AtomicU64::new(0)),
            search_done: Arc::new(AtomicBool::new(true)),
        })
    }

    fn is_stdin(&self) -> bool {
        self.input.path == "-"
    }

    fn reads_loaded(&self) -> u64 {
        self.buffer.reads.read().unwrap().len() as u64
    }

    /// Reads in the input and whether that number is final. For stdin it is
    /// the number loaded so far until the end of the input has been reached.
    fn total_reads(&self) -> (u64, bool) {
        if self.is_stdin() {
            (self.reads_loaded(), self.buffer.is_fully_loaded())
        } else {
            (
                self.total_reads.load(Ordering::Relaxed),
                self.count_done.load(Ordering::Relaxed),
            )
        }
    }

    /// Count the reads of a file in the background. Counting does not depend
    /// on any view setting, so this runs once.
    fn start_counter(&mut self) {
        if self.is_stdin() {
            return;
        }
        let input = self.input.clone();
        let total = Arc::clone(&self.total_reads);
        let done = Arc::clone(&self.count_done);
        self.counter = Worker::spawn(move |stop| count_reads(&input, &total, &done, stop));
    }

    /// Start computing statistics from scratch with `phred_range`.
    fn start_stats(&mut self, phred_range: PhredRange) {
        self.stats_worker.stop();
        *self.stats.lock().unwrap() = FastqStats::default();
        if self.is_stdin() {
            self.stdin_stats = Some(StatsAccumulator::new(phred_range));
            self.update_stdin_stats();
            return;
        }
        let input = self.input.clone();
        let stats = Arc::clone(&self.stats);
        self.stats_worker =
            Worker::spawn(move |stop| scan_stats(&input, phred_range, &stats, stop));
    }

    /// Add reads loaded from stdin since the last call to its statistics.
    /// Only the new reads are visited, so this is cheap to call every frame.
    fn update_stdin_stats(&mut self) {
        let Some(acc) = self.stdin_stats.as_mut() else {
            return;
        };
        let reads = self.buffer.reads.read().unwrap();
        let seen = acc.reads() as usize;
        let complete = self.buffer.is_fully_loaded();
        if seen == reads.len() && self.stats.lock().unwrap().scanned_all == complete {
            return;
        }
        for record in &reads[seen..] {
            acc.add(record);
        }
        acc.snapshot_into(&mut self.stats.lock().unwrap(), complete);
    }

    fn reset_search(&mut self) {
        // Stop the worker before clearing the index; otherwise it can set
        // bits from the old pattern/scope after the reset.
        self.search_worker.stop();
        self.search_index.lock().unwrap().reset();
        self.search_indexed_up_to = 0;
        self.search_done.store(true, Ordering::SeqCst);
        self.search_scanned.store(0, Ordering::SeqCst);
    }

    /// Mark matches among the reads loaded since the last call.
    fn index_loaded(&mut self, pattern: Option<&Pattern>, scope: SearchScope) {
        let reads = self.buffer.reads.read().unwrap();
        let len = reads.len() as u64;
        if let Some(pattern) = pattern {
            let mut index = self.search_index.lock().unwrap();
            let mut scratch = MatchScratch::new();
            for record_idx in self.search_indexed_up_to..len {
                if scope.matches_with(pattern, &reads[record_idx as usize], &mut scratch) {
                    index.set(record_idx);
                }
            }
        }
        self.search_indexed_up_to = len;
    }

    /// Start the look-ahead search over the whole file. Skipped for stdin,
    /// which cannot be re-read; its loaded reads are indexed by
    /// [`Mate::index_loaded`].
    fn start_search(&mut self, pattern: &Pattern, scope: SearchScope) {
        self.search_worker.stop();
        if self.is_stdin() {
            self.search_done.store(true, Ordering::SeqCst);
            return;
        }
        self.search_done.store(false, Ordering::SeqCst);
        self.search_scanned.store(0, Ordering::SeqCst);

        let input = self.input.clone();
        let pattern = pattern.clone();
        let index = Arc::clone(&self.search_index);
        let scanned = Arc::clone(&self.search_scanned);
        let done = Arc::clone(&self.search_done);
        self.search_worker = Worker::spawn(move |stop| {
            scan_search(&input, &pattern, scope, &index, &scanned, &done, stop)
        });
    }

    /// True while any background scan of this file is still running.
    fn is_busy(&self) -> bool {
        self.counter.is_running()
            || self.stats_worker.is_running()
            || self.search_worker.is_running()
    }

    fn worker_error(&self) -> Option<String> {
        self.counter
            .error()
            .or_else(|| self.stats_worker.error())
            .or_else(|| self.search_worker.error())
    }
}

/// An input file and the reference needed to decode it (CRAM only), so
/// background scans can open their own reader on it.
#[derive(Clone)]
struct Input {
    path: String,
    reference: Option<String>,
}

impl Input {
    fn open(&self) -> Result<Box<dyn RecordSource>> {
        open_source(&self.path, self.reference.as_deref())
    }
}

/// Count the reads of `path`, publishing the running total.
fn count_reads(
    input: &Input,
    total: &AtomicU64,
    done: &AtomicBool,
    stop: &AtomicBool,
) -> Result<()> {
    let mut reader = input.open()?;
    let mut count = 0u64;
    while reader.next_record()?.is_some() {
        if stop.load(Ordering::Relaxed) {
            return Ok(());
        }
        count += 1;
        if count.is_multiple_of(COUNT_PROGRESS_EVERY) {
            total.store(count, Ordering::Relaxed);
        }
    }
    total.store(count, Ordering::Relaxed);
    done.store(true, Ordering::Relaxed);
    Ok(())
}

/// Compute statistics over the whole file, publishing a snapshot every
/// [`STATS_PROGRESS_EVERY`] reads and once at the end.
fn scan_stats(
    input: &Input,
    phred_range: PhredRange,
    stats: &Mutex<FastqStats>,
    stop: &AtomicBool,
) -> Result<()> {
    let mut reader = input.open()?;
    let mut acc = StatsAccumulator::new(phred_range);
    while let Some(record) = reader.next_record()? {
        if stop.load(Ordering::Relaxed) {
            return Ok(());
        }
        acc.add(record);
        if acc.reads().is_multiple_of(STATS_PROGRESS_EVERY) {
            acc.snapshot_into(&mut stats.lock().unwrap(), false);
        }
    }
    if !stop.load(Ordering::Relaxed) {
        acc.snapshot_into(&mut stats.lock().unwrap(), true);
    }
    Ok(())
}

/// Stream the file once, setting index bits for every record containing a
/// match. This provides look-ahead past the loaded buffer without storing
/// anything but the one-bit-per-record index.
fn scan_search(
    input: &Input,
    pattern: &Pattern,
    scope: SearchScope,
    index: &Mutex<SearchIndex>,
    scanned: &AtomicU64,
    done: &AtomicBool,
    stop: &AtomicBool,
) -> Result<()> {
    let mut reader = input.open()?;
    let mut record_idx = 0u64;
    let mut scratch = MatchScratch::new();
    while let Some(record) = reader.next_record()? {
        if stop.load(Ordering::Relaxed) {
            return Ok(());
        }
        if scope.matches_with(pattern, record, &mut scratch) {
            index.lock().unwrap().set(record_idx);
        }
        record_idx += 1;
        if record_idx.is_multiple_of(SEARCH_PROGRESS_EVERY) {
            scanned.store(record_idx, Ordering::Relaxed);
        }
    }
    scanned.store(record_idx, Ordering::SeqCst);
    done.store(true, Ordering::SeqCst);
    Ok(())
}

/// Calculate the layout for the stats page with blocks
fn calculate_stats_layout(area: Rect) -> (Rect, Vec<Vec<Rect>>) {
    let main_layout = Layout::vertical([Constraint::Length(1), Constraint::Min(0)]);
    let [title_area, main_area] = main_layout.areas(area);

    // Create a 2x3 grid for stats blocks
    let block_layout = Layout::vertical([
        Constraint::Length(12), // Basic stats and Adapter stats (side by side)
        Constraint::Length(6),  // Quality histogram and Average Quality histogram (side by side)
        Constraint::Min(6),     // Position quality chart
    ]);

    let main_areas = block_layout
        .split(main_area)
        .iter()
        .enumerate()
        .map(|(i, &area)| {
            match i {
                0 | 1 => {
                    // First and second rows - split horizontally for two charts
                    Layout::horizontal([Constraint::Percentage(50), Constraint::Percentage(50)])
                        .split(area)
                        .to_vec()
                }
                _ => {
                    // Third row - full width
                    Layout::horizontal([Constraint::Percentage(100)])
                        .split(area)
                        .to_vec()
                }
            }
        })
        .collect();
    (title_area, main_areas)
}

/// Placeholder for a quality panel when the input has no base qualities
/// (FASTA), so the chart is not mistaken for an all-zero distribution.
fn no_quality_message(title: &'static str) -> Paragraph<'static> {
    Paragraph::new(Line::from("No quality data in input (FASTA)"))
        .block(Block::default().title(title).borders(Borders::ALL))
        .alignment(Alignment::Center)
        .wrap(Wrap { trim: true })
}

/// Create quality histogram bar chart
fn create_quality_histogram<'a>(stats: &FastqStats, color_scheme: &ColorScheme) -> BarChart<'a> {
    let bin_width = 5;
    let bars: Vec<Bar> = (0..=40)
        .step_by(bin_width)
        .map(|q| {
            // Sum counts for this bin (e.g., Q0-4, Q5-9, etc.)
            let bin_count: u64 = (q..(q + bin_width).min(stats.quality_histogram.len()))
                .map(|idx| stats.quality_histogram.get(idx).unwrap_or(&0))
                .sum();

            Bar::default()
                .value(bin_count)
                .label(Line::from(format!("Q{}", q)))
                .style(quality_score_style(q as u8, color_scheme))
        })
        .collect();

    let title = Line::from(Span::styled(
        "Per-Base Quality Score Distribution",
        Style::default(),
    ));
    BarChart::default()
        .data(BarGroup::default().bars(&bars))
        .block(Block::new().title(title).borders(Borders::ALL))
        .bar_width(3)
        .bar_gap(1)
}
/// Create average read quality histogram bar chart
fn create_average_read_quality_histogram<'a>(
    stats: &FastqStats,
    color_scheme: &ColorScheme,
) -> BarChart<'a> {
    let bars: Vec<Bar> = stats
        .average_read_quality_histogram
        .iter()
        .enumerate()
        .map(|(i, &count)| {
            let range_start = i * AVG_QUAL_BIN_WIDTH;

            let label = if i == AVG_QUAL_NUM_BINS - 1 {
                format!("Q{:.0}+", range_start) // Last bin is "40+"
            } else {
                format!("Q{:.0}", range_start)
            };

            Bar::default()
                .value(count)
                .label(Line::from(label))
                .style(quality_score_style(range_start as u8, color_scheme))
        })
        .collect();

    let title = Line::from(Span::styled(
        "Average Quality Per Read Distribution",
        Style::default(),
    ));
    BarChart::default()
        .data(BarGroup::default().bars(&bars))
        .block(Block::new().title(title).borders(Borders::ALL))
        .bar_width(3)
        .bar_gap(1)
}

fn create_position_quality_chart<'a>(
    stats: &FastqStats,
    color_scheme: &ColorScheme,
) -> BarChart<'a> {
    let display_positions = stats.position_quality.len().min(50);
    let nbars = min(display_positions, 30);
    let bin_width = display_positions / nbars.max(1);
    // bars are averaged in bins so position 0-5, 6-10 etc.
    let bars = (0..display_positions)
        .step_by(bin_width) // Show every 2nd position to fit better
        .map(|pos| {
            // need to get all average quality scores for positions pos to pos + bin_width and average them
            let avg_qual = stats
                .position_quality
                .iter()
                .skip(pos)
                .take(bin_width)
                .sum::<f64>()
                / bin_width as f64;
            Bar::default()
                .value((avg_qual) as u64)
                .style(quality_score_style(avg_qual as u8, color_scheme))
        })
        .collect::<Vec<_>>();

    let title = Line::from(Span::styled(
        "Average Quality by Position",
        Style::default(),
    ));
    BarChart::default()
        .data(BarGroup::default().bars(&bars))
        .block(Block::new().title(title).borders(Borders::ALL))
        .bar_width(2)
        .bar_gap(1)
}

/// Quality character at `pos`, or the lowest score when the record has no
/// quality for that base.
fn quality_at(quality: &[u8], pos: usize) -> u8 {
    quality.get(pos).copied().unwrap_or(0)
}

/// Get color style based on quality score
fn quality_score_style(quality: u8, color_scheme: &ColorScheme) -> Style {
    Style::default().fg(color_scheme.quality_to_color(quality))
}

fn stdin_has_data() -> bool {
    let fd = stdin().as_raw_fd();
    let fd = unsafe { BorrowedFd::borrow_raw(fd) };
    let mut fds = [PollFd::new(fd, PollFlags::POLLIN)];
    match poll(&mut fds, nix::poll::PollTimeout::from(Some(250_u16))) {
        Ok(n) if n > 0 => fds[0].revents().unwrap().contains(PollFlags::POLLIN),
        _ => false,
    }
}

fn spawn_input_thread<I>(key_iter: I) -> Receiver<Key>
where
    I: Iterator<Item = std::io::Result<Key>> + Send + 'static,
{
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        for k in key_iter.flatten() {
            if tx.send(k).is_err() {
                break;
            }
        }
    });
    rx
}

/// Create adapter statistics display
fn create_adapter_stats_display(stats: &FastqStats) -> Paragraph<'_> {
    let adapter_stats = &stats.adapter_stats;

    let contamination_rate = if stats.processed_reads > 0 {
        (adapter_stats.contaminated_reads as f64 / stats.processed_reads as f64) * 100.0
    } else {
        0.0
    };

    let avg_adapters_per_read = if adapter_stats.contaminated_reads > 0 {
        adapter_stats.total_adapters_found as f64 / adapter_stats.contaminated_reads as f64
    } else {
        0.0
    };

    let mut content = vec![
        Line::from(format!(
            "Contaminated reads: {} ({:.1}%)",
            adapter_stats
                .contaminated_reads
                .to_formatted_string(&Locale::en),
            contamination_rate
        )),
        Line::from(format!(
            "Total adapter hits: {} (avg {:.1} types per contaminated read)",
            adapter_stats
                .total_adapters_found
                .to_formatted_string(&Locale::en),
            avg_adapters_per_read
        )),
        Line::from(""),
    ];

    // Show top 3 most common adapters
    let mut adapter_pairs: Vec<_> = adapter_stats.adapters_detected.iter().collect();
    adapter_pairs.sort_by_key(|(_, adapter_match)| std::cmp::Reverse(adapter_match.count));

    content.push(Line::from("Top adapters detected:"));

    if adapter_pairs.is_empty() {
        content.push(Line::from("  No adapters detected"));
    } else {
        for (i, (_, adapter_match)) in adapter_pairs.iter().take(3).enumerate() {
            content.push(Line::from(format!(
                "  {}. {}: {} reads (avg pos: {:.1})",
                i + 1,
                adapter_match.name,
                adapter_match.count.to_formatted_string(&Locale::en),
                adapter_match.avg_position
            )));
        }
    }

    let title = Line::from(Span::styled("Adapter Contamination", Style::default()));

    Paragraph::new(content)
        .block(Block::default().title(title).borders(Borders::ALL))
        .wrap(Wrap { trim: true })
}

impl TuiViewer {
    /// Calculate how many lines a record takes on screen
    fn calculate_record_lines(&self, record: &Record, terminal_width: usize) -> usize {
        calculate_record_lines(
            record,
            terminal_width,
            self.no_wrap,
            self.show_quality,
            self.show_translation,
        )
    }

    /// Calculate how many lines a paired record (two mates) takes on screen.
    /// For paired mode, each "record" is a pair: 1 shared header + 1 R1 seq
    /// (+ quality) + 1 R2 seq (+ quality), plus the gutter prefix.
    fn calculate_paired_record_lines(
        &self,
        r1: &Record,
        r2: &Record,
        terminal_width: usize,
    ) -> usize {
        let header_lines = 1;

        if self.show_translation {
            // 6 frames for R1 + 6 frames for R2, each with gutter
            let frames_r1 = translate_frames(r1.seq());
            let frames_r2 = translate_frames(r2.seq());
            let frame_lines: usize = frames_r1
                .iter()
                .chain(frames_r2.iter())
                .map(|f| {
                    if self.no_wrap {
                        1
                    } else {
                        f.len().max(1).div_ceil(terminal_width.saturating_sub(2))
                    }
                })
                .sum();
            return header_lines + frame_lines;
        }

        let seq_len_r1 = r1.seq().len();
        let seq_len_r2 = r2.seq().len();

        let seq_lines_r1 = if self.no_wrap {
            1
        } else {
            seq_len_r1.max(1).div_ceil(terminal_width.saturating_sub(2))
        };
        let seq_lines_r2 = if self.no_wrap {
            1
        } else {
            seq_len_r2.max(1).div_ceil(terminal_width.saturating_sub(2))
        };

        let quality_lines = if self.show_quality && r1.has_quality() {
            seq_lines_r1
        } else {
            0
        } + if self.show_quality && r2.has_quality() {
            seq_lines_r2
        } else {
            0
        };

        header_lines + seq_lines_r1 + seq_lines_r2 + quality_lines
    }

    /// Calculate how many records fit on screen starting from given position
    /// Returns (records_that_fit, total_lines_used)
    fn calculate_records_per_page(
        &self,
        start_position: u64,
        terminal_height: usize,
        terminal_width: usize,
    ) -> Result<(usize, usize)> {
        let available_height = terminal_height.saturating_sub(2); // Reserve space for header/footer

        if self.is_paired() {
            return self.calculate_paired_records_per_page(
                start_position,
                available_height,
                terminal_width,
            );
        }

        // Read from the already loaded buffer without triggering additional loads
        let reads_guard = self.mates[0].buffer.reads.read().unwrap();
        let start_idx = start_position as usize;

        if start_idx >= reads_guard.len() {
            return Ok((0, 0)); // No records available at this position
        }

        let mut lines_used = 0;
        let mut records_count = 0;

        // Check records from start_position onwards
        for record in reads_guard.iter().skip(start_idx) {
            let record_lines = self.calculate_record_lines(record, terminal_width);

            if lines_used + record_lines > available_height {
                break; // This record wouldn't fit
            }

            lines_used += record_lines;
            records_count += 1;
        }

        Ok((records_count, lines_used))
    }

    /// For paired mode: calculate how many pairs fit on screen
    fn calculate_paired_records_per_page(
        &self,
        start_position: u64,
        available_height: usize,
        terminal_width: usize,
    ) -> Result<(usize, usize)> {
        let reads_guard = self.mates[0].buffer.reads.read().unwrap();
        let reads_guard_r2 = self.mates[1].buffer.reads.read().unwrap();
        let start_idx = start_position as usize;

        if start_idx >= reads_guard.len() || start_idx >= reads_guard_r2.len() {
            return Ok((0, 0));
        }

        let mut lines_used = 0;
        let mut pairs_count = 0;

        for (r1, r2) in reads_guard
            .iter()
            .skip(start_idx)
            .zip(reads_guard_r2.iter().skip(start_idx))
        {
            let pair_lines = self.calculate_paired_record_lines(r1, r2, terminal_width);
            if lines_used + pair_lines > available_height {
                break;
            }
            lines_used += pair_lines;
            pairs_count += 1;
        }

        Ok((pairs_count, lines_used))
    }

    /// Calculate page down: find next position that fills the screen
    fn calculate_page_down(&self, terminal_height: usize, terminal_width: usize) -> Result<u64> {
        let (records_per_page, _) = self.calculate_records_per_page(
            self.current_position,
            terminal_height,
            terminal_width,
        )?;

        // Move by at least 1 record, or by the calculated page size
        let records_to_move = records_per_page.max(1);

        Ok(self.current_position + records_to_move as u64)
    }

    /// Calculate page up: find previous position that would fill the screen when moving forward
    fn calculate_page_up(&self, terminal_height: usize, terminal_width: usize) -> Result<u64> {
        if self.current_position == 0 {
            return Ok(0);
        }

        if self.is_paired() {
            return self.calculate_paired_page_up(terminal_height, terminal_width);
        }

        // We need to work backwards to find where to start so that moving forward
        // From current position we move bacjwards -1 and sum up the lines
        // needed until we reach the terminal height
        let mut lines_used = 0;
        let mut records_to_move = 0;
        let reads_guard = self.mates[0].buffer.reads.read().unwrap();
        let start_idx = self.current_position as usize;
        for record in reads_guard.iter().take(start_idx).rev() {
            let record_lines = self.calculate_record_lines(record, terminal_width);
            if lines_used + record_lines > terminal_height {
                break; // This record would not fit
            }
            lines_used += record_lines;
            records_to_move += 1;
        }
        // return the position that would fill the screen
        let new_position = if records_to_move > 0 {
            self.current_position.saturating_sub(records_to_move as u64)
        } else {
            0 // If no records fit, stay at the start
        };
        Ok(new_position)
    }

    /// Page up for paired mode: count backwards in pairs
    fn calculate_paired_page_up(
        &self,
        terminal_height: usize,
        terminal_width: usize,
    ) -> Result<u64> {
        if self.current_position == 0 {
            return Ok(0);
        }

        let mut lines_used = 0;
        let mut pairs_to_move = 0;
        let reads_guard = self.mates[0].buffer.reads.read().unwrap();
        let reads_guard_r2 = self.mates[1].buffer.reads.read().unwrap();
        let start_idx = self.current_position as usize;

        // Iterate backwards through pairs
        for i in (0..start_idx).rev() {
            if i >= reads_guard.len() || i >= reads_guard_r2.len() {
                break;
            }
            let pair_lines = self.calculate_paired_record_lines(
                &reads_guard[i],
                &reads_guard_r2[i],
                terminal_width,
            );
            if lines_used + pair_lines > terminal_height {
                break;
            }
            lines_used += pair_lines;
            pairs_to_move += 1;
        }

        let new_position = if pairs_to_move > 0 {
            self.current_position.saturating_sub(pairs_to_move as u64)
        } else {
            0
        };
        Ok(new_position)
    }

    pub fn new(
        file_path: String,
        file_path_opt2: Option<String>,
        reference: Option<String>,
    ) -> Result<Self> {
        let running = Arc::new(AtomicBool::new(true));
        flag::register(SIGINT, Arc::clone(&running))?;

        if file_path == "-" && !stdin_has_data() {
            return Err(anyhow::anyhow!(
                "No data available on stdin. Maybe you forgot to pipe input?"
            ));
        }

        let stdout = stdout().into_raw_mode()?;
        let stdout = stdout.into_alternate_screen()?;
        let backend = TermionBackend::new(stdout);
        let terminal = Terminal::new(backend)?;

        let mut mates = vec![Mate::open(file_path, reference.clone())?];
        if let Some(path) = file_path_opt2 {
            mates.push(Mate::open(path, reference)?);
        }

        let mut viewer = TuiViewer {
            terminal,
            mates,
            current_position: 0,
            horizontal_offset: 0,
            no_wrap: false,
            show_quality: false,
            show_stats: false,
            show_help: false,
            orientation: ReadOrientation::AsIs,
            show_translation: false,
            stats_scroll: 0,
            stats_tab: 0,
            help_scroll: 0,
            phred_range: PhredRange::Default,
            color_scheme: ColorScheme::RedGreen,
            search_mode: false,
            search_input: String::new(),
            search_query: String::new(),
            search_pattern: None,
            search_include_rc: true,
            search_scope: SearchScope::default(),
            search_scope_auto: false,
            search_match_k: 0,
        };

        // Load the first window so the phred range can be detected before
        // any statistics are computed with it.
        for mate in &mut viewer.mates {
            mate.buffer.load_window(0, BUFFER_WINDOW_SIZE)?;
        }
        // Inputs without qualities (FASTA) keep the default range.
        if let Some((phred_min, phred_max)) =
            determine_min_max_phred(&viewer.mates[0].buffer.reads.read().unwrap())
        {
            viewer.phred_range = PhredRange::from_min_max(phred_min, phred_max);
        }

        for mate in &mut viewer.mates {
            mate.start_counter();
        }
        viewer.restart_stats();

        Ok(viewer)
    }

    fn is_paired(&self) -> bool {
        self.mates.len() == 2
    }

    /// (Re)compute the statistics of every mate, e.g. after the phred range
    /// changed.
    fn restart_stats(&mut self) {
        for mate in &mut self.mates {
            mate.start_stats(self.phred_range);
        }
    }

    /// Handle a key while the search prompt is active.
    fn handle_search_key(&mut self, key: Key) {
        match key {
            Key::Char('\n') => {
                self.search_input = self.search_input.trim().to_string();
                self.search_mode = false;
                self.search_query = self.search_input.clone();
                self.apply_search_pattern();
                // Like less: jump to the first match at or after the current
                // position once the pattern is committed.
                self.jump_to_first_match();
            }
            Key::Char(c) => {
                if !c.is_control() {
                    self.search_input.push(c);
                }
            }
            Key::Backspace => {
                self.search_input.pop();
            }
            Key::Esc | Key::Ctrl('c') => {
                self.search_mode = false;
                self.search_input.clear();
            }
            _ => {}
        }
    }

    /// Build the search pattern for the current view: a nucleotide pattern
    /// (with optional reverse complement) in the sequence view, a peptide
    /// pattern in the 6-frame translation view.
    fn build_search_pattern(&self) -> Option<Pattern> {
        if self.show_translation {
            Pattern::new_peptide(&self.search_query)
        } else {
            Pattern::new(&self.search_query, self.search_include_rc)
        }
    }

    /// (Re)build the active search pattern from the committed query, reset
    /// the index, index already-loaded records and start the look-ahead worker.
    fn apply_search_pattern(&mut self) {
        if self.search_query.trim().is_empty() {
            self.clear_search();
            return;
        }
        let pattern = match self.build_search_pattern() {
            Some(p) => p,
            None => {
                self.clear_search();
                return;
            }
        };
        // Name-like queries contain characters that are not nucleotide codes.
        // Search the read name automatically instead of silently finding no
        // sequence matches. A manual scope choice is left alone.
        if !self.show_translation {
            if self.search_scope_auto {
                self.search_scope = SearchScope::Sequence;
                self.search_scope_auto = false;
            }
            if self.search_scope == SearchScope::Sequence && !pattern.is_nucleotide_literal() {
                self.search_scope = SearchScope::Id;
                self.search_scope_auto = true;
            }
        }
        self.search_match_k = 0;
        for mate in &mut self.mates {
            mate.reset_search();
        }
        self.search_pattern = Some(pattern.clone());
        self.index_loaded_records();
        for mate in &mut self.mates {
            mate.start_search(&pattern, self.search_scope);
        }
    }

    /// Clear the active search and stop the look-ahead worker.
    fn clear_search(&mut self) {
        self.search_pattern = None;
        self.search_query.clear();
        // Reset the per-search reverse-complement toggle, but keep the scope:
        // it is a sticky preference for how the user likes to search.
        self.search_include_rc = true;
        self.search_match_k = 0;
        for mate in &mut self.mates {
            mate.reset_search();
        }
    }

    /// Scan newly loaded records on the main thread and mark matches in the
    /// index. Cheap per frame: only the delta since the last scan is visited.
    fn index_loaded_records(&mut self) {
        let pattern = self.search_pattern.as_ref();
        for mate in &mut self.mates {
            mate.index_loaded(pattern, self.search_scope);
        }
    }

    /// Jump to a record index, streaming forward if it is not loaded yet.
    fn jump_to_record(&mut self, record: u64) -> Result<()> {
        self.current_position = record;
        for mate in &mut self.mates {
            mate.buffer.load_window(record, BUFFER_WINDOW_SIZE)?;
        }
        let len = self.loaded_len();
        if len > 0 && self.current_position >= len {
            self.current_position = len - 1;
        }
        Ok(())
    }

    /// Number of records that can be shown. In paired mode the view zips R1
    /// and R2, so only positions loaded in both files count.
    fn loaded_len(&self) -> u64 {
        self.mates.iter().map(Mate::reads_loaded).min().unwrap_or(0)
    }

    /// Move forward to `target` if that record exists; otherwise stay put, so
    /// scrolling or paging past the last record keeps the last page on screen
    /// instead of an empty one (and an empty input cannot underflow).
    fn move_forward_to(&mut self, target: u64) -> Result<()> {
        for mate in &mut self.mates {
            mate.buffer.load_window(target, BUFFER_WINDOW_SIZE)?;
        }
        if target < self.loaded_len() {
            self.current_position = target;
        }
        Ok(())
    }

    /// Find a match relative to `from`. In paired mode a match in either mate
    /// counts. Returns the record and its 1-based number among all matches.
    fn find_match(&self, from: u64, direction: MatchDirection) -> Option<(u64, u64)> {
        let index = self.mates[0].search_index.lock().unwrap();
        if let Some(r2) = self.mates.get(1) {
            let index_r2 = r2.search_index.lock().unwrap();
            let target = match direction {
                MatchDirection::AtOrAfter => index.next_union(&index_r2, from),
                MatchDirection::Before => index.prev_union(&index_r2, from),
            }?;
            Some((target, index.count_union_up_to(&index_r2, target)))
        } else {
            let target = match direction {
                MatchDirection::AtOrAfter => index.next(from),
                MatchDirection::Before => index.prev(from),
            }?;
            Some((target, index.count_up_to(target)))
        }
    }

    /// Total number of matching records (pairs count once).
    fn match_count(&self) -> u64 {
        let index = self.mates[0].search_index.lock().unwrap();
        match self.mates.get(1) {
            Some(r2) => index.count_union_up_to(&r2.search_index.lock().unwrap(), u64::MAX),
            None => index.count(),
        }
    }

    /// Jump to the match at or after the current position, including the
    /// current record. Used when a search is committed or its scope changes.
    fn jump_to_first_match(&mut self) {
        if let Some((target, k)) = self.find_match(self.current_position, MatchDirection::AtOrAfter)
        {
            self.search_match_k = k;
            let _ = self.jump_to_record(target);
        }
    }

    /// Move to the next (or previous) matching record relative to the current
    /// position. In paired mode, a match in either R1 or R2 counts.
    fn next_match(&mut self, forward: bool) -> Result<()> {
        if self.search_pattern.is_none() {
            return Ok(());
        }
        let found = if forward {
            self.find_match(self.current_position + 1, MatchDirection::AtOrAfter)
        } else {
            self.find_match(self.current_position, MatchDirection::Before)
        };
        if let Some((target, k)) = found {
            self.search_match_k = k;
            self.jump_to_record(target)?;
        }
        Ok(())
    }

    /// The first error any background worker reported, for the status line.
    fn worker_error(&self) -> Option<String> {
        let paired = self.is_paired();
        self.mates.iter().enumerate().find_map(|(i, mate)| {
            let error = mate.worker_error()?;
            Some(if paired {
                format!("R{}: {}", i + 1, error)
            } else {
                error
            })
        })
    }

    /// The pattern to use for highlighting right now: while typing a new
    /// pattern in the search prompt it updates live, otherwise the committed
    /// pattern applies. The pattern is interpreted for the current view
    /// (nucleotide vs. amino acid).
    fn highlight_pattern(&self) -> Option<Pattern> {
        if self.search_mode && !self.search_input.trim().is_empty() {
            if self.show_translation {
                Pattern::new_peptide(self.search_input.trim())
            } else {
                Pattern::new(self.search_input.trim(), self.search_include_rc)
            }
        } else {
            self.search_pattern.clone()
        }
    }

    /// Build a record header line. When the active search scope includes read
    /// names, matching parts of the ID and description are highlighted.
    fn header_line<'a>(&self, id: &str, desc: Option<&str>) -> Line<'a> {
        let mut spans: Vec<Span<'a>> = Vec::new();
        let pattern = if self.search_scope.includes_id() {
            self.highlight_pattern()
        } else {
            None
        };

        match pattern {
            Some(pattern) => {
                spans.extend(Self::highlighted_spans(
                    id,
                    &pattern.find_id_matches(id.as_bytes()),
                ));
                if let Some(desc) = desc {
                    spans.push(Span::raw(" "));
                    spans.extend(Self::highlighted_spans(
                        desc,
                        &pattern.find_id_matches(desc.as_bytes()),
                    ));
                }
            }
            None => {
                spans.push(Span::raw(id.to_string()));
                if let Some(desc) = desc {
                    spans.push(Span::raw(format!(" {}", desc)));
                }
            }
        }

        Line::from(spans)
    }

    /// Split `text` into plain and yellow-highlighted spans using byte ranges.
    fn highlighted_spans<'a>(text: &str, ranges: &[(usize, usize)]) -> Vec<Span<'a>> {
        if ranges.is_empty() {
            return vec![Span::raw(text.to_string())];
        }

        let mut spans = Vec::new();
        let mut pos = 0;
        for &(start, end) in ranges {
            if start > pos {
                spans.push(Span::raw(text[pos..start].to_string()));
            }
            spans.push(Span::styled(
                text[start..end].to_string(),
                Style::default().fg(Color::Black).bg(Color::Yellow),
            ));
            pos = end;
        }
        if pos < text.len() {
            spans.push(Span::raw(text[pos..].to_string()));
        }
        spans
    }

    /// Render a base as either its quality color or the search-match
    /// highlight (yellow background) when it falls inside a match range.
    fn base_span<'a>(&self, base: u8, quality: u8, highlighted: bool) -> Span<'a> {
        if highlighted {
            Span::styled(
                char::from(base).to_string(),
                Style::default().fg(Color::Black).bg(Color::Yellow),
            )
        } else {
            let quality_score = quality.saturating_sub(self.phred_range.base_phred());
            let color = self.color_scheme.quality_to_color(quality_score);
            Span::styled(char::from(base).to_string(), Style::default().fg(color))
        }
    }

    /// Render an amino acid as either its quality color (stops in red) or the
    /// search-match highlight (yellow background) when inside a match range.
    /// Quality scores here are absolute (0-40), so no phred offset is applied.
    fn aa_span<'a>(&self, ch: u8, quality: Option<u8>, highlighted: bool) -> Span<'a> {
        if highlighted {
            Span::styled(
                char::from(ch).to_string(),
                Style::default().fg(Color::Black).bg(Color::Yellow),
            )
        } else if ch == b'*' {
            Span::styled(char::from(ch).to_string(), Style::default().fg(Color::Red))
        } else {
            let style = match quality {
                Some(quality) => Style::default().fg(self.color_scheme.quality_to_color(quality)),
                // Without base qualities (FASTA) the frame is not shaded.
                None => Style::default(),
            };
            Span::styled(char::from(ch).to_string(), style)
        }
    }

    /// Colour bases by quality. A malformed record whose quality string is
    /// shorter than its sequence is shown rather than crashing the viewer:
    /// bases without a quality are drawn as Q0.
    fn colorize_sequence<'a>(&self, sequence: &[u8], quality: &[u8]) -> Vec<Span<'a>> {
        let ranges = self
            .highlight_pattern()
            .map(|p| p.find_merged_matches(sequence))
            .unwrap_or_default();
        let mut spans = Vec::with_capacity(sequence.len());
        let mut pos = 0;
        for (start, end) in &ranges {
            for (i, &base) in sequence[pos..*start].iter().enumerate() {
                spans.push(self.base_span(base, quality_at(quality, pos + i), false));
            }
            for &base in &sequence[*start..*end] {
                spans.push(self.base_span(base, 0, true));
            }
            pos = *end;
        }
        for (i, &base) in sequence[pos..].iter().enumerate() {
            spans.push(self.base_span(base, quality_at(quality, pos + i), false));
        }
        spans
    }

    /// Plain-text sequence line (shown when quality is displayed separately),
    /// with search matches highlighted.
    fn raw_sequence_with_highlight<'a>(&self, sequence: &[u8]) -> Vec<Span<'a>> {
        let ranges = self
            .highlight_pattern()
            .map(|p| p.find_merged_matches(sequence))
            .unwrap_or_default();
        let mut spans = Vec::with_capacity(sequence.len());
        let mut pos = 0;
        for (start, end) in &ranges {
            for &base in &sequence[pos..*start] {
                spans.push(Span::raw(char::from(base).to_string()));
            }
            for &base in &sequence[*start..*end] {
                spans.push(self.base_span(base, 0, true));
            }
            pos = *end;
        }
        for &base in &sequence[pos..] {
            spans.push(Span::raw(char::from(base).to_string()));
        }
        spans
    }

    /// Six-frame translation of `seq`, with per-amino-acid quality scores when
    /// the record has qualities and zeros otherwise (rendered unshaded).
    fn translation_frames(&self, seq: &[u8], qual: Option<&[u8]>) -> Vec<(String, Vec<u8>)> {
        match qual {
            Some(qual) => translate_frames_with_quality(seq, qual, self.phred_range.base_phred()),
            None => translate_frames(seq)
                .into_iter()
                .map(|frame| {
                    let len = frame.len();
                    (frame, vec![0u8; len])
                })
                .collect(),
        }
    }

    /// Return the sequence and quality of a record in the currently selected
    /// orientation. Quality stays paired with its base, so it is reversed
    /// together with the sequence for the reverse complement. Records without
    /// qualities (FASTA) yield an empty quality vector.
    fn oriented_seq_qual(&self, record: &Record) -> (Vec<u8>, Vec<u8>) {
        match self.orientation {
            ReadOrientation::AsIs => (
                record.seq().to_vec(),
                record.qual().unwrap_or_default().to_vec(),
            ),
            ReadOrientation::ReverseComplement => (
                dna::revcomp(record.seq()),
                record
                    .qual()
                    .unwrap_or_default()
                    .iter()
                    .rev()
                    .copied()
                    .collect(),
            ),
        }
    }

    pub fn run(&mut self) -> Result<()> {
        // Determine where to read keyboard input from:
        // - If reading from a file, stdin is available for keyboard input
        // - If reading from stdin (piped), we need to use /dev/tty for keyboard input
        let use_tty = self.mates.iter().any(Mate::is_stdin);
        let key_source: Box<dyn Iterator<Item = std::io::Result<Key>> + Send> = if use_tty {
            Box::new(File::open("/dev/tty")?.keys())
        } else {
            Box::new(stdin().keys())
        };

        let rx = spawn_input_thread(key_source);

        self.terminal.clear()?;
        let mut needs_redraw = true;
        let mut was_busy = true;
        let mut last_terminal_size = self.terminal.size()?;
        loop {
            // Check for terminal resize
            let current_size = self.terminal.size()?;
            if current_size != last_terminal_size {
                needs_redraw = true;
                last_terminal_size = current_size;
            }

            // Only draw if something changed
            if needs_redraw {
                self.draw()?;
                needs_redraw = false;
            }

            // Read keyboard input - from stdin for files, /dev/tty for piped input
            let mut had_input = false;
            while let Ok(key) = rx.try_recv() {
                had_input = true;
                needs_redraw = true;
                if self.search_mode {
                    self.handle_search_key(key);
                    continue;
                }
                match key {
                    Key::Char('q') | Key::Ctrl('c') => return Ok(()),
                    Key::Char('/') => {
                        if !self.show_stats && !self.show_help {
                            self.search_mode = true;
                            self.search_input.clear();
                        }
                    }
                    Key::Char('n') => {
                        if !self.show_stats && !self.show_help {
                            self.next_match(true)?;
                        }
                    }
                    Key::Char('N') => {
                        if !self.show_stats && !self.show_help {
                            self.next_match(false)?;
                        }
                    }
                    Key::Char('x') => {
                        if self.search_pattern.is_some() && !self.show_translation {
                            self.search_include_rc = !self.search_include_rc;
                            self.apply_search_pattern();
                        }
                    }
                    Key::Char('i') => {
                        if self.search_pattern.is_some() {
                            self.search_scope = self.search_scope.next();
                            self.search_scope_auto = false;
                            self.apply_search_pattern();
                            self.jump_to_first_match();
                        }
                    }
                    Key::Char('c') => {
                        // iterate color schemes
                        self.color_scheme = self.color_scheme.next();
                    }
                    Key::Char('s') => {
                        self.show_stats = !self.show_stats;
                        self.stats_scroll = 0; // Reset scroll when toggling
                    }
                    Key::Char('\t') => {
                        if self.show_stats && self.is_paired() {
                            self.stats_tab = (self.stats_tab + 1) % 2;
                            self.stats_scroll = 0;
                        }
                    }
                    Key::Char('p') => {
                        self.show_quality = !self.show_quality;
                    }
                    Key::Char('r') => {
                        // Cycle read orientation: 5'->3', reverse complement
                        self.orientation = self.orientation.next();
                        self.horizontal_offset = 0;
                    }
                    Key::Char('t') => {
                        self.show_translation = !self.show_translation;
                        self.horizontal_offset = 0;
                        // Reinterpret the active query for the new view:
                        // nucleotide pattern when leaving, peptide pattern
                        // when entering the 6-frame translation view.
                        if self.search_pattern.is_some() {
                            self.apply_search_pattern();
                        }
                    }
                    Key::Char('h') => {
                        self.show_help = !self.show_help;
                        self.help_scroll = 0; // Reset scroll when toggling
                    }
                    Key::Char('S') => {
                        if !self.show_stats && !self.show_help {
                            self.no_wrap = !self.no_wrap;
                            self.horizontal_offset = 0;
                        }
                    }
                    Key::Char('e') => {
                        // Adjust phred encoding range by iterating through known ranges
                        self.phred_range = match self.phred_range {
                            PhredRange::Solexa => PhredRange::Illumina1_3,
                            PhredRange::Illumina1_3 => PhredRange::Illumina1_5,
                            PhredRange::Illumina1_5 => PhredRange::Default,
                            PhredRange::Default => PhredRange::Solexa,
                        };
                        // Recompute the statistics with the new offset
                        self.restart_stats();
                    }
                    Key::Left => {
                        if self.no_wrap
                            && self.horizontal_offset > 0
                            && !self.show_stats
                            && !self.show_help
                        {
                            self.horizontal_offset = self.horizontal_offset.saturating_sub(10);
                        }
                    }
                    Key::Right => {
                        if self.no_wrap && !self.show_stats && !self.show_help {
                            self.horizontal_offset += 10;
                        }
                    }
                    Key::Down | Key::Char('j') => {
                        if self.show_stats {
                            self.stats_scroll += 1;
                        } else if self.show_help {
                            self.help_scroll += 1;
                        } else {
                            self.move_forward_to(self.current_position + 1)?;
                        }
                    }
                    Key::Up | Key::Char('k') => {
                        if self.show_stats {
                            self.stats_scroll = self.stats_scroll.saturating_sub(1);
                        } else if self.show_help {
                            self.help_scroll = self.help_scroll.saturating_sub(1);
                        } else if self.current_position > 0 {
                            self.current_position -= 1;
                        }
                    }
                    Key::PageDown | Key::Char('J') | Key::Char(' ') => {
                        if self.show_stats {
                            self.stats_scroll += PAGE_SCROLL_AMOUNT;
                        } else if self.show_help {
                            self.help_scroll += PAGE_SCROLL_AMOUNT;
                        } else {
                            let terminal_size = self.terminal.size()?;
                            let new_position = self.calculate_page_down(
                                terminal_size.height as usize,
                                terminal_size.width as usize,
                            )?;
                            self.move_forward_to(new_position)?;
                        }
                    }
                    Key::PageUp | Key::Char('K') => {
                        if self.show_stats {
                            self.stats_scroll =
                                self.stats_scroll.saturating_sub(PAGE_SCROLL_AMOUNT);
                        } else if self.show_help {
                            self.help_scroll = self.help_scroll.saturating_sub(PAGE_SCROLL_AMOUNT);
                        } else {
                            let terminal_size = self.terminal.size()?;
                            let new_position = self.calculate_page_up(
                                terminal_size.height as usize,
                                terminal_size.width as usize,
                            )?;
                            self.current_position = new_position;
                        }
                    }
                    Key::Home => {
                        if self.show_stats {
                            self.stats_scroll = 0;
                        } else if self.show_help {
                            self.help_scroll = 0;
                        } else {
                            self.current_position = 0;
                        }
                    }
                    _ => {}
                }
            }

            // Stdin statistics follow the reads loaded so far.
            for mate in &mut self.mates {
                mate.update_stdin_stats();
            }

            // Keep redrawing while background scans make progress (counts,
            // stats, search look-ahead), plus one frame after they finish.
            let busy = self.mates.iter().any(Mate::is_busy);
            if !had_input && (busy || was_busy) {
                needs_redraw = true;
            }
            was_busy = busy;

            // Small delay to prevent excessive CPU usage
            std::thread::sleep(Duration::from_millis(if needs_redraw { 100 } else { 250 }));
        }
    }

    fn draw(&mut self) -> Result<()> {
        // Keep the match index up to date for any records loaded since the
        // last frame (also covers the initial window).
        self.index_loaded_records();

        let terminal_size = self.terminal.size()?;
        let terminal_width = terminal_size.width as usize;

        let no_wrap = self.no_wrap;
        let show_stats = self.show_stats;
        let show_help = self.show_help;
        let horizontal_offset = self.horizontal_offset;
        let current_position = self.current_position;
        let help_scroll = self.help_scroll;

        // Prepare content for main area
        let mut prepared_lines = Vec::new();

        let available_height = terminal_size.height.saturating_sub(2) as usize;
        let max_visible = available_height / 2;

        if self.is_paired() {
            // --- PAIRED MODE: interleaved R1/R2 with gutter bar ---
            let records_r1 = self.mates[0]
                .buffer
                .get_window(self.current_position, max_visible + 5)?;
            let records_r2 = self.mates[1]
                .buffer
                .get_window(self.current_position, max_visible + 5)?;
            for (r1, r2) in records_r1.iter().zip(records_r2.iter()) {
                prepared_lines.push(self.header_line(r1.id(), r1.desc()));

                let (seq1, qual1) = self.oriented_seq_qual(r1);
                let (seq2, qual2) = self.oriented_seq_qual(r2);
                let has_qual1 = r1.has_quality();
                let has_qual2 = r2.has_quality();

                if self.show_translation {
                    let frames1 =
                        self.translation_frames(&seq1, has_qual1.then_some(qual1.as_slice()));
                    let frames2 =
                        self.translation_frames(&seq2, has_qual2.then_some(qual2.as_slice()));
                    let highlight = self.highlight_pattern();
                    let labels: Vec<String> = (0..6)
                        .map(|i| {
                            if i < 3 {
                                format!("F{}: ", i + 1)
                            } else {
                                format!("F{} (rc): ", i + 1)
                            }
                        })
                        .collect();
                    let label_width = labels.iter().map(|l| l.chars().count()).max().unwrap_or(0);

                    for (i, (frame, scores)) in frames1.iter().enumerate() {
                        let label = format!("{:<width$}1:\u{00a0}", labels[i], width = label_width);
                        let (vs, ve) = if no_wrap {
                            let s = horizontal_offset.min(frame.len());
                            (s, (s + terminal_width.saturating_sub(4)).min(frame.len()))
                        } else {
                            (0, frame.len())
                        };
                        let visible = &frame.as_bytes()[vs..ve];
                        let v_scores = &scores[vs..ve];
                        let ranges = highlight
                            .as_ref()
                            .map(|p| p.find_peptide_matches(frame.as_bytes()))
                            .unwrap_or_default();
                        let mut spans =
                            vec![Span::styled(label, Style::default().fg(Color::DarkGray))];
                        let mut pos = 0;
                        for (start, end) in &ranges {
                            let s = (*start).max(vs);
                            let e = (*end).min(ve);
                            if s >= e {
                                continue;
                            }
                            for (k, &ch) in visible[pos..s - vs].iter().enumerate() {
                                spans.push(self.aa_span(
                                    ch,
                                    has_qual1.then_some(v_scores[pos + k]),
                                    false,
                                ));
                            }
                            for &ch in &visible[s - vs..e - vs] {
                                spans.push(self.aa_span(ch, None, true));
                            }
                            pos = e - vs;
                        }
                        for (k, &ch) in visible[pos..].iter().enumerate() {
                            spans.push(self.aa_span(
                                ch,
                                has_qual1.then_some(v_scores[pos + k]),
                                false,
                            ));
                        }
                        prepared_lines.push(Line::from(spans));
                    }
                    for (i, (frame, scores)) in frames2.iter().enumerate() {
                        let label = format!("{:<width$}2:\u{00a0}", labels[i], width = label_width);
                        let (vs, ve) = if no_wrap {
                            let s = horizontal_offset.min(frame.len());
                            (s, (s + terminal_width.saturating_sub(4)).min(frame.len()))
                        } else {
                            (0, frame.len())
                        };
                        let visible = &frame.as_bytes()[vs..ve];
                        let v_scores = &scores[vs..ve];
                        let ranges = highlight
                            .as_ref()
                            .map(|p| p.find_peptide_matches(frame.as_bytes()))
                            .unwrap_or_default();
                        let mut spans =
                            vec![Span::styled(label, Style::default().fg(Color::DarkGray))];
                        let mut pos = 0;
                        for (start, end) in &ranges {
                            let s = (*start).max(vs);
                            let e = (*end).min(ve);
                            if s >= e {
                                continue;
                            }
                            for (k, &ch) in visible[pos..s - vs].iter().enumerate() {
                                spans.push(self.aa_span(
                                    ch,
                                    has_qual2.then_some(v_scores[pos + k]),
                                    false,
                                ));
                            }
                            for &ch in &visible[s - vs..e - vs] {
                                spans.push(self.aa_span(ch, None, true));
                            }
                            pos = e - vs;
                        }
                        for (k, &ch) in visible[pos..].iter().enumerate() {
                            spans.push(self.aa_span(
                                ch,
                                has_qual2.then_some(v_scores[pos + k]),
                                false,
                            ));
                        }
                        prepared_lines.push(Line::from(spans));
                    }
                    continue;
                }

                // R1 sequence
                let (vs1, vq1) = if no_wrap {
                    let s = horizontal_offset.min(seq1.len());
                    let e = (s + terminal_width.saturating_sub(3)).min(seq1.len());
                    (&seq1[s..e], &qual1[s.min(qual1.len())..e.min(qual1.len())])
                } else {
                    (&seq1[..], &qual1[..])
                };
                let prefix_r1 = Span::styled("1:\u{00a0}", Style::default());
                let mut r1s = vec![prefix_r1];
                r1s.extend(if !self.show_quality && has_qual1 {
                    self.colorize_sequence(vs1, vq1)
                } else {
                    self.raw_sequence_with_highlight(vs1)
                });
                prepared_lines.push(Line::from(r1s));
                if self.show_quality && has_qual1 {
                    let mut r1q = vec![Span::styled("1:\u{00a0}", Style::default())];
                    r1q.extend(self.colorize_sequence(vq1, vq1));
                    prepared_lines.push(Line::from(r1q));
                }

                // R2 sequence
                let (vs2, vq2) = if no_wrap {
                    let s = horizontal_offset.min(seq2.len());
                    let e = (s + terminal_width.saturating_sub(3)).min(seq2.len());
                    (&seq2[s..e], &qual2[s.min(qual2.len())..e.min(qual2.len())])
                } else {
                    (&seq2[..], &qual2[..])
                };
                let mut r2s = vec![Span::styled("2:\u{00a0}", Style::default())];
                r2s.extend(if !self.show_quality && has_qual2 {
                    self.colorize_sequence(vs2, vq2)
                } else {
                    self.raw_sequence_with_highlight(vs2)
                });
                prepared_lines.push(Line::from(r2s));
                if self.show_quality && has_qual2 {
                    let mut r2q = vec![Span::styled("2:\u{00a0}", Style::default())];
                    r2q.extend(self.colorize_sequence(vq2, vq2));
                    prepared_lines.push(Line::from(r2q));
                }
            }
        } else {
            // --- SINGLE FILE MODE ---
            let records = self.mates[0]
                .buffer
                .get_window(self.current_position, max_visible + 5)?;
            for record in records.iter() {
                prepared_lines.push(self.header_line(record.id(), record.desc()));
                let (oriented_seq, oriented_qual) = self.oriented_seq_qual(record);
                let has_qual = record.has_quality();

                if self.show_translation {
                    let frames = self.translation_frames(
                        &oriented_seq,
                        has_qual.then_some(oriented_qual.as_slice()),
                    );
                    let highlight = self.highlight_pattern();
                    let labels: Vec<String> = frames
                        .iter()
                        .enumerate()
                        .map(|(i, _)| {
                            if i < 3 {
                                format!("  F{}: ", i + 1)
                            } else {
                                format!("  F{} (rc): ", i + 1)
                            }
                        })
                        .collect();
                    let lw = labels.iter().map(|l| l.chars().count()).max().unwrap_or(0);
                    for (i, (frame, scores)) in frames.iter().enumerate() {
                        let label = format!("{:<width$}", labels[i], width = lw);
                        let (vs, ve) = if no_wrap {
                            let s = horizontal_offset.min(frame.len());
                            (s, (s + terminal_width.saturating_sub(1)).min(frame.len()))
                        } else {
                            (0, frame.len())
                        };
                        let visible = &frame.as_bytes()[vs..ve];
                        let v_scores = &scores[vs..ve];
                        let ranges = highlight
                            .as_ref()
                            .map(|p| p.find_peptide_matches(frame.as_bytes()))
                            .unwrap_or_default();
                        let mut spans =
                            vec![Span::styled(label, Style::default().fg(Color::DarkGray))];
                        let mut pos = 0;
                        for (start, end) in &ranges {
                            let s = (*start).max(vs);
                            let e = (*end).min(ve);
                            if s >= e {
                                continue;
                            }
                            for (k, &ch) in visible[pos..s - vs].iter().enumerate() {
                                spans.push(self.aa_span(
                                    ch,
                                    has_qual.then_some(v_scores[pos + k]),
                                    false,
                                ));
                            }
                            for &ch in &visible[s - vs..e - vs] {
                                spans.push(self.aa_span(ch, None, true));
                            }
                            pos = e - vs;
                        }
                        for (k, &ch) in visible[pos..].iter().enumerate() {
                            spans.push(self.aa_span(
                                ch,
                                has_qual.then_some(v_scores[pos + k]),
                                false,
                            ));
                        }
                        prepared_lines.push(Line::from(spans));
                    }
                    continue;
                }

                let (vs, vq) = if no_wrap {
                    let s = horizontal_offset.min(oriented_seq.len());
                    let e = (s + terminal_width.saturating_sub(1)).min(oriented_seq.len());
                    (
                        &oriented_seq[s..e],
                        &oriented_qual[s.min(oriented_qual.len())..e.min(oriented_qual.len())],
                    )
                } else {
                    (&oriented_seq[..], &oriented_qual[..])
                };
                let ss = if !self.show_quality && has_qual {
                    self.colorize_sequence(vs, vq)
                } else {
                    self.raw_sequence_with_highlight(vs)
                };
                prepared_lines.push(Line::from(ss));
                if self.show_quality && has_qual {
                    prepared_lines.push(Line::from(self.colorize_sequence(vq, vq)));
                }
            }
        }

        let wrap_status = if no_wrap { "NO-WRAP" } else { "WRAP" };
        let mode_display = if self.is_paired() {
            "PAIRED"
        } else {
            wrap_status
        };

        let help_text = if no_wrap {
            "↑/k: Up | ↓/j: Down | ←/→: Scroll | PgUp/PgDn: Page | S: Wrap | /: Search | s: Stats | c: Colors | q: Quit"
        } else {
            "↑/k: Up | ↓/j: Down | PgUp/PgDn: Page | S: No-Wrap | /: Search | s: Stats | c: Colors | q: Quit"
        };

        let search_status = if self.search_pattern.is_some() {
            let count = self.match_count();
            let mode = match &self.search_pattern {
                Some(p) if p.is_peptide() => "aa",
                Some(_) => {
                    if self.search_include_rc {
                        "+RC"
                    } else {
                        "-RC"
                    }
                }
                None => "",
            };
            let scope_prefix = match self.search_scope {
                SearchScope::Sequence => "",
                SearchScope::Id => "name ",
                SearchScope::SequenceAndId => "seq+name ",
            };
            // Look-ahead is complete once every mate has been scanned.
            let done = self
                .mates
                .iter()
                .all(|m| m.search_done.load(Ordering::SeqCst));
            let scanning = if done {
                String::new()
            } else {
                let scanned = self
                    .mates
                    .iter()
                    .map(|m| m.search_scanned.load(Ordering::SeqCst))
                    .min()
                    .unwrap_or(0);
                format!(" scanning…({})", scanned.to_formatted_string(&Locale::en))
            };
            format!(
                " | /{} {}{} {}/{}{}",
                self.search_query, scope_prefix, mode, self.search_match_k, count, scanning
            )
        } else {
            String::new()
        };

        let is_paired = self.is_paired();
        // Stats page tab: R1, or R2 in paired mode.
        let stats_mate = if is_paired { self.stats_tab.min(1) } else { 0 };
        let (stats_total, stats_total_exact) = self.mates[stats_mate].total_reads();
        let stats_from_stdin = self.mates[stats_mate].is_stdin();

        let file_label = self
            .mates
            .iter()
            .map(|mate| {
                std::path::Path::new(&mate.input.path)
                    .file_name()
                    .and_then(|name| name.to_str())
                    .unwrap_or(&mate.input.path)
                    .to_string()
            })
            .collect::<Vec<_>>()
            .join(" | ");
        let records_label = {
            let (total, exact) = self.mates[0].total_reads();
            let total = total.to_formatted_string(&Locale::en);
            if exact {
                total
            } else if self.mates[0].is_stdin() {
                format!("{total}+")
            } else {
                format!("{total}…")
            }
        };
        let worker_error = self.worker_error();

        let scope_hint = self.search_scope.label();
        let footer = if self.search_mode {
            format!("/{}▌", self.search_input)
        } else if self.search_pattern.as_ref().is_some_and(|p| p.is_peptide()) {
            format!(
                "{} | n/N: prev/next | i: {} | /: Search | q: Quit",
                help_text, scope_hint
            )
        } else if self.search_pattern.is_some() {
            format!(
                "{} | n/N: prev/next | x: {} | i: {} | /: Search | q: Quit",
                help_text,
                if self.search_include_rc {
                    "+RC"
                } else {
                    "exact"
                },
                scope_hint
            )
        } else {
            help_text.to_string()
        };

        self.terminal.draw(|f| {
            let full_area = f.area();

            // If showing stats or help, use full screen for them
            if show_stats {
                // Select stats based on active tab
                let stats_lock = self.mates[stats_mate]
                    .stats
                    .lock()
                    .unwrap_or_else(|e| e.into_inner());
                let label = if stats_mate == 1 { "R2" } else { "R1" };

                // Calculate layout for stats blocks
                let (title_area, main_areas) = calculate_stats_layout(full_area);

                // Tab bar + title
                let status_indicator = if stats_lock.scanned_all { "" } else { "..." };

                let tab_bar = if is_paired {
                    let r1_tab = if self.stats_tab == 0 {
                        Span::styled(
                            " [R1] ",
                            Style::default().fg(Color::White).bg(Color::DarkGray),
                        )
                    } else {
                        Span::styled("  R1  ", Style::default().fg(Color::DarkGray))
                    };
                    let r2_tab = if self.stats_tab == 1 {
                        Span::styled(
                            " [R2] ",
                            Style::default().fg(Color::White).bg(Color::DarkGray),
                        )
                    } else {
                        Span::styled("  R2  ", Style::default().fg(Color::DarkGray))
                    };
                    Line::from(vec![
                        r1_tab,
                        Span::raw("  "),
                        r2_tab,
                        Span::raw(format!(
                            "  Sequencing Statistics ({}) {}",
                            label, status_indicator
                        )),
                    ])
                } else {
                    Line::from(Span::styled(
                        format!("Sequencing Statistics {}", status_indicator),
                        Style::default(),
                    ))
                };
                f.render_widget(tab_bar, title_area);

                // Basic stats block
                if let Some(basic_area) = main_areas.first().and_then(|areas| areas.first()) {
                    let basic_stats_content = vec![
                        Line::from(format!(
                            "Total reads: {}{}",
                            stats_total.to_formatted_string(&Locale::en),
                            match (stats_total_exact, stats_from_stdin) {
                                (true, _) => "",
                                (false, true) => "+ (loaded from stdin so far)",
                                (false, false) => " (counting...)",
                            }
                        )),
                        Line::from(format!(
                            "Processed: {} {}",
                            stats_lock.processed_reads.to_formatted_string(&Locale::en),
                            if stats_lock.scanned_all {
                                "(complete)"
                            } else if stats_from_stdin {
                                "(loaded reads)"
                            } else if stats_lock.processed_reads > 0 {
                                "(updating...)"
                            } else {
                                "(starting...)"
                            }
                        )),
                        Line::from(format!("Avg length: {:.1}", stats_lock.avg_length)),
                        Line::from(format!("Min length: {}", stats_lock.min_length)),
                        Line::from(format!("Max length: {}", stats_lock.max_length)),
                        Line::from(format!("GC content: {:.1}%", stats_lock.gc_content)),
                        Line::from(format!("N content: {:.2}%", stats_lock.n_content)),
                    ];
                    let title = Line::from(Span::styled("Basic Statistics", Style::default()));
                    let basic_stats_block = Paragraph::new(basic_stats_content)
                        .block(Block::default().title(title).borders(Borders::ALL))
                        .wrap(Wrap { trim: true });

                    f.render_widget(basic_stats_block, *basic_area);
                }

                // Quality histogram block (left side)
                if let Some(histogram_area) = main_areas.get(1).and_then(|areas| areas.first()) {
                    if !stats_lock.quality_available && stats_lock.processed_reads > 0 {
                        f.render_widget(
                            no_quality_message("Quality Score Distribution"),
                            *histogram_area,
                        );
                    } else if !stats_lock.quality_histogram.is_empty() {
                        let quality_chart =
                            create_quality_histogram(&stats_lock, &self.color_scheme);
                        f.render_widget(quality_chart, *histogram_area);
                    } else {
                        let title = Line::from(Span::styled(
                            "Quality Score Distribution",
                            Style::default(),
                        ));
                        let calculating = Paragraph::new("Calculating quality histogram...")
                            .block(Block::default().title(title).borders(Borders::ALL))
                            .alignment(Alignment::Center);
                        f.render_widget(calculating, *histogram_area);
                    }
                }

                // Average read quality histogram block (right side)
                if let Some(average_histogram_area) =
                    main_areas.get(1).and_then(|areas| areas.get(1))
                {
                    if !stats_lock.quality_available && stats_lock.processed_reads > 0 {
                        f.render_widget(
                            no_quality_message("Average Read Quality Distribution"),
                            *average_histogram_area,
                        );
                    } else if stats_lock.processed_reads > 0 {
                        let average_quality_chart =
                            create_average_read_quality_histogram(&stats_lock, &self.color_scheme);
                        f.render_widget(average_quality_chart, *average_histogram_area);
                    } else {
                        let title = Line::from(Span::styled(
                            "Average Read Quality Distribution",
                            Style::default(),
                        ));
                        let calculating =
                            Paragraph::new("Calculating average quality histogram...")
                                .block(Block::default().title(title).borders(Borders::ALL))
                                .alignment(Alignment::Center);
                        f.render_widget(calculating, *average_histogram_area);
                    }
                }

                // Adapter stats block
                if let Some(adapter_area) = main_areas.first().and_then(|areas| areas.get(1)) {
                    let adapter_stats_block = create_adapter_stats_display(&stats_lock);
                    f.render_widget(adapter_stats_block, *adapter_area);
                }

                // Position quality block
                if let Some(position_area) = main_areas.get(2).and_then(|areas| areas.first()) {
                    if !stats_lock.quality_available && stats_lock.processed_reads > 0 {
                        f.render_widget(
                            no_quality_message("Average Quality by Position"),
                            *position_area,
                        );
                    } else if !stats_lock.position_quality.is_empty() {
                        let position_chart =
                            create_position_quality_chart(&stats_lock, &self.color_scheme);
                        f.render_widget(position_chart, *position_area);
                    } else {
                        let title = Line::from(Span::styled(
                            "Average Quality by Position",
                            Style::default(),
                        ));
                        let calculating = Paragraph::new("Calculating position quality...")
                            .block(Block::default().title(title).borders(Borders::ALL))
                            .alignment(Alignment::Center);
                        f.render_widget(calculating, *position_area);
                    }
                }

                return; // Exit early, don't render main content
            }

            if show_help {
                // Full screen help view
                let help_chunks = Layout::default()
                    .direction(Direction::Vertical)
                    .constraints([
                        Constraint::Length(1), // Status line
                        Constraint::Min(0),    // Help content
                        Constraint::Length(1), // Help line
                    ])
                    .split(full_area);

                // Status line
                let status = Paragraph::new("fqless - Help Screen (Press 'h' to exit)".to_string())
                    .style(Style::default());

                f.render_widget(status, help_chunks[0]);

                // Full screen help content
                let help_content = vec![
                    Line::from(Span::styled(
                        format!("FQLESS - FastQ File Viewer v{}", env!("CARGO_PKG_VERSION")),
                        Style::default(),
                    )),
                    Line::from(""),
                    Line::from("Navigation:"),
                    Line::from("  ↑/k        - Move up one record"),
                    Line::from("  ↓/j        - Move down one record"),
                    Line::from("  PgUp/K     - Move up one page"),
                    Line::from("  PgDn/J/Spc - Move down one page"),
                    Line::from("  Home       - Go to first record"),
                    Line::from(""),
                    Line::from("Display Options:"),
                    Line::from("  S          - Toggle no-wrap mode (horizontal scrolling)"),
                    Line::from("  ←/→        - Scroll left/right (no-wrap mode only)"),
                    Line::from("  c          - Cycle color schemes"),
                    Line::from("  s          - Toggle statistics panel"),
                    Line::from("  p          - Toggle base quality display"),
                    Line::from("  h          - Toggle this help screen"),
                    Line::from(""),
                    Line::from("Orientation & Translation:"),
                    Line::from(
                        "  r          - Cycle read orientation (5'->3', reverse complement)",
                    ),
                    Line::from("  t          - Toggle 6-frame translation"),
                    Line::from(""),
                    Line::from("Search:"),
                    Line::from(
                        "  /          - Search (exact, case-insensitive, matches in yellow)",
                    ),
                    Line::from("              In translation view: amino acid sequence (e.g. CAR)"),
                    Line::from("  n          - Jump to next match"),
                    Line::from("  N          - Jump to previous match"),
                    Line::from("  x          - Toggle reverse-complement search (sequence view)"),
                    Line::from("  i          - Cycle search scope (sequence / read name / both)"),
                    Line::from("  Esc        - Cancel search"),
                    Line::from("  Enter      - Commit search and jump to first match"),
                    Line::from(""),
                    Line::from("Quality & Color:"),
                    Line::from("  e          - Cycle phred score encoding range"),
                    Line::from(""),
                    Line::from("Other:"),
                    Line::from("  q          - Quit"),
                    Line::from(""),
                    Line::from(""),
                    Line::from("Statistics Panel:"),
                    Line::from(
                        " Displays read statistics, quality histograms, and adapter contamination",
                    ),
                    Line::from(" In paired mode, TAB switches between R1 and R2 stats"),
                ];

                let visible_height = help_chunks[1].height.saturating_sub(2) as usize; // Account for borders
                let total_lines = help_content.len();
                let max_scroll = total_lines.saturating_sub(visible_height);
                let actual_scroll = help_scroll.min(max_scroll);

                let visible_content: Vec<Line> = help_content
                    .into_iter()
                    .skip(actual_scroll)
                    .take(visible_height)
                    .collect();

                let help_panel = Paragraph::new(visible_content).block(
                    Block::default()
                        .borders(Borders::ALL)
                        .title(format!("Help (Scroll: {}/{})", actual_scroll, max_scroll))
                        .border_style(Style::default()),
                );

                f.render_widget(help_panel, help_chunks[1]);

                // Help line
                let help_footer: Paragraph<'_> = Paragraph::new(
                    "↑/↓: Scroll | PgUp/PgDn: Page | Home: Top | h: Hide Help | q: Quit",
                )
                .style(Style::default().fg(Color::DarkGray));

                f.render_widget(help_footer, help_chunks[2]);

                return; // Exit early, don't render main content
            }

            // Normal main content view
            let main_area = Layout::default()
                .direction(Direction::Horizontal)
                .constraints([Constraint::Percentage(100)])
                .split(full_area);

            // Main content area
            let main_chunks = Layout::default()
                .direction(Direction::Vertical)
                .constraints([
                    Constraint::Length(1), // Status line
                    Constraint::Min(0),    // Main content
                    Constraint::Length(1), // Help line
                ])
                .split(main_area[0]);

            // Status line
            let status = Paragraph::new(format!(
                "File: {} | Records: {} | Pos: {} | {} | {} | {} {}-{}{}{}",
                file_label,
                records_label,
                current_position,
                mode_display,
                self.orientation.name(),
                self.phred_range.name(),
                self.phred_range.base_phred(),
                self.phred_range.top_phred(),
                if self.show_translation {
                    " | 6-FRAME"
                } else {
                    ""
                },
                search_status
            ))
            .style(Style::default().bg(Color::DarkGray).fg(Color::White));

            f.render_widget(status, main_chunks[0]);

            // Main content
            let main_content = if no_wrap {
                Paragraph::new(prepared_lines)
            } else {
                Paragraph::new(prepared_lines).wrap(Wrap { trim: false })
            };

            f.render_widget(main_content, main_chunks[1]);

            // Help line, or the error of a failed background scan: workers
            // cannot print while the TUI owns the screen.
            let help: Paragraph<'_> = match &worker_error {
                Some(error) => Paragraph::new(format!("Error: {error}"))
                    .style(Style::default().fg(Color::White).bg(Color::Red)),
                None => Paragraph::new(footer.clone()).style(Style::default().fg(Color::DarkGray)),
            };

            f.render_widget(help, main_chunks[2]);
        })?;

        Ok(())
    }
}

impl Drop for TuiViewer {
    fn drop(&mut self) {
        // Workers stop and join when the mates are dropped; stop them before
        // restoring the terminal so the exit is not delayed afterwards.
        self.mates.clear();

        let _ = write!(
            self.terminal.backend_mut(),
            "{}{}{}",
            termion::clear::All,
            termion::cursor::Goto(1, 1),
            termion::cursor::Show
        );
        let _ = self.terminal.show_cursor();
        let _ = self.terminal.flush();
    }
}

#[cfg(test)]
mod bench {
    //! Throughput of the background scans on a real file:
    //! `FQLESS_BENCH=reads.fq cargo test --release bench -- --ignored --nocapture`
    use super::*;
    use std::time::Instant;

    fn path() -> Option<String> {
        std::env::var("FQLESS_BENCH").ok()
    }

    fn report(name: &str, start: Instant, path: &str) {
        let secs = start.elapsed().as_secs_f64();
        let mb = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0) as f64 / 1e6;
        eprintln!("{name:<14} {secs:>6.2}s  {:>7.0} MB/s", mb / secs);
    }

    #[test]
    #[ignore]
    fn bench_scans() {
        let Some(path) = path() else { return };
        let input = Input {
            path: path.clone(),
            reference: std::env::var("FQLESS_BENCH_REFERENCE").ok(),
        };
        let stop = AtomicBool::new(false);

        let start = Instant::now();
        let (total, done) = (AtomicU64::new(0), AtomicBool::new(false));
        count_reads(&input, &total, &done, &stop).unwrap();
        report("count", start, &path);

        let start = Instant::now();
        let stats = Mutex::new(FastqStats::default());
        scan_stats(&input, PhredRange::Default, &stats, &stop).unwrap();
        report("stats", start, &path);

        for query in ["ACGTACGTAC", "read12345"] {
            let pattern = Pattern::new(query, true).unwrap();
            let scope = if pattern.is_nucleotide_literal() {
                SearchScope::Sequence
            } else {
                SearchScope::Id
            };
            let index = Mutex::new(SearchIndex::new());
            let (scanned, done) = (AtomicU64::new(0), AtomicBool::new(false));
            let start = Instant::now();
            scan_search(&input, &pattern, scope, &index, &scanned, &done, &stop).unwrap();
            report(&format!("search {query}"), start, &path);
        }
        eprintln!("reads: {}", total.load(Ordering::Relaxed));
    }
}
