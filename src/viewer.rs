use crate::adapter::{AdapterDetector, AdapterStats};
use crate::buffer::DisplayBuffer;
use crate::color::ColorScheme;
use crate::reader::FastqReader;
use crate::search::{Pattern, SearchIndex};
use crate::utils::{
    calculate_record_lines, determine_min_max_phred, translate_frames_with_quality, PhredRange,
    ReadOrientation,
};
use anyhow::Result;
use bio::alphabets::dna;
use bio::io::fastq;
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
use std::io::stdin;
use std::os::unix::io::AsRawFd;
use std::sync::mpsc::{self, Receiver};
use std::sync::{
    atomic::{AtomicBool, AtomicU64, Ordering},
    Arc, Mutex,
};
use std::thread::{self, JoinHandle};
use std::time::Duration;
use std::{cmp::min, sync::RwLock};
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

#[derive(Debug, Clone)]
pub struct FastqStats {
    pub total_reads: u64,
    pub processed_reads: u64,
    pub avg_length: f64,
    pub min_length: usize,
    pub max_length: usize,
    pub quality_histogram: Vec<u64>, // Histogram of individual base quality scores (0-93)
    pub position_quality: Vec<f64>,  // Average quality at each position
    pub gc_content: f64,
    pub n_content: f64,
    pub scanned_all: bool,
    pub average_read_qualities: Vec<f64>, // Average quality score for each read (one value per read)
    pub adapter_stats: AdapterStats,      // Adapter contamination statistics
}

impl Default for FastqStats {
    fn default() -> Self {
        Self {
            total_reads: 0,
            processed_reads: 0,
            avg_length: 0.0,
            min_length: usize::MAX,
            max_length: 0,
            quality_histogram: vec![0; 94], // Quality scores 0-93
            position_quality: Vec::new(),
            gc_content: 0.0,
            n_content: 0.0,
            scanned_all: false,
            average_read_qualities: Vec::new(),
            adapter_stats: AdapterStats::default(),
        }
    }
}

pub struct TuiViewer {
    terminal: Terminal<TermionBackend<AlternateScreen<RawTerminal<std::io::Stdout>>>>,
    buffer: DisplayBuffer,
    file_path: String,
    current_position: u64,
    horizontal_offset: usize,
    no_wrap: bool,
    show_quality: bool,
    show_stats: bool,
    show_help: bool,
    orientation: ReadOrientation,
    show_translation: bool,
    stats_scroll: usize,
    help_scroll: usize,
    phred_range: PhredRange,
    stats: Arc<Mutex<FastqStats>>,
    stats_worker_handle: Option<JoinHandle<()>>,
    stats_stop_flag: Arc<AtomicBool>,
    color_scheme: ColorScheme,
    search_mode: bool,
    search_input: String,
    search_query: String,
    search_pattern: Option<Pattern>,
    search_include_rc: bool,
    search_index: Arc<Mutex<SearchIndex>>,
    search_indexed_up_to: u64,
    search_match_k: u64,
    search_worker_handle: Option<JoinHandle<()>>,
    search_stop_flag: Arc<AtomicBool>,
    search_worker_scanned: Arc<AtomicU64>,
    search_worker_done: Arc<AtomicBool>,
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
    // Create bins for average read qualities with bin width of 2
    let bin_width = 5;
    let max_qscore = 40; // Maximum quality score to consider
    let num_bins = (max_qscore / bin_width) + 1; // +1 for the last "40+" bin
    let mut bins = vec![0u64; num_bins];

    for &avg_qual in &stats.average_read_qualities {
        let bin_index = ((avg_qual as u8) / bin_width as u8).min((num_bins - 1) as u8) as usize;
        bins[bin_index] += 1;
    }

    let bars: Vec<Bar> = (0..num_bins)
        .map(|i| {
            let range_start = i * bin_width;
            let count = bins[i];

            let label = if i == num_bins - 1 {
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
            "Total adapter instances: {} (avg {:.1} per contaminated read)",
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
                "  {}: {} ({}x, pos: {:.1})",
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
    fn calculate_record_lines(&self, record: &fastq::Record, terminal_width: usize) -> usize {
        calculate_record_lines(
            record,
            terminal_width,
            self.no_wrap,
            self.show_quality,
            self.show_translation,
        )
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

        // Read from the already loaded buffer without triggering additional loads
        let reads_guard = self.buffer.reads.read().unwrap();
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

        // We need to work backwards to find where to start so that moving forward
        // From current position we move bacjwards -1 and sum up the lines
        // needed until we reach the terminal height
        let mut lines_used = 0;
        let mut records_to_move = 0;
        let reads_guard = self.buffer.reads.read().unwrap();
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

    pub fn new(file_path: String) -> Result<Self> {
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

        let buffer = DisplayBuffer::new(&file_path)?;

        let mut viewer = TuiViewer {
            terminal,
            buffer,
            file_path: file_path.clone(),
            current_position: 0,
            horizontal_offset: 0,
            no_wrap: false,
            show_quality: false,
            show_stats: false,
            show_help: false,
            orientation: ReadOrientation::AsIs,
            show_translation: false,
            stats_scroll: 0,
            help_scroll: 0,
            phred_range: PhredRange::Default,
            stats: Arc::new(Mutex::new(FastqStats::default())),
            stats_worker_handle: None,
            stats_stop_flag: Arc::new(AtomicBool::new(false)),
            color_scheme: ColorScheme::RedGreen,
            search_mode: false,
            search_input: String::new(),
            search_query: String::new(),
            search_pattern: None,
            search_include_rc: true,
            search_index: Arc::new(Mutex::new(SearchIndex::new())),
            search_indexed_up_to: 0,
            search_match_k: 0,
            search_worker_handle: None,
            search_stop_flag: Arc::new(AtomicBool::new(false)),
            search_worker_scanned: Arc::new(AtomicU64::new(0)),
            search_worker_done: Arc::new(AtomicBool::new(true)),
        };

        // Start background statistics calculation
        viewer.start_stats_worker();

        Ok(viewer)
    }

    // Start background statistics calculation
    pub fn start_stats_worker(&mut self) {
        // Stop any existing worker first
        self.stop_stats_worker();

        // Create new stop flag for the new worker
        self.stats_stop_flag = Arc::new(AtomicBool::new(false));

        let stats_clone = Arc::clone(&self.stats);
        let phred_range = self.phred_range;
        let stop_flag = Arc::clone(&self.stats_stop_flag);
        let reads = Arc::clone(&self.buffer.reads);
        let file_path = self.file_path.clone();

        // Start fast read counter thread first
        let stats_for_counter = Arc::clone(&self.stats);
        let file_path_for_counter = file_path.clone();
        let stop_flag_for_counter = Arc::clone(&self.stats_stop_flag);
        thread::spawn(move || {
            if let Err(e) = Self::count_reads_fast(
                stats_for_counter,
                file_path_for_counter.as_str(),
                stop_flag_for_counter,
            ) {
                eprintln!("Read counting error: {}", e);
            }
        });

        // Start main stats worker
        let handle = thread::spawn(move || {
            if let Err(e) = Self::calculate_stats_background(
                stats_clone,
                file_path.as_str(),
                reads,
                phred_range,
                stop_flag,
            ) {
                eprintln!("Stats calculation error: {}", e);
            }
        });

        self.stats_worker_handle = Some(handle);
    }

    // Stop the stats worker if it's running
    fn stop_stats_worker(&mut self) {
        if let Some(handle) = self.stats_worker_handle.take() {
            // Signal the worker to stop
            self.stats_stop_flag.store(true, Ordering::SeqCst);
            // Wait for it to finish
            let _ = handle.join();
        }
    }

    // Reset stats to default (useful when changing phred range)
    fn reset_stats(&mut self) {
        if let Ok(mut stats_lock) = self.stats.lock() {
            *stats_lock = FastqStats::default();
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
                let target = self
                    .search_index
                    .lock()
                    .unwrap()
                    .next(self.current_position);
                if let Some(m) = target {
                    let k = self.search_index.lock().unwrap().count_up_to(m);
                    self.search_match_k = k;
                    let _ = self.jump_to_record(m);
                }
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
        self.search_index.lock().unwrap().reset();
        self.search_indexed_up_to = 0;
        self.search_match_k = 0;
        self.search_pattern = Some(pattern);
        self.index_loaded_records();
        self.start_search_worker();
    }

    /// Clear the active search and stop the look-ahead worker.
    fn clear_search(&mut self) {
        self.stop_search_worker();
        self.search_pattern = None;
        self.search_query.clear();
        self.search_include_rc = true;
        self.search_index.lock().unwrap().reset();
        self.search_indexed_up_to = 0;
        self.search_match_k = 0;
    }

    /// Scan newly loaded records on the main thread and mark matches in the
    /// index. Cheap per frame: only the delta since the last scan is visited.
    fn index_loaded_records(&mut self) {
        let Some(pattern) = self.search_pattern.clone() else {
            self.search_indexed_up_to = self.buffer.reads.read().unwrap().len() as u64;
            return;
        };
        let reads = self.buffer.reads.read().unwrap();
        let len = reads.len() as u64;
        let mut index = self.search_index.lock().unwrap();
        for record_idx in self.search_indexed_up_to..len {
            if pattern.has_match(reads[record_idx as usize].seq()) {
                index.set(record_idx);
            }
        }
        self.search_indexed_up_to = len;
    }

    /// Jump to a record index, streaming forward if it is not loaded yet.
    fn jump_to_record(&mut self, record: u64) -> Result<()> {
        self.current_position = record;
        self.buffer
            .load_window(self.current_position, BUFFER_WINDOW_SIZE)?;
        let len = self.buffer.reads.read().unwrap().len() as u64;
        if len > 0 && self.current_position >= len {
            self.current_position = len - 1;
        }
        Ok(())
    }

    /// Move to the next (or previous) matching record relative to the current
    /// position. Returns Ok(false) if no further match is known.
    fn next_match(&mut self, forward: bool) -> Result<()> {
        if self.search_pattern.is_none() {
            return Ok(());
        }
        let index = self.search_index.lock().unwrap();
        let target = if forward {
            index.next(self.current_position + 1)
        } else {
            index.prev(self.current_position)
        };
        let k = target.map(|m| index.count_up_to(m));
        drop(index);
        if let Some(m) = target {
            self.search_match_k = k.unwrap_or(0);
            self.jump_to_record(m)?;
        }
        Ok(())
    }

    /// Start (or restart) the background look-ahead search worker. Scans the
    /// whole file on its own reader, setting index bits as it goes. Skipped
    /// for stdin, which cannot be re-read.
    fn start_search_worker(&mut self) {
        self.stop_search_worker();
        let Some(pattern) = self.search_pattern.clone() else {
            return;
        };
        if self.file_path == "-" {
            self.search_worker_done.store(true, Ordering::SeqCst);
            return;
        }

        self.search_stop_flag = Arc::new(AtomicBool::new(false));
        self.search_worker_done.store(false, Ordering::SeqCst);
        self.search_worker_scanned.store(0, Ordering::SeqCst);

        let index = Arc::clone(&self.search_index);
        let scanned = Arc::clone(&self.search_worker_scanned);
        let done = Arc::clone(&self.search_worker_done);
        let stop_flag = Arc::clone(&self.search_stop_flag);
        let file_path = self.file_path.clone();

        let handle = thread::spawn(move || {
            if let Err(e) =
                Self::search_background(&pattern, &file_path, index, scanned, done, stop_flag)
            {
                eprintln!("Search error: {}", e);
            }
        });
        self.search_worker_handle = Some(handle);
    }

    // Stop the search worker if it is running
    fn stop_search_worker(&mut self) {
        if let Some(handle) = self.search_worker_handle.take() {
            self.search_stop_flag.store(true, Ordering::SeqCst);
            let _ = handle.join();
        }
    }

    /// Background search: stream the file once, setting index bits for every
    /// record containing a match. This provides look-ahead past the loaded
    /// buffer without storing anything but the one-bit-per-record index.
    fn search_background(
        pattern: &Pattern,
        file_path: &str,
        index: Arc<Mutex<SearchIndex>>,
        scanned: Arc<AtomicU64>,
        done: Arc<AtomicBool>,
        stop_flag: Arc<AtomicBool>,
    ) -> Result<()> {
        let reader = FastqReader::new(file_path)?;
        let mut record_idx = 0u64;
        for record in reader.into_fastq_reader().records() {
            if stop_flag.load(Ordering::Relaxed) {
                return Ok(());
            }
            let record = record?;
            if pattern.has_match(record.seq()) {
                index.lock().unwrap().set(record_idx);
            }
            record_idx += 1;
            if record_idx.is_multiple_of(10000) {
                scanned.store(record_idx, Ordering::Relaxed);
            }
        }
        scanned.store(record_idx, Ordering::SeqCst);
        done.store(true, Ordering::SeqCst);
        Ok(())
    }

    // Fast read counter - just counts total reads without detailed analysis
    fn count_reads_fast(
        stats: Arc<Mutex<FastqStats>>,
        file_path: &str,
        stop_flag: Arc<AtomicBool>,
    ) -> Result<()> {
        // Skip counting for stdin (can't read twice)
        if file_path == "-" {
            return Ok(());
        }

        let reader = FastqReader::new(file_path)?;
        let mut total_reads = 0u64;

        for record in reader.into_fastq_reader().records() {
            let _ = record?;

            if stop_flag.load(Ordering::Relaxed) {
                return Ok(());
            }

            total_reads += 1;

            // Update count every n reads for faster feedback
            if total_reads.is_multiple_of(100000) {
                if let Ok(mut stats_lock) = stats.lock() {
                    stats_lock.total_reads = total_reads;
                }
            }
        }

        // Final update
        if !stop_flag.load(Ordering::Relaxed) {
            if let Ok(mut stats_lock) = stats.lock() {
                stats_lock.total_reads = total_reads;
            }
        }

        Ok(())
    }

    fn calculate_stats_background(
        stats: Arc<Mutex<FastqStats>>,
        file_path: &str,
        reads: Arc<RwLock<Vec<fastq::Record>>>,
        phred_range: PhredRange,
        stop_flag: Arc<AtomicBool>,
    ) -> Result<()> {
        let mut total_reads = 0u64;
        let mut total_length = 0u64;
        let mut min_length = usize::MAX;
        let mut max_length = 0usize;
        let mut quality_histogram =
            vec![0u64; (phred_range.top_phred() - phred_range.base_phred() + 1) as usize];
        let mut position_quality_sums = Vec::new();
        let mut position_counts = Vec::new();
        let mut total_gc = 0u64;
        let mut total_n = 0u64;
        let mut average_read_qualities = Vec::new();

        // Initialize adapter detector
        let adapter_detector = AdapterDetector::new();
        let mut local_adapter_stats = AdapterStats::default();

        // if file_path is "-", we are reading from stdin so we use the arc reads
        // otherwise we make a new FastqReader
        let is_stdin = file_path == "-";

        let reads: Box<dyn Iterator<Item = Result<bio::io::fastq::Record, bio::io::fastq::Error>>> =
            if is_stdin {
                let guard = reads.read().unwrap();
                let records = guard.clone().into_iter().map(Ok);
                Box::new(records)
            } else {
                let reader = FastqReader::new(file_path)?;
                Box::new(reader.into_fastq_reader().records())
            };
        for record in reads {
            let record = record?;
            if stop_flag.load(Ordering::Relaxed) {
                return Ok(());
            }

            total_reads += 1;
            let seq_len = record.seq().len();
            total_length += seq_len as u64;
            min_length = min_length.min(seq_len);
            max_length = max_length.max(seq_len);

            let seq = record.seq();
            let qual = record.qual();

            // Combined pass through sequence for GC/N content and adapter detection
            // This is more cache-friendly than separate passes
            let mut gc_count = 0u64;
            let mut n_count = 0u64;

            // Fast GC/N counting
            for &base in seq {
                match base {
                    b'G' | b'g' | b'C' | b'c' => gc_count += 1,
                    b'N' | b'n' => n_count += 1,
                    _ => {}
                }
            }

            total_gc += gc_count;
            total_n += n_count;

            // Detect adapters in every read
            let adapter_detections = adapter_detector.detect_adapters(seq);
            adapter_detector.update_stats(&mut local_adapter_stats, &adapter_detections);

            // Process quality scores with optimized loop
            let base_phred = phred_range.base_phred();
            let mut total_qual = 0u32;

            // Ensure position vectors are large enough upfront
            if seq_len > position_quality_sums.len() {
                position_quality_sums.resize(seq_len, 0.0);
                position_counts.resize(seq_len, 0);
            }

            for (pos, &qual_val) in qual.iter().enumerate() {
                let quality_score = qual_val.saturating_sub(base_phred);
                total_qual += quality_score as u32;

                // Update histogram
                if (quality_score as usize) < quality_histogram.len() {
                    quality_histogram[quality_score as usize] += 1;
                }

                // Update position-wise quality (no bounds checking needed now)
                position_quality_sums[pos] += quality_score as f64;
                position_counts[pos] += 1;
            }

            let average_quality = if seq_len > 0 {
                total_qual as f64 / seq_len as f64
            } else {
                0.0
            };

            average_read_qualities.push(average_quality);

            // Update stats every so often and check stop flag
            if total_reads.is_multiple_of(100000) {
                if stop_flag.load(Ordering::Relaxed) {
                    return Ok(());
                }

                // Update all stats including histograms, but only every 100k reads
                if let Ok(mut stats_lock) = stats.lock() {
                    stats_lock.processed_reads = total_reads;
                    stats_lock.avg_length = total_length as f64 / total_reads as f64;
                    stats_lock.min_length = if min_length == usize::MAX {
                        0
                    } else {
                        min_length
                    };
                    stats_lock.max_length = max_length;
                    stats_lock.gc_content = (total_gc as f64 / total_length as f64) * 100.0;
                    stats_lock.n_content = (total_n as f64 / total_length as f64) * 100.0;

                    // Update histograms (still expensive but only every 100k)
                    stats_lock.quality_histogram = quality_histogram.clone();
                    stats_lock.average_read_qualities = average_read_qualities.clone();
                    stats_lock.adapter_stats = local_adapter_stats.clone();

                    // Calculate position quality
                    stats_lock.position_quality = position_quality_sums
                        .iter()
                        .zip(position_counts.iter())
                        .map(|(sum, count)| {
                            if *count > 0 {
                                sum / (*count as f64)
                            } else {
                                0.0
                            }
                        })
                        .collect();
                }
            }
        }

        // Final update (only if not stopped)
        if !stop_flag.load(Ordering::Relaxed) {
            if let Ok(mut stats_lock) = stats.lock() {
                // Only update total_reads if fast counter didn't run (stdin case)
                if file_path == "-" {
                    stats_lock.total_reads = total_reads;
                }
                stats_lock.processed_reads = total_reads;
                stats_lock.avg_length = if total_reads > 0 {
                    total_length as f64 / total_reads as f64
                } else {
                    0.0
                };
                stats_lock.min_length = if min_length == usize::MAX {
                    0
                } else {
                    min_length
                };
                stats_lock.max_length = max_length;
                stats_lock.quality_histogram = quality_histogram;
                stats_lock.gc_content = if total_length > 0 {
                    (total_gc as f64 / total_length as f64) * 100.0
                } else {
                    0.0
                };
                stats_lock.n_content = if total_length > 0 {
                    (total_n as f64 / total_length as f64) * 100.0
                } else {
                    0.0
                };

                stats_lock.position_quality = position_quality_sums
                    .iter()
                    .zip(position_counts.iter())
                    .map(|(sum, count)| {
                        if *count > 0 {
                            sum / (*count as f64)
                        } else {
                            0.0
                        }
                    })
                    .collect();
                stats_lock.average_read_qualities = average_read_qualities;
                stats_lock.adapter_stats = local_adapter_stats;
                stats_lock.scanned_all = true;
            }
        }
        Ok(())
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
    fn aa_span<'a>(&self, ch: u8, quality: u8, highlighted: bool) -> Span<'a> {
        if highlighted {
            Span::styled(
                char::from(ch).to_string(),
                Style::default().fg(Color::Black).bg(Color::Yellow),
            )
        } else {
            let color = if ch == b'*' {
                Color::Red
            } else {
                self.color_scheme.quality_to_color(quality)
            };
            Span::styled(char::from(ch).to_string(), Style::default().fg(color))
        }
    }

    fn colorize_sequence<'a>(&self, sequence: &[u8], quality: &[u8]) -> Vec<Span<'a>> {
        let ranges = self
            .highlight_pattern()
            .map(|p| p.find_merged_matches(sequence))
            .unwrap_or_default();
        let mut spans = Vec::with_capacity(sequence.len());
        let mut pos = 0;
        for (start, end) in &ranges {
            for (i, &base) in sequence[pos..*start].iter().enumerate() {
                spans.push(self.base_span(base, quality[pos + i], false));
            }
            for &base in &sequence[*start..*end] {
                spans.push(self.base_span(base, 0, true));
            }
            pos = *end;
        }
        for (i, &base) in sequence[pos..].iter().enumerate() {
            spans.push(self.base_span(base, quality[pos + i], false));
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

    /// Return the sequence and quality of a record in the currently selected
    /// orientation. Quality stays paired with its base, so it is reversed
    /// together with the sequence for the reverse complement.
    fn oriented_seq_qual(&self, record: &fastq::Record) -> (Vec<u8>, Vec<u8>) {
        match self.orientation {
            ReadOrientation::AsIs => (record.seq().to_vec(), record.qual().to_vec()),
            ReadOrientation::ReverseComplement => (
                dna::revcomp(record.seq()),
                record.qual().iter().rev().copied().collect(),
            ),
        }
    }

    pub fn run(&mut self) -> Result<()> {
        // Load initial record
        self.buffer
            .load_window(self.current_position, BUFFER_WINDOW_SIZE)?;

        // Determine where to read keyboard input from:
        // - If reading from a file, stdin is available for keyboard input
        // - If reading from stdin (piped), we need to use /dev/tty for keyboard input
        let use_tty = self.file_path == "-";
        let key_source: Box<dyn Iterator<Item = std::io::Result<Key>> + Send> = if use_tty {
            Box::new(File::open("/dev/tty")?.keys())
        } else {
            Box::new(stdin().keys())
        };

        let rx = spawn_input_thread(key_source);

        // use the buffered reads to determine the phred range
        let (phred_min, phred_max) = determine_min_max_phred(&self.buffer.reads.read().unwrap());
        self.phred_range = PhredRange::from_min_max(phred_min, phred_max);

        self.terminal.clear()?;
        let mut last_stats_update = 0;
        let mut needs_redraw = true;
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
                    Key::Char('c') => {
                        // iterate color schemes
                        self.color_scheme = self.color_scheme.next();
                    }
                    Key::Char('s') => {
                        self.show_stats = !self.show_stats;
                        self.stats_scroll = 0; // Reset scroll when toggling
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
                        // Reset stats and restart worker with new phred range
                        self.reset_stats();
                        self.start_stats_worker();
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
                            self.current_position += 1;
                            self.buffer
                                .load_window(self.current_position, BUFFER_WINDOW_SIZE)?;
                            // if records are empty or we reached the end, don't go further
                            if self.current_position
                                >= self.buffer.reads.read().unwrap().len() as u64
                            {
                                self.current_position =
                                    self.buffer.reads.read().unwrap().len() as u64 - 1;
                            }
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
                            self.current_position = new_position;
                            self.buffer
                                .load_window(self.current_position, BUFFER_WINDOW_SIZE)?;
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

            // Redraw if stats are updating (for both stats view and status line in normal view)
            if !had_input {
                if let Ok(stats_lock) = self.stats.try_lock() {
                    if !stats_lock.scanned_all {
                        needs_redraw = true;
                    }
                }
            }

            // Small delay to prevent excessive CPU usage
            std::thread::sleep(Duration::from_millis(if needs_redraw { 100 } else { 250 }));
            last_stats_update += if needs_redraw { 100 } else { 250 };

            // For stdin input, check if stats worker needs to be restarted
            // Only restart if it has finished or crashed, not repeatedly
            if use_tty && last_stats_update >= 200 {
                let worker_finished = self
                    .stats_worker_handle
                    .as_ref()
                    .is_none_or(|h| h.is_finished());

                if worker_finished {
                    self.start_stats_worker();
                }
                last_stats_update = 0;
            }
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

        // Calculate available space first
        let available_height = terminal_size.height.saturating_sub(2) as usize;
        // assuming 2 lines for each record as minimum
        let max_visible = available_height / 2;
        let records = self
            .buffer
            .get_window(self.current_position, max_visible + 5)?; // +5 to ensure we have enough lines loaded
                                                                  // print current position and horizontal offset for debug

        //for i in 0..max_visible {
        for record in records.iter() {
            // remove @ from the ID if it exists
            let name = if let Some(desc) = record.desc() {
                format!("{} {}", record.id(), desc)
            } else {
                record.id().to_string()
            };
            // Header line in cyan
            prepared_lines.push(Line::from(Span::styled(name.to_string(), Style::default())));

            let (oriented_seq, oriented_qual) = self.oriented_seq_qual(record);

            if self.show_translation {
                // 6-frame translation: one line per reading frame, colored by
                // the effective quality of the amino acid. Amino acid search
                // matches are highlighted in yellow.
                let frames = translate_frames_with_quality(
                    &oriented_seq,
                    &oriented_qual,
                    self.phred_range.base_phred(),
                );
                let highlight = self.highlight_pattern();
                // Left-pad the labels so the amino acid columns line up
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
                let label_width = labels
                    .iter()
                    .map(|label| label.chars().count())
                    .max()
                    .unwrap_or(0);
                for (i, (frame, scores)) in frames.iter().enumerate() {
                    let label = format!("{:<width$}", labels[i], width = label_width);
                    let (visible_start, visible_end) = if no_wrap {
                        let start = horizontal_offset.min(frame.len());
                        (
                            start,
                            (start + terminal_width.saturating_sub(1)).min(frame.len()),
                        )
                    } else {
                        (0, frame.len())
                    };
                    let visible = &frame.as_bytes()[visible_start..visible_end];
                    let visible_scores = &scores[visible_start..visible_end];
                    // Amino acid frames are ASCII, so byte offsets line up
                    // with columns.
                    let ranges: Vec<(usize, usize)> = highlight
                        .as_ref()
                        .map(|p| p.find_peptide_matches(frame.as_bytes()))
                        .unwrap_or_default();
                    let mut spans = vec![Span::styled(label, Style::default().fg(Color::DarkGray))];
                    let mut pos = 0;
                    for (start, end) in &ranges {
                        let s = (*start).max(visible_start);
                        let e = (*end).min(visible_end);
                        if s >= e {
                            continue;
                        }
                        for (k, &ch) in visible[pos..s - visible_start].iter().enumerate() {
                            spans.push(self.aa_span(ch, visible_scores[pos + k], false));
                        }
                        for &ch in &visible[s - visible_start..e - visible_start] {
                            spans.push(self.aa_span(ch, 0, true));
                        }
                        pos = e - visible_start;
                    }
                    for (k, &ch) in visible[pos..].iter().enumerate() {
                        spans.push(self.aa_span(ch, visible_scores[pos + k], false));
                    }
                    prepared_lines.push(Line::from(spans));
                }
                continue;
            }

            // Handle sequence display based on wrap mode
            let visible_sequence = if no_wrap {
                let start = horizontal_offset.min(oriented_seq.len());
                let end = (start + terminal_width.saturating_sub(1)).min(oriented_seq.len());
                &oriented_seq[start..end]
            } else {
                &oriented_seq
            };

            let visible_quality = if no_wrap && oriented_qual.len() > horizontal_offset {
                let quality_end =
                    (horizontal_offset + terminal_width.saturating_sub(1)).min(oriented_qual.len());
                &oriented_qual[horizontal_offset..quality_end]
            } else {
                &oriented_qual
            };

            let sequence_spans = if !self.show_quality {
                self.colorize_sequence(visible_sequence, visible_quality)
            } else {
                self.raw_sequence_with_highlight(visible_sequence)
            };

            prepared_lines.push(Line::from(sequence_spans));

            // if show_quality is true, add quality line showing the quality characters
            if self.show_quality {
                let quality_spans = self.colorize_sequence(visible_quality, visible_quality);
                prepared_lines.push(Line::from(quality_spans));
            }
        }

        let wrap_status = if no_wrap { "NO-WRAP" } else { "WRAP" };

        let help_text = if no_wrap {
            "↑/k: Up | ↓/j: Down | ←/→: Scroll | PgUp/PgDn: Page | S: Wrap | /: Search | s: Stats | c: Colors | q: Quit"
        } else {
            "↑/k: Up | ↓/j: Down | PgUp/PgDn: Page | S: No-Wrap | /: Search | s: Stats | c: Colors | q: Quit"
        };

        let search_status = if self.search_pattern.is_some() {
            let count = self.search_index.lock().unwrap().count();
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
            let scanning = if self.search_worker_done.load(Ordering::SeqCst) {
                String::new()
            } else {
                format!(
                    " scanning…({})",
                    self.search_worker_scanned
                        .load(Ordering::SeqCst)
                        .to_formatted_string(&Locale::en)
                )
            };
            format!(
                " | /{} {mode} {}/{}{}",
                self.search_query, self.search_match_k, count, scanning
            )
        } else {
            String::new()
        };

        let footer = if self.search_mode {
            format!("/{}▌", self.search_input)
        } else if self.search_pattern.as_ref().is_some_and(|p| p.is_peptide()) {
            format!("{} | n/N: prev/next | /: Search | q: Quit", help_text)
        } else if self.search_pattern.is_some() {
            format!(
                "{} | n/N: prev/next | x: {} | /: Search | q: Quit",
                help_text,
                if self.search_include_rc {
                    "+RC"
                } else {
                    "exact"
                }
            )
        } else {
            help_text.to_string()
        };

        self.terminal.draw(|f| {
            let full_area = f.area();

            // If showing stats or help, use full screen for them
            if show_stats {
                let stats_lock = self.stats.lock().unwrap_or_else(|e| e.into_inner());

                // Calculate layout for stats blocks
                let (title_area, main_areas) = calculate_stats_layout(full_area);

                // Title
                let is_processing = stats_lock.processed_reads > 0;

                let status_indicator = if is_processing && !stats_lock.scanned_all {
                    "..."
                } else if stats_lock.scanned_all {
                    ""
                } else {
                    "[stalled]"
                };

                let title = Line::from(Span::styled(
                    format!("FASTQ Statistics {}", status_indicator),
                    Style::default(),
                ));
                f.render_widget(title, title_area);

                // Basic stats block
                if let Some(basic_area) = main_areas.first().and_then(|areas| areas.first()) {
                    let basic_stats_content = vec![
                        // format the total reads with spaces for better readability for high numbers
                        // eg. "Total reads: 1 234 567"
                        Line::from(format!(
                            "Total reads: {}",
                            stats_lock.total_reads.to_formatted_string(&Locale::en)
                        )),
                        Line::from(format!(
                            "Processed: {} {}",
                            stats_lock.processed_reads.to_formatted_string(&Locale::en),
                            if stats_lock.processed_reads > 0 && !stats_lock.scanned_all {
                                "(updating...)"
                            } else if stats_lock.scanned_all {
                                "(complete)"
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
                    if !stats_lock.quality_histogram.is_empty() {
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
                    if !stats_lock.average_read_qualities.is_empty() {
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
                    if !stats_lock.position_quality.is_empty() {
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
                    Line::from("  PgUp       - Move up 10 records"),
                    Line::from("  PgDn       - Move down 10 records"),
                    Line::from("  Home       - Go to first record"),
                    Line::from("  End        - Go to last record"),
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
            let filename = std::path::Path::new(&self.file_path)
                .file_name()
                .and_then(|name| name.to_str())
                .unwrap_or(&self.file_path);

            let status = Paragraph::new(format!(
                "File: {} | Records: {} | Pos: {} | {} | {} | {} {}-{}{}{}",
                filename,
                self.stats
                    .lock()
                    .unwrap()
                    .total_reads
                    .to_formatted_string(&Locale::en),
                current_position,
                wrap_status,
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

            // Help line
            let help: Paragraph<'_> =
                Paragraph::new(footer.clone()).style(Style::default().fg(Color::DarkGray));

            f.render_widget(help, main_chunks[2]);
        })?;

        Ok(())
    }
}

impl Drop for TuiViewer {
    fn drop(&mut self) {
        self.stop_stats_worker();
        self.stop_search_worker();

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
