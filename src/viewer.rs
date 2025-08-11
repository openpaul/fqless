use crate::adapter::{AdapterDetector, AdapterStats};
use crate::buffer::DisplayBuffer;
use crate::color::ColorScheme;
use crate::reader::FastqReader;
use anyhow::Result;
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
    atomic::{AtomicBool, Ordering},
    Arc, Mutex,
};
use std::thread::{self, JoinHandle};
use std::time::Duration;
use std::{cmp::min, sync::RwLock};
use std::{fs::File, io::stdout};
use std::{io::Write, os::fd::BorrowedFd};
use termion::{
    event::Key,
    input::{MouseTerminal, TermRead},
    raw::{IntoRawMode, RawTerminal},
    screen::{AlternateScreen, IntoAlternateScreen},
};

#[derive(Debug, Clone)]
pub struct FastqStats {
    pub total_reads: u64,
    pub processed_reads: u64,
    pub avg_length: f64,
    pub min_length: usize,
    pub max_length: usize,
    pub quality_histogram: Vec<u64>, // Histogram of quality scores (0-93)
    pub position_quality: Vec<f64>,  // Average quality at each position
    pub gc_content: f64,
    pub n_content: f64,
    pub scanned_all: bool,
    pub average_read_qualities: Vec<f64>, // Average quality per read
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
    terminal:
        Terminal<TermionBackend<AlternateScreen<MouseTerminal<RawTerminal<std::io::Stdout>>>>>,
    buffer: DisplayBuffer,
    file_path: String,
    current_position: u64,
    horizontal_offset: usize,
    no_wrap: bool,
    show_quality: bool,
    show_stats: bool,
    show_help: bool,
    stats_scroll: usize,
    help_scroll: usize,
    phred_range: PhredRange,
    stats: Arc<Mutex<FastqStats>>,
    stats_worker_handle: Option<JoinHandle<()>>,
    stats_stop_flag: Arc<AtomicBool>,
    color_scheme: ColorScheme,
}

/// Calculate the layout for the stats page with blocks
fn calculate_stats_layout(area: Rect) -> (Rect, Vec<Vec<Rect>>) {
    let main_layout = Layout::vertical([Constraint::Length(1), Constraint::Min(0)]);
    let [title_area, main_area] = main_layout.areas(area);

    // Create a 2x3 grid for stats blocks
    let block_layout = Layout::vertical([
        Constraint::Length(10), // Basic stats and Adapter stats (side by side)
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
    let bars: Vec<Bar> = (0..=40)
        .step_by(2)
        .map(|q| {
            let count = stats.quality_histogram.get(q).unwrap_or(&0);
            Bar::default()
                .value(*count)
                .label(Line::from(format!("Q{}", q)))
                .style(quality_score_style(q as u8, color_scheme))
        })
        .collect();

    let title = Line::from(Span::styled("Quality Score Distribution", Style::default()));
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
            let _range_end = range_start + bin_width;
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
        "Average Read Quality Distribution",
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
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PhredRange {
    Solexa,
    Illumina1_3,
    Illumina1_5,
    Default,
}

impl PhredRange {
    fn range(&self) -> Option<(u8, u8)> {
        match self {
            Self::Default => Some((33, 126)),
            Self::Solexa => Some((59, 104)),
            Self::Illumina1_3 => Some((64, 126)),
            Self::Illumina1_5 => Some((64, 126)),
        }
    }

    fn all() -> &'static [Self] {
        &[
            Self::Solexa,
            Self::Illumina1_3,
            Self::Illumina1_5,
            Self::Default,
        ]
    }

    fn from_min_max(min: u8, max: u8) -> Self {
        for variant in Self::all() {
            if let Some((rmin, rmax)) = variant.range() {
                if min >= rmin && max <= rmax {
                    return *variant;
                }
            }
        }
        Self::Default
    }

    fn base_phred(&self) -> u8 {
        self.range().map(|(min, _)| min).unwrap_or(0)
    }

    fn top_phred(&self) -> u8 {
        self.range().map(|(_, max)| max).unwrap_or(0)
    }

    fn name(&self) -> &'static str {
        match self {
            Self::Solexa => "Solexa",
            Self::Illumina1_3 => "Illumina 1.3",
            Self::Illumina1_5 => "Illumina 1.8",
            Self::Default => "Default",
        }
    }
}

fn determine_min_max_phred(records: &[fastq::Record]) -> (u8, u8) {
    let mut min_phred = u8::MAX;
    let mut max_phred = u8::MIN;

    for record in records {
        for &q in record.qual() {
            min_phred = min_phred.min(q);
            max_phred = max_phred.max(q);
        }
    }

    (min_phred, max_phred)
}

fn stdin_has_data() -> bool {
    let fd = stdin().as_raw_fd();
    let fd = unsafe { BorrowedFd::borrow_raw(fd) };
    let mut fds = [PollFd::new(fd, PollFlags::POLLIN)];
    match poll(&mut fds, nix::poll::PollTimeout::from(Some(250 as u16))) {
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
        for key in key_iter {
            if let Ok(k) = key {
                if tx.send(k).is_err() {
                    break;
                }
            }
        }
    });
    rx
}

/// Create adapter statistics display
fn create_adapter_stats_display(stats: &FastqStats) -> Paragraph {
    let adapter_stats = &stats.adapter_stats;

    let contamination_rate = if stats.total_reads > 0 {
        (adapter_stats.contaminated_reads as f64 / stats.total_reads as f64) * 100.0
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
            "Total adapters found: {}",
            adapter_stats
                .total_adapters_found
                .to_formatted_string(&Locale::en)
        )),
        Line::from(""),
    ];

    // Show top 3 most common adapters
    let mut adapter_pairs: Vec<_> = adapter_stats.adapters_detected.iter().collect();
    adapter_pairs.sort_by(|a, b| b.1.count.cmp(&a.1.count));

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
    pub fn new(file_path: String) -> Result<Self> {
        let running = Arc::new(AtomicBool::new(true));
        flag::register(SIGINT, Arc::clone(&running))?;

        if file_path == "-" && !stdin_has_data() {
            return Err(anyhow::anyhow!(
                "No data available on stdin. Maybe you forgot to pipe input?"
            ));
        }

        let stdout = stdout().into_raw_mode()?;
        let stdout = MouseTerminal::from(stdout);
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
            stats_scroll: 0,
            help_scroll: 0,
            phred_range: PhredRange::Default,
            stats: Arc::new(Mutex::new(FastqStats::default())),
            stats_worker_handle: None,
            stats_stop_flag: Arc::new(AtomicBool::new(false)),
            color_scheme: ColorScheme::RedGreen,
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
                let records = guard.clone().into_iter().map(|r| Ok(r));
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

            // Fast GC/N counting using unsafe for performance
            unsafe {
                let seq_ptr = seq.as_ptr();
                for i in 0..seq.len() {
                    match *seq_ptr.add(i) {
                        b'G' | b'g' | b'C' | b'c' => gc_count += 1,
                        b'N' | b'n' => n_count += 1,
                        _ => {}
                    }
                }
            }

            total_gc += gc_count;
            total_n += n_count;

            // Detect adapters (only every 10th read to reduce overhead)
            if total_reads % 10 == 0 {
                let adapter_detections = adapter_detector.detect_adapters(seq);
                adapter_detector.update_stats(&mut local_adapter_stats, &adapter_detections);
            }

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
                total_qual += qual_val as u32;

                // Update histogram
                if (quality_score as usize) < quality_histogram.len() {
                    quality_histogram[quality_score as usize] += 1;
                }

                // Update position-wise quality (no bounds checking needed now)
                position_quality_sums[pos] += quality_score as f64;
                position_counts[pos] += 1;
            }

            let average_quality = {
                let base_qual_total = seq_len as u32 * base_phred as u32;
                if total_qual >= base_qual_total {
                    (total_qual - base_qual_total) as f64 / seq_len as f64
                } else {
                    0.0 // Handle underflow case
                }
            };

            average_read_qualities.push(average_quality);

            // Update stats every so often and check stop flag
            if total_reads % 100000 == 0 {
                if stop_flag.load(Ordering::Relaxed) {
                    return Ok(());
                }

                if let Ok(mut stats_lock) = stats.lock() {
                    stats_lock.total_reads = total_reads;
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
                    stats_lock.quality_histogram = quality_histogram.clone();
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

                    // Calculate average quality per position
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
                    stats_lock.average_read_qualities = average_read_qualities.clone();
                    stats_lock.adapter_stats = local_adapter_stats.clone();
                }
            }
        }

        // Final update (only if not stopped)
        if !stop_flag.load(Ordering::Relaxed) {
            if let Ok(mut stats_lock) = stats.lock() {
                stats_lock.total_reads = total_reads;
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

    fn colorize_sequence<'a>(&self, sequence: &[u8], quality: &[u8]) -> Vec<Span<'a>> {
        sequence
            .iter()
            .zip(quality.iter())
            .map(|(nucleotide, qual_char)| {
                let quality_score = qual_char.saturating_sub(self.phred_range.base_phred());
                let color = self.color_scheme.quality_to_color(quality_score);
                Span::styled(
                    char::from(*nucleotide).to_string(),
                    Style::default().fg(color),
                )
            })
            .collect()
    }

    pub fn run(&mut self) -> Result<()> {
        // Load initial record
        self.buffer.load_window(self.current_position, 1000)?;

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
        loop {
            // Always draw first to update display
            self.draw()?;

            // Read keyboard input - from stdin for files, /dev/tty for piped input
            while let Ok(key) = rx.try_recv() {
                match key {
                    Key::Char('q') | Key::Ctrl('c') => return Ok(()),
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
                    Key::Char('r') => {
                        // Adjust coloring range by iterating through known ranges
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
                            self.buffer.load_window(self.current_position, 1000)?;
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
                            self.stats_scroll += 10;
                        } else if self.show_help {
                            self.help_scroll += 10;
                        } else {
                            self.current_position += 10;
                            self.buffer.load_window(self.current_position, 1000)?;
                        }
                    }
                    Key::PageUp | Key::Char('K') => {
                        if self.show_stats {
                            self.stats_scroll = self.stats_scroll.saturating_sub(10);
                        } else if self.show_help {
                            self.help_scroll = self.help_scroll.saturating_sub(10);
                        } else {
                            self.current_position = self.current_position.saturating_sub(10);
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

            // Small delay to prevent excessive CPU usage
            std::thread::sleep(Duration::from_millis(16)); // ~60 FPS
            last_stats_update += 16;
            if use_tty && last_stats_update >= 200 {
                self.start_stats_worker();
                last_stats_update = 0;
            }
        }
    }

    fn draw(&mut self) -> Result<()> {
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
            prepared_lines.push(Line::from(Span::styled(
                format!("{}", name),
                Style::default(),
            )));
            // Handle sequence display based on wrap mode
            let visible_sequence = if no_wrap {
                let start = horizontal_offset.min(record.seq().len());
                let end = (start + terminal_width.saturating_sub(1)).min(record.seq().len());
                &record.seq()[start..end]
            } else {
                &record.seq()
            };

            let visible_quality = if no_wrap && record.qual().len() > horizontal_offset {
                let quality_end =
                    (horizontal_offset + terminal_width.saturating_sub(1)).min(record.qual().len());
                &record.qual()[horizontal_offset..quality_end]
            } else {
                &record.qual()
            };

            let sequence_spans = if !self.show_quality {
                self.colorize_sequence(visible_sequence, visible_quality)
            } else {
                vec![Span::raw(
                    String::from_utf8_lossy(visible_sequence).to_string(),
                )]
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
            "↑/k: Up | ↓/j: Down | ←/→: Scroll | PgUp/PgDn: Page | S: Wrap | s: Stats | c: Colors | q: Quit"
        } else {
            "↑/k: Up | ↓/j: Down | PgUp/PgDn: Page | S: No-Wrap | s: Stats | c: Colors | q: Quit"
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
                    "🔄"
                } else if stats_lock.scanned_all {
                    "✅"
                } else {
                    "⏳"
                };

                let title = Line::from(Span::styled(
                    format!("FASTQ Statistics {}", status_indicator),
                    Style::default(),
                ));
                f.render_widget(title, title_area);

                // Basic stats block
                if let Some(basic_area) = main_areas.get(0).and_then(|areas| areas.get(0)) {
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
                if let Some(histogram_area) = main_areas.get(1).and_then(|areas| areas.get(0)) {
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
                if let Some(adapter_area) = main_areas.get(0).and_then(|areas| areas.get(1)) {
                    let adapter_stats_block = create_adapter_stats_display(&stats_lock);
                    f.render_widget(adapter_stats_block, *adapter_area);
                }

                // Position quality block
                if let Some(position_area) = main_areas.get(2).and_then(|areas| areas.get(0)) {
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

                // Adapter stats block
                if let Some(adapter_area) = main_areas.get(0).and_then(|areas| areas.get(1)) {
                    let adapter_stats_block = create_adapter_stats_display(&stats_lock);
                    f.render_widget(adapter_stats_block, *adapter_area);
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
                let status = Paragraph::new(format!("fqless - Help Screen (Press 'h' to exit)"))
                    .style(Style::default());

                f.render_widget(status, help_chunks[0]);

                // Full screen help content
                let help_content = vec![
                    Line::from(Span::styled("FQLESS - FastQ File Viewer", Style::default())),
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
                    Line::from("  c          - Toggle quality-based color coding"),
                    Line::from("  s          - Toggle statistics panel"),
                    Line::from("  h          - Toggle this help screen"),
                    Line::from(""),
                    Line::from("Other:"),
                    Line::from("  q          - Quit"),
                    Line::from(""),
                    Line::from("Color Coding:"),
                    Line::from("  High quality bases appear in green"),
                    Line::from("  Medium quality bases appear in yellow/orange"),
                    Line::from("  Low quality bases appear in red"),
                    Line::from(""),
                    Line::from("Statistics Panel:"),
                    Line::from("  Shows real-time analysis of FastQ file"),
                    Line::from("  Quality histogram shows distribution of quality scores"),
                    Line::from("  Position quality shows per-base quality across reads"),
                    Line::from("  GC content and N content percentages"),
                    Line::from(""),
                    Line::from("Scrolling:"),
                    Line::from("  In stats/help mode: ↑/↓ scroll line by line"),
                    Line::from("  PgUp/PgDn scroll by pages"),
                    Line::from("  Home goes to top"),
                    Line::from(""),
                    Line::from("File Support:"),
                    Line::from("  Supports both compressed (.gz) and uncompressed FastQ"),
                    Line::from("  Streaming mode for large files (doesn't load everything)"),
                    Line::from("  Indexed mode for fast random access in uncompressed files"),
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
                "File: {} | Records: {} | Pos: {} | {} | {} {}-{} ",
                filename,
                self.stats
                    .lock()
                    .unwrap()
                    .total_reads
                    .to_formatted_string(&Locale::en),
                current_position,
                wrap_status,
                self.phred_range.name(),
                self.phred_range.base_phred(),
                self.phred_range.top_phred()
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
                Paragraph::new(help_text).style(Style::default().fg(Color::DarkGray));

            f.render_widget(help, main_chunks[2]);
        })?;

        Ok(())
    }
}

impl Drop for TuiViewer {
    fn drop(&mut self) {
        self.stop_stats_worker();

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
