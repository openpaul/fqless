# fqless - FastQ Less

![Coverage](coverage-badge.svg) ![Coverage (no TUI)](coverage-no-tui-badge.svg)

fqless is a small less-like viewer for FastQ sequencing files. It allows the user to get a quick glance at their sequencing data without the need for a heavy GUI.

It displays the name and the color-coded sequence. It hides the quality line from the user, as this is machine information for machines and not for human beings.

It has no problem opening many gigabytes of FASTQ files, as it will not load everything but will stream the reads.
It also spawns workers to build a FastQC-like stats page (access via pressing `s`). As such, after opening a file you will see high CPU usage as reads are counted and adapters are detected.

fqless is released under the GPLv2 or any newer version of the GPL. It comes with no warranty. Use it as you wish.

Use it via

```sh
fqless file.fastq.gz
fqless file.fastq
```

Or pipe into it (although stats are then only shown for all loaded reads, not the whole file):

```sh
cat file.fastq | fqless
cat file.fastq.gz | fqless
```

See all options via `-h` flag.

## How to build

### On UNIX:

```
git clone https://github.com/openpaul/fqless
cd fqless
cargo build --release
```

## What it does

![a screenshot of fqless](https://raw.githubusercontent.com/openpaul/fqless/master/fqless.png)

- Opens the file, shows the sequence color coded, so one can decide if the run quality is nice and fits the expectations.
- Shows number of reads and some basic statistics
- Press `r` to cycle read orientation: as-is (`5'->3'`) or reverse complemented. Quality scores stay paired with their bases.
- Press `t` to show all six reading frames translated to amino acids (stops shown in red as `*`).
- Press `/` to search: matches are highlighted in yellow as you type. `n`/`N` jump to the next/previous match, and `x` toggles reverse-complement search (reads are unoriented, so by default a pattern also matches its reverse complement). Searches are exact and case-insensitive, and a background thread indexes the whole file (one bit per read) so you can jump to matches beyond what is loaded.

## Statistics (press 's')

fqless calculates some statistics in the background:

### Basic Statistics
- **Total reads**: Total number of reads in the file (counted quickly in a separate thread)
- **Processed**: Number of reads analyzed for detailed stats (will eventually be equal to total reads)
- **Avg/Min/Max length**: Read length statistics
- **GC content**: Percentage of G and C bases
- **N content**: Percentage of ambiguous (N) bases

### Quality Metrics
- **Per-Base Quality Score Distribution**: Shows the distribution of quality scores for individual bases across all reads. Quality scores are binned in groups of 5 (Q0-4, Q5-9, etc.).
- **Average Quality Per Read Distribution**: Shows the distribution of average quality scores per read. Each read contributes one value (its mean quality).
- **Average Quality by Position**: Quality scores averaged across all reads at each position in the read.

### Adapter Contamination
- **Contaminated reads**: Number and percentage of reads containing adapter sequences
- **Total adapter instances**: Total adapter sequences found (a single read may contain multiple adapters)
- **Top adapters detected**: Lists the most common adapter sequences found with their frequency and average position

Supported adapter types:
- Illumina (TruSeq Universal, TruSeq Index, Nextera)
- Oxford Nanopore (Rapid and Ligation adapters)
- PacBio (SMRTbell and Iso-Seq)
- PolyA tails

## What it does not
- Write to disk: Nothing is written to disk
- Fuzzy search: Only exact substring matches are supported, not wildcards or regular expressions.


# Bugs

This software has most certainly some bugs. 
I am testing a lot, but will probably not catch everything. 
So I am very glad if you submit issues. Best would be to supply a minimal working, or not working, fastq file, that causes the problem.

![fqless stat feature](https://raw.githubusercontent.com/openpaul/fqless/master/fqless_stats.png)


## LLM disclaimer

Claude Sonnet 4.5 was used intermittely during development.