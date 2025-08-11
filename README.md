# fqless - FastQ Less

fqless is intended to be a small less-like viewer for FastQ sequencing files. Allowing the user to get a quick glance at their sequencing data without the need of a heavy GUI/Java-bloatware.

It displays the name and the color coded sequence. It hides the quality line from the user, as this is machine information for machines and not for human beings.

It has no problem opening many gigabytes of fastq files, as it will not load everything, but will be streaming the reads.
It also spawns a worker to build a fastqc like stats page (access via pressing `s`).

fqless is released under the GPLv2 or any newer version of the GPL. It comes with no warranty. Use it as you wish.

Use it via

```sh
fqless file.fastq.gz
fqless file.fastq
```

Or pipe into it (although stats are then not supported and only shown for all loaded reads, not the whole file):

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
  It tries to guess the quality encoding of the file but you can also change that live, if it guessed the wrong one.
- Move by arrow keys. Change quality encoding by right and left arrow.

## What it does not

## Features

- Parse FastQ [done]
- read plain text as well as fastq.gz [done]
- Detect encoding [done]
- Color code sequence [done]
- allow switching between encodings [done]
- allow STDIN [done]
- compute fastq quality metrics

# Bugs

This software has most certainly some heavy bugs. I am testing a lot, but will probably not catch everything. So I am very glad if you submit issues. Best would be to supply a minimal working, or not working, fastq file, that causes the problem.
