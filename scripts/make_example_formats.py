#!/usr/bin/env python3
"""Generate the FASTA and SAM examples used to test fqless.

The FASTQ files in example/ are the source of truth. From them this script
writes:

* example/<name>.fasta          - the same reads as FASTA (no qualities)
* example/fake_reference.fasta  - one fake chromosome built from the reads,
                                  with 10 bp N spacers between them
* example/<name>.sam            - the reads "aligned" to that fake chromosome

Each read is placed at its own offset, so every alignment is a perfect
`<len>M` match; the last read of each file is stored reverse-complemented
(FLAG 16) so the reverse-strand handling is exercised. The first read also
gets a secondary alignment (FLAG 256), which fqless must skip.

Run from the repository root:

    python3 scripts/make_example_formats.py
"""

from pathlib import Path

EXAMPLE_DIR = Path(__file__).resolve().parent.parent / "example"
FASTQ_FILES = ["simple.fastq", "fqless.fastq"]
SPACER = "N" * 10
SPACER_LEN = len(SPACER)
REFERENCE_NAME = "chr1"

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def read_fastq(path: Path) -> list[tuple[str, str, str]]:
    """Return (header, sequence, quality) tuples from a 4-line FASTQ file."""
    lines = path.read_text().splitlines()
    records = []
    for i in range(0, len(lines), 4):
        header, seq, plus, qual = lines[i : i + 4]
        assert header.startswith("@"), f"{path}: bad header at line {i + 1}"
        assert plus.startswith("+"), f"{path}: missing '+' at line {i + 3}"
        assert len(seq) == len(qual), f"{path}: seq/qual length mismatch"
        records.append((header[1:], seq, qual))
    return records


def main() -> None:
    reference_parts: list[str] = []
    layouts: dict[str, list[tuple[str, str, str, int, bool]]] = {}

    for name in FASTQ_FILES:
        records = read_fastq(EXAMPLE_DIR / name)
        layout = []

        for index, (header, seq, qual) in enumerate(records):
            reverse = index == len(records) - 1
            stored = revcomp(seq) if reverse else seq
            offset = sum(len(part) for part in reference_parts) + SPACER_LEN * max(
                len(reference_parts) - 1, 0
            )
            reference_parts.append(stored)
            layout.append((header, stored, qual, offset, reverse))

        layouts[name] = layout

        fasta = "".join(f">{header}\n{seq}\n" for header, seq, _ in records)
        (EXAMPLE_DIR / name.replace(".fastq", ".fasta")).write_text(fasta)

    reference = SPACER.join(reference_parts)
    (EXAMPLE_DIR / "fake_reference.fasta").write_text(
        f">{REFERENCE_NAME}\n{reference}\n"
    )

    for name, layout in layouts.items():
        sam = [
            "@HD\tVN:1.6\tSO:coordinate",
            f"@SQ\tSN:{REFERENCE_NAME}\tLN:{len(reference)}",
            "@PG\tID:fqless-example\tPN:fqless-example\tVN:1.0",
        ]

        for header, stored, qual, offset, reverse in layout:
            qname = header.split()[0]
            flag = 16 if reverse else 0
            sam.append(
                "\t".join(
                    [
                        qname,
                        str(flag),
                        REFERENCE_NAME,
                        str(offset + 1),
                        "60",
                        f"{len(stored)}M",
                        "*",
                        "0",
                        "0",
                        stored,
                        qual,
                    ]
                )
            )

        # One secondary alignment per file: fqless must ignore it.
        header, stored, qual, offset, reverse = layout[0]
        qname = header.split()[0]
        sam.append(
            "\t".join(
                [
                    qname,
                    str(256 | (16 if reverse else 0)),
                    REFERENCE_NAME,
                    str(offset + 1),
                    "0",
                    f"{len(stored)}M",
                    "*",
                    "0",
                    "0",
                    stored,
                    qual,
                ]
            )
        )

        (EXAMPLE_DIR / name.replace(".fastq", ".sam")).write_text("\n".join(sam) + "\n")


if __name__ == "__main__":
    main()
