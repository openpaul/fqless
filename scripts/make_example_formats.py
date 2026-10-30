#!/usr/bin/env python3
"""Generate the FASTA and SAM examples used to test fqless.

The FASTQ files in example/ are the source of truth. From them this script
writes:

* example/<name>.fasta          - the same reads as FASTA (no qualities)
* example/fake_reference.fasta  - one fake chromosome built from the reads,
                                  with 10 bp N spacers between them
* example/<name>.sam            - the reads "aligned" to that fake chromosome
* example/<name>.bam            - the same alignments, BGZF-compressed

Each read is placed at its own offset, so every alignment is a perfect
`<len>M` match; the last read of each file is stored reverse-complemented
(FLAG 16) so the reverse-strand handling is exercised. The first read also
gets a secondary alignment (FLAG 256), which fqless must skip.

Run from the repository root:

    python3 scripts/make_example_formats.py
"""

import re
import struct
import zlib
from pathlib import Path

EXAMPLE_DIR = Path(__file__).resolve().parent.parent / "example"
FASTQ_FILES = ["simple.fastq", "fqless.fastq"]
SPACER = "N" * 10
SPACER_LEN = len(SPACER)
REFERENCE_NAME = "chr1"

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# BAM packs bases four to a nibble, using this table from the specification.
SEQ_CODES = {base: code for code, base in enumerate("=ACMGRSVTWYHKDBN")}
CIGAR_OPS = "MIDNSHP=XB"
# Canonical 28-byte BGZF end-of-file marker.
BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def reg2bin(beg: int, end: int) -> int:
    """Binning index used by BAM, from the SAM specification."""
    end -= 1
    if beg >> 14 == end >> 14:
        return 4681 + (beg >> 14)
    if beg >> 17 == end >> 17:
        return 585 + (beg >> 17)
    if beg >> 20 == end >> 20:
        return 73 + (beg >> 20)
    if beg >> 23 == end >> 23:
        return 9 + (beg >> 23)
    if beg >> 26 == end >> 26:
        return 1 + (beg >> 26)
    return 0


def bgzf_compress(data: bytes) -> bytes:
    """Compress data as BGZF, the block-gzip variant BAM uses."""
    out = bytearray()

    for start in range(0, len(data), 60000):
        block = data[start : start + 60000]
        compressor = zlib.compressobj(6, zlib.DEFLATED, -15)
        payload = compressor.compress(block) + compressor.flush()
        block_size = len(payload) + 26
        out += struct.pack(
            "<BBBBIBBHBBHH",
            31,
            139,
            8,
            4,
            0,
            0,
            255,
            6,
            ord("B"),
            ord("C"),
            2,
            block_size - 1,
        )
        out += payload
        out += struct.pack("<II", zlib.crc32(block) & 0xFFFFFFFF, len(block))

    out += BGZF_EOF
    return bytes(out)


def encode_sequence(seq: str) -> bytes:
    codes = [SEQ_CODES[base] for base in seq]
    if len(codes) % 2:
        codes.append(0)
    return bytes((codes[i] << 4) | codes[i + 1] for i in range(0, len(codes), 2))


def parse_cigar(cigar: str) -> list[tuple[int, int]]:
    if cigar == "*":
        return []
    return [(int(n), CIGAR_OPS.index(op)) for n, op in re.findall(r"(\d+)([MIDNSHP=XB])", cigar)]


def write_bam(path: Path, sam_text: str) -> None:
    """Encode the generated SAM text as BGZF-compressed BAM."""
    header_lines = [line for line in sam_text.splitlines() if line.startswith("@")]
    body_lines = [line for line in sam_text.splitlines() if not line.startswith("@")]

    references = []
    for line in header_lines:
        if line.startswith("@SQ"):
            fields = dict(field.split(":", 1) for field in line.split("\t")[1:])
            references.append((fields["SN"], int(fields["LN"])))
    reference_ids = {name: i for i, (name, _) in enumerate(references)}

    text = ("\n".join(header_lines) + "\n").encode()
    out = bytearray(b"BAM\x01")
    out += struct.pack("<i", len(text) + 1) + text + b"\x00"
    out += struct.pack("<i", len(references))
    for name, length in references:
        encoded_name = name.encode() + b"\x00"
        out += struct.pack("<i", len(encoded_name)) + encoded_name
        out += struct.pack("<i", length)

    for line in body_lines:
        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.split(
            "\t"
        )[:11]
        ref_id = reference_ids.get(rname, -1)
        pos0 = int(pos) - 1 if pos not in ("0", "*") else -1
        ops = parse_cigar(cigar)
        span = sum(n for n, op in ops if op in (0, 2, 3, 7, 8))
        bin_ = reg2bin(pos0, pos0 + max(span, 1)) if pos0 >= 0 else 4680
        next_ref = reference_ids.get(rnext, -1)
        next_pos = int(pnext) - 1 if pnext not in ("0", "*") else -1
        l_seq = 0 if seq == "*" else len(seq)
        qual_bytes = b"\xff" * l_seq if qual == "*" else bytes(ord(c) - 33 for c in qual)
        read_name = qname.encode() + b"\x00"

        block = struct.pack(
            "<iiBBHHHiiii",
            ref_id,
            pos0,
            len(read_name),
            int(mapq),
            bin_,
            len(ops),
            int(flag),
            l_seq,
            next_ref,
            next_pos,
            int(tlen),
        )
        block += read_name
        block += b"".join(struct.pack("<I", (n << 4) | op) for n, op in ops)
        block += encode_sequence(seq if seq != "*" else "")
        block += qual_bytes

        out += struct.pack("<i", len(block)) + block

    path.write_bytes(bgzf_compress(bytes(out)))


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
            stored_qual = qual[::-1] if reverse else qual
            offset = sum(len(part) for part in reference_parts) + SPACER_LEN * max(
                len(reference_parts) - 1, 0
            )
            reference_parts.append(stored)
            layout.append((header, stored, stored_qual, offset, reverse))

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

        sam_text = "\n".join(sam) + "\n"
        sam_path = EXAMPLE_DIR / name.replace(".fastq", ".sam")
        sam_path.write_text(sam_text)
        write_bam(sam_path.with_suffix(".bam"), sam_text)


if __name__ == "__main__":
    main()
