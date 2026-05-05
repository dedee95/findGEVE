#!/usr/bin/env python3
# Usage: python trim_tir.py input.geve.fna
# Output: input.geve.clean.fna

import re
import sys
from pathlib import Path

def parse_fna(path):
    records = []
    header = None
    seq_lines = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_lines)))
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            records.append((header, "".join(seq_lines)))
    return records

def parse_int(header, key):
    m = re.search(rf"{key}=(\d+)", header)
    return int(m.group(1)) if m else None

def parse_tir(header):
    m = re.search(r"tirL=(\d+)\.\.(\d+)", header)
    if not m:
        return None, None, None, None
    tl_s, tl_e = int(m.group(1)), int(m.group(2))
    m = re.search(r"tirR=(\d+)\.\.(\d+)", header)
    if not m:
        return None, None, None, None
    tr_s, tr_e = int(m.group(1)), int(m.group(2))
    return tl_s, tl_e, tr_s, tr_e

def write_fasta(fh, header, seq):
    fh.write(header + "\n")
    for i in range(0, len(seq), 80):
        fh.write(seq[i:i + 80] + "\n")

def trim_tir(seq, geve_start, tl_s, tl_e, tr_s, tr_e):
    # Convert genome-absolute TIR coords to local 0-based indices.
    # geve_start is 1-based; seq[0] == genome position geve_start.
    inner_start = tl_e - geve_start + 1   # first base after left TIR
    inner_end   = tr_s - geve_start        # one past last base before right TIR
    return seq[inner_start:inner_end]

def update_header(header, new_seq, geve_start):
    new_end    = geve_start + len(new_seq) - 1
    new_length = len(new_seq)
    header = re.sub(r"end=\d+",    f"end={new_end}",       header)
    header = re.sub(r"length=\d+", f"length={new_length}", header)
    header = re.sub(r" tirL=\S+",  "",                     header)
    header = re.sub(r" tirR=\S+",  "",                     header)
    header = re.sub(r" tir_id=\S+","",                     header)
    return header

def main():
    if len(sys.argv) < 2:
        print("Usage: trim_tir.py input.geve.fna")
        sys.exit(1)

    in_path  = Path(sys.argv[1])
    out_path = in_path.with_suffix("").with_suffix(".clean.fna")
    if in_path.suffix == ".fna":
        stem = in_path.stem
        if stem.endswith(".geve"):
            out_path = in_path.with_name(stem + ".clean.fna")
        else:
            out_path = in_path.with_name(stem + ".clean.fna")

    records = parse_fna(in_path)
    n_trimmed = 0
    n_kept    = 0

    with open(out_path, "w") as fh:
        for header, seq in records:
            geve_start = parse_int(header, "start")
            tl_s, tl_e, tr_s, tr_e = parse_tir(header)

            if tl_s is not None and geve_start is not None:
                trimmed = trim_tir(seq, geve_start, tl_s, tl_e, tr_s, tr_e)
                new_header = update_header(header, trimmed, geve_start)
                write_fasta(fh, new_header, trimmed)
                n_trimmed += 1
            else:
                write_fasta(fh, header, seq)
                n_kept += 1

    print(f"Done: {n_trimmed} TIR trimmed, {n_kept} no-TIR (kept as-is)")
    print(f"Output -> {out_path}")

if __name__ == "__main__":
    main()
