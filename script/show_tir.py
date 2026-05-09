#!/usr/bin/env python3
# Usage: python show_tir.py input.geve.fna
# Output: input.geve.tir.fna

import re
import sys
from pathlib import Path

REVCOMP_TABLE = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

def revcomp(seq):
    return seq.translate(REVCOMP_TABLE)[::-1]

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
    ml = re.search(r"tirL=(\d+)\.\.(\d+)", header)
    mr = re.search(r"tirR=(\d+)\.\.(\d+)", header)
    if not ml or not mr:
        return None
    return int(ml.group(1)), int(ml.group(2)), int(mr.group(1)), int(mr.group(2))

def get_geve_name(header):
    m = re.match(r">(\S+)", header)
    return m.group(1) if m else "UNKNOWN"

def write_fasta(fh, name, seq):
    fh.write(f">{name}\n")
    for i in range(0, len(seq), 80):
        fh.write(seq[i:i + 80] + "\n")

def main():
    if len(sys.argv) < 2:
        print("Usage: show_tir.py input.geve.fna")
        sys.exit(1)

    in_path  = Path(sys.argv[1])
    stem     = in_path.name.replace(".fna", "")
    out_path = in_path.with_name(stem + ".tir.fna")

    records   = parse_fna(in_path)
    n_written = 0
    n_skipped = 0

    with open(out_path, "w") as fh:
        for header, seq in records:
            tir = parse_tir(header)
            if tir is None:
                n_skipped += 1
                continue

            geve_start          = parse_int(header, "start")
            tl_s, tl_e, tr_s, tr_e = tir

            # Convert genome-absolute coords to local 0-based indices
            tir_l = seq[tl_s - geve_start : tl_e - geve_start + 1]
            tir_r = seq[tr_s - geve_start : tr_e - geve_start + 1]

            geve_name = get_geve_name(header)
            write_fasta(fh, f"{geve_name}_L", tir_l)
            write_fasta(fh, f"{geve_name}_R", tir_r)
            n_written += 1

    print(f"Done: {n_written} TIR pair(s) written, {n_skipped} no-TIR skipped")
    print(f"Output -> {out_path}")

if __name__ == "__main__":
    main()
