#!/usr/bin/env python3
"""
run_pyrodigal_v2.py - Predict prokaryotic ORFs in complete genome assemblies with Pyrodigal-GV.
Author: Dede Kurniawan (dedekurniawan@genomics.cn)
"""

from __future__ import annotations

import argparse
import gzip
import importlib
import logging
import sys
import tempfile
import time
from bisect import bisect_right
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


HELP_TEXT = """\
run_pyrodigal_v2.py - Predict prokaryotic ORFs in a whole genome with Pyrodigal-GV.

Usage: run_pyrodigal_v2.py --prefix <prefix> <genome.fa> [OPTIONS]

Mandatory:
  --prefix             Prefix used for ORF IDs (prefix_ORF1, prefix_ORF2, ...)
  genome               Input genome assembly FASTA (gzip is acceptable)

Optionals:
  -o, --outdir         Output directory                              [default: .]
  -t, --threads        CPU workers for ORF prediction                [default: 4]
  -g, --gff            Host eukaryotic annotation in GFF/GFF3;
                       predicted ORFs overlapping gene/transcript/CDS/exon
                       spans are removed (gzip is acceptable)
  -h, --help           Show this help and exit

Outputs:
  <prefix>.orf.pep     Predicted protein sequences
  <prefix>.orf.bed     Four columns: contig, start, end, ORF ID
                       (standard BED coordinates: 0-based, end-exclusive)
"""

USAGE_TEXT = "Usage: run_pyrodigal_v2.py --prefix <prefix> <genome.fa> [OPTIONS]\n"

_LOG = logging.getLogger("run_pyrodigal_v2")


@dataclass
class Orf:
    contig_index: int
    contig: str
    start: int
    end: int
    strand: int
    protein: str
    orf_id: str = ""


@dataclass
class HostIntervalIndex:
    intervals: List[Tuple[int, int]]
    starts: List[int]


_GFF_GENE_LIKE_FEATURES = {
    "gene",
    "mrna", "transcript", "primary_transcript",
    "lncrna", "ncrna", "rrna", "trna", "snorna", "snrna", "mirna",
    "pseudogene", "pseudogenic_transcript",
}
_GFF_PART_FEATURES = {"cds", "exon"}
_GFF_IGNORED_FEATURES = {
    "intron", "start_codon", "stop_codon",
    "five_prime_utr", "three_prime_utr", "5utr", "3utr", "utr",
}


def setup_logging() -> None:
    """Configure concise console logging without creating a log file."""
    _LOG.setLevel(logging.INFO)
    _LOG.handlers.clear()
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    ))
    _LOG.addHandler(handler)
    _LOG.propagate = False


def _is_gzip(path: Path) -> bool:
    """Detect gzip input by its magic bytes, independent of file extension."""
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def _open_text(path: Path):
    """Open a plain-text or gzip-compressed file for reading."""
    if _is_gzip(path):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def read_fasta(genome_path: Path) -> List[Tuple[int, str, str]]:
    """Read FASTA records using only the first whitespace-delimited header word.

    Header parsing stops only at whitespace. Underscores, periods, pipes,
    hyphens, colons, and other non-whitespace symbols remain part of the contig
    ID. Records are returned in input order as (record_index, contig, sequence).

    Some assemblies contain repeated first-word IDs with different descriptions.
    Version 2 accepts these records instead of aborting. They remain distinct
    internally through ``record_index``, while the PEP/BED contig field retains
    the exact first-word ID requested by the user.
    """
    records: List[Tuple[int, str, str]] = []
    id_counts: Dict[str, int] = defaultdict(int)
    current_name: Optional[str] = None
    sequence_parts: List[str] = []

    def flush_record() -> None:
        nonlocal current_name, sequence_parts
        if current_name is None:
            return
        sequence = "".join(sequence_parts).upper()
        if not sequence:
            raise ValueError(f"FASTA record '{current_name}' has an empty sequence")
        records.append((len(records), current_name, sequence))
        sequence_parts = []

    with _open_text(genome_path) as fh:
        for line_number, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush_record()
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header at line {line_number}")

                # Split on whitespace only. Symbols such as underscore, pipe,
                # period, colon, slash, and hyphen are intentionally retained.
                current_name = header.split(maxsplit=1)[0]
                id_counts[current_name] += 1
                continue

            if current_name is None:
                raise ValueError(
                    f"FASTA sequence encountered before the first header at line {line_number}"
                )
            sequence_parts.append("".join(line.split()))

    flush_record()
    if not records:
        raise ValueError("No FASTA records were found")

    duplicated = {name: count for name, count in id_counts.items() if count > 1}
    if duplicated:
        examples = ", ".join(
            f"{name!r} ({count} records)"
            for name, count in list(duplicated.items())[:5]
        )
        if len(duplicated) > 5:
            examples += f", ... and {len(duplicated) - 5:,} more"
        _LOG.warning(
            f"Repeated first-word FASTA contig ID(s) detected and accepted: {examples}. "
            "Records are processed separately, and BED contig names retain the first word."
        )

    return records


def _merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge overlapping or directly adjacent 1-based inclusive intervals."""
    if not intervals:
        return []
    normalized = sorted((min(start, end), max(start, end)) for start, end in intervals)
    merged: List[Tuple[int, int]] = []
    for start, end in normalized:
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def parse_host_gff_intervals(gff_path: Path) -> Dict[str, HostIntervalIndex]:
    """Parse host annotation spans from a plain or gzip-compressed GFF/GFF3.

    Gene/transcript-like rows and CDS/exon rows are accepted. Coordinates are
    stored as 1-based inclusive intervals, matching Pyrodigal-GV Gene.begin and
    Gene.end. Any predicted ORF touching one of these intervals is removed.
    """
    intervals_by_contig: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    n_rows = 0
    n_used = 0

    with _open_text(gff_path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            columns = line.split("\t")
            if len(columns) < 5:
                continue
            n_rows += 1
            contig = columns[0].split(None, 1)[0]
            feature = columns[2].strip().lower() if len(columns) > 2 else ""
            try:
                start = int(columns[3])
                end = int(columns[4])
            except (TypeError, ValueError):
                continue
            if start <= 0 or end <= 0:
                continue

            if feature in _GFF_IGNORED_FEATURES:
                continue

            use_feature = feature in _GFF_GENE_LIKE_FEATURES or feature in _GFF_PART_FEATURES
            if not use_feature:
                use_feature = any(
                    token in feature
                    for token in ("gene", "transcript", "mrna", "cds", "exon")
                )
            if not use_feature:
                continue

            intervals_by_contig[contig].append((min(start, end), max(start, end)))
            n_used += 1

    indexes: Dict[str, HostIntervalIndex] = {}
    for contig, intervals in intervals_by_contig.items():
        merged = _merge_intervals(intervals)
        indexes[contig] = HostIntervalIndex(
            intervals=merged,
            starts=[start for start, _end in merged],
        )

    n_intervals = sum(len(index.intervals) for index in indexes.values())
    _LOG.info(
        f"Host GFF mask: parsed {n_rows:,} feature row(s); used {n_used:,}; "
        f"collapsed to {n_intervals:,} interval(s) on {len(indexes):,} contig(s)"
    )
    if not indexes:
        _LOG.warning(
            "Host GFF mask contained no usable gene/transcript/CDS/exon intervals; "
            "no predicted ORFs will be removed"
        )
    return indexes


def _overlaps_host_interval(index: Optional[HostIntervalIndex], start: int, end: int) -> bool:
    """Return True when a 1-based inclusive ORF overlaps any host interval."""
    if index is None or not index.intervals:
        return False
    position = bisect_right(index.starts, end) - 1
    return position >= 0 and index.intervals[position][1] >= start


def _load_pyrodigal_gv():
    """Import Pyrodigal-GV lazily so --help works even before installation."""
    try:
        return importlib.import_module("pyrodigal_gv")
    except ImportError as exc:
        raise RuntimeError(
            "pyrodigal-gv is not installed. Install it with: "
            "conda install -c bioconda pyrodigal-gv  (or pip install pyrodigal-gv)"
        ) from exc


def _predict_orfs_on_contig(
    args: Tuple[int, str, str]
) -> Tuple[int, str, List[Tuple[int, int, int, str]], Optional[str]]:
    """Run Pyrodigal-GV on one complete contig."""
    contig_index, contig_name, sequence = args
    try:
        pyrodigal_gv = _load_pyrodigal_gv()
        finder = pyrodigal_gv.ViralGeneFinder(meta=True)
        genes = finder.find_genes(sequence.encode("ascii"))
    except Exception as exc:
        return contig_index, contig_name, [], f"pyrodigal-gv failed on {contig_name}: {exc}"

    records: List[Tuple[int, int, int, str]] = []
    for gene in genes:
        try:
            start = int(gene.begin)
            end = int(gene.end)
            strand = int(gene.strand)
            protein = str(gene.translate()).rstrip("*")
        except Exception as exc:
            return (
                contig_index,
                contig_name,
                [],
                f"failed to read a predicted ORF on {contig_name}: {exc}",
            )
        if not protein:
            continue
        records.append((min(start, end), max(start, end), strand, protein))

    records.sort(key=lambda row: (row[0], row[1], row[2]))
    return contig_index, contig_name, records, None


def predict_orfs(
    fasta_records: List[Tuple[int, str, str]],
    threads: int,
) -> List[Orf]:
    """Run Pyrodigal-GV in meta mode on every complete contig."""
    executor = None
    if threads > 1 and len(fasta_records) > 1:
        executor = ProcessPoolExecutor(max_workers=threads)
        results_iter = executor.map(_predict_orfs_on_contig, fasta_records, chunksize=1)
    else:
        results_iter = (_predict_orfs_on_contig(record) for record in fasta_records)

    predicted: List[Orf] = []
    errors: List[str] = []
    try:
        for contig_index, contig_name, records, error in results_iter:
            if error is not None:
                errors.append(error)
                continue
            for start, end, strand, protein in records:
                predicted.append(Orf(
                    contig_index=contig_index,
                    contig=contig_name,
                    start=start,
                    end=end,
                    strand=strand,
                    protein=protein,
                ))
    finally:
        if executor is not None:
            executor.shutdown(wait=True)

    if errors:
        for error in errors[:10]:
            _LOG.error(error)
        if len(errors) > 10:
            _LOG.error(f"... and {len(errors) - 10:,} additional contig error(s)")
        raise RuntimeError(
            f"Pyrodigal-GV failed on {len(errors):,} of {len(fasta_records):,} contig(s)"
        )

    predicted.sort(key=lambda orf: (orf.contig_index, orf.start, orf.end, orf.strand))
    _LOG.info(
        f"ORF prediction: {len(predicted):,} ORF(s) on {len(fasta_records):,} contig(s)"
    )
    return predicted


def filter_orfs_by_host_gff(
    orfs: List[Orf],
    host_indexes: Dict[str, HostIntervalIndex],
) -> Tuple[List[Orf], int]:
    """Remove every predicted ORF that overlaps a host annotation interval."""
    retained: List[Orf] = []
    removed = 0
    for orf in orfs:
        if _overlaps_host_interval(host_indexes.get(orf.contig), orf.start, orf.end):
            removed += 1
        else:
            retained.append(orf)
    _LOG.info(
        f"Host GFF mask: removed {removed:,} overlapping ORF(s); "
        f"{len(retained):,} ORF(s) retained"
    )
    return retained, removed


def assign_orf_ids(orfs: List[Orf], prefix: str) -> None:
    """Assign deterministic, contiguous genome-wide ORF IDs."""
    for number, orf in enumerate(orfs, start=1):
        orf.orf_id = f"{prefix}_ORF{number}"


def _write_wrapped_sequence(fh, sequence: str, width: int = 80) -> None:
    for start in range(0, len(sequence), width):
        fh.write(sequence[start:start + width] + "\n")


def write_outputs(orfs: List[Orf], pep_path: Path, bed_path: Path) -> None:
    """Write protein FASTA and four-column standard BED files atomically."""
    pep_tmp: Optional[Path] = None
    bed_tmp: Optional[Path] = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=pep_path.parent,
            prefix=f".{pep_path.name}.", suffix=".tmp", delete=False,
        ) as pep_fh:
            pep_tmp = Path(pep_fh.name)
            for orf in orfs:
                pep_fh.write(f">{orf.orf_id}\n")
                _write_wrapped_sequence(pep_fh, orf.protein)

        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=bed_path.parent,
            prefix=f".{bed_path.name}.", suffix=".tmp", delete=False,
        ) as bed_fh:
            bed_tmp = Path(bed_fh.name)
            for orf in orfs:
                # Convert Pyrodigal/GFF 1-based inclusive coordinates to BED.
                bed_fh.write(
                    f"{orf.contig}\t{orf.start - 1}\t{orf.end}\t{orf.orf_id}\n"
                )

        pep_tmp.replace(pep_path)
        bed_tmp.replace(bed_path)
    except Exception:
        if pep_tmp is not None:
            pep_tmp.unlink(missing_ok=True)
        if bed_tmp is not None:
            bed_tmp.unlink(missing_ok=True)
        raise


class _Parser(argparse.ArgumentParser):
    def format_help(self) -> str:
        return HELP_TEXT

    def format_usage(self) -> str:
        return USAGE_TEXT


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = _Parser(prog="run_pyrodigal_v2.py", add_help=False)
    parser.add_argument("-h", "--help", action="help", help="Show this help and exit")
    parser.add_argument("genome", type=Path)
    parser.add_argument("--prefix", type=str, required=True)
    parser.add_argument("-o", "--outdir", type=Path, default=Path("."))
    parser.add_argument("-t", "--threads", type=int, default=4)
    parser.add_argument(
        "-g", "--gff", type=Path, default=None,
        help="Optional host eukaryotic GFF/GFF3 mask",
    )
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    setup_logging()
    started = time.time()

    if not args.genome.is_file():
        _LOG.error(f"Genome file not found: {args.genome}")
        return 2
    if args.gff is not None and not args.gff.is_file():
        _LOG.error(f"GFF annotation file not found: {args.gff}")
        return 2
    if args.threads < 1:
        _LOG.error("--threads must be at least 1")
        return 2
    if not args.prefix or any(char.isspace() for char in args.prefix):
        _LOG.error("--prefix must be non-empty and must not contain whitespace")
        return 2

    try:
        _load_pyrodigal_gv()
    except RuntimeError as exc:
        _LOG.error(str(exc))
        return 2

    try:
        args.outdir.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        _LOG.error(f"Cannot create output directory {args.outdir}: {exc}")
        return 2

    pep_path = args.outdir / f"{args.prefix}.orf.pep"
    bed_path = args.outdir / f"{args.prefix}.orf.bed"

    _LOG.info(
        f"run_pyrodigal_v2 started | prefix='{args.prefix}' | "
        f"threads={args.threads} | genome={args.genome}"
    )
    if args.gff is not None:
        _LOG.info(f"Host GFF mask enabled: {args.gff}")

    try:
        fasta_records = read_fasta(args.genome)
        _LOG.info(f"Input FASTA: {len(fasta_records):,} contig(s)")
        orfs = predict_orfs(fasta_records, args.threads)
        if args.gff is not None:
            host_indexes = parse_host_gff_intervals(args.gff)
            orfs, _removed = filter_orfs_by_host_gff(orfs, host_indexes)
        assign_orf_ids(orfs, args.prefix)
        write_outputs(orfs, pep_path, bed_path)
    except (OSError, UnicodeError, ValueError, RuntimeError) as exc:
        _LOG.error(str(exc))
        return 2

    _LOG.info(f"Proteins -> {pep_path}")
    _LOG.info(f"BED      -> {bed_path}")
    _LOG.info(
        f"run_pyrodigal_v2 completed: {len(orfs):,} retained ORF(s) "
        f"in {time.time() - started:.1f} s"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
