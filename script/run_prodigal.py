#!/usr/bin/env python3
"""
run_prodigal.py - Predict Pyrodigal-GV ORFs outside identified GEVE regions.
Author: Dede Kurniawan (dedekurniawan@genomics.cn)
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import re
import sys
import time
from bisect import bisect_right
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from datetime import datetime
from importlib.metadata import PackageNotFoundError, version as _pkg_version
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

try:
    import pyrodigal_gv
except ImportError:  # handled explicitly in main so --help still works
    pyrodigal_gv = None


HELP_TEXT = """\
run_prodigal.py - Predict Pyrodigal-GV ORFs outside identified GEVE regions.

Usage: run_prodigal.py --prefix <prefix> <genome.fa> <summary.tsv> [OPTIONS]

Mandatory:
  --prefix             Prefix for predicted ORF IDs and output file names
  genome               Input whole-genome assembly FASTA (gzip is acceptable)
  summary              findGEVE summary.tsv containing contig_id, start, and end

Optionals:
  -o, --outdir         Output directory                    [default: ./Prodigal_<YYYYMMDD>]
  -t, --threads        Parallel contig workers             [default: 4]
  -h, --help           Show this help and exit

Outputs:
  <prefix>.nonGEVE.pep   Predicted protein sequences
  <prefix>.nonGEVE.gff3  Predicted CDS coordinates
  run.log                Detailed execution log
"""

USAGE_TEXT = (
    "Usage: run_prodigal.py --prefix <prefix> "
    "<genome.fa> <summary.tsv> [OPTIONS]\n"
)

_LOG = logging.getLogger("run_prodigal")
OUTPUT = 25
logging.addLevelName(OUTPUT, "OUTPUT")


def _output(self, message, *args, **kwargs):
    if self.isEnabledFor(OUTPUT):
        self._log(OUTPUT, message, args, **kwargs)


logging.Logger.output = _output

_NATKEY_RE = re.compile(r"(\d+)")
_PREFIX_RE = re.compile(r"^[^\s/\\;=]+$")


def _natural_key(value: str):
    return [
        int(token) if token.isdigit() else token.lower()
        for token in _NATKEY_RE.split(str(value))
    ]


def setup_logging(log_path: Optional[Path] = None) -> None:
    """Configure terminal and file logging in the same style as findGEVE.py."""
    _LOG.handlers.clear()
    _LOG.setLevel(logging.DEBUG)
    fmt = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(fmt)
    sh.setLevel(logging.INFO)
    _LOG.addHandler(sh)

    if log_path is not None:
        fh = logging.FileHandler(log_path, mode="w")
        fh.setFormatter(fmt)
        fh.setLevel(logging.DEBUG)
        _LOG.addHandler(fh)


@dataclass(frozen=True)
class Interval:
    """A 1-based inclusive genomic interval."""

    start: int
    end: int


@dataclass
class OrfRecord:
    """A retained ORF predicted outside all GEVE intervals."""

    contig: str
    start: int
    end: int
    strand: int
    protein: str
    partial_begin: bool = False
    partial_end: bool = False
    translation_table: Optional[int] = None
    orf_id: str = ""


@dataclass
class PredictionResult:
    contig: str
    orfs: List[OrfRecord]
    n_raw: int
    n_overlap_removed: int
    error: Optional[str] = None


class _Parser(argparse.ArgumentParser):
    def format_help(self) -> str:
        return HELP_TEXT

    def format_usage(self) -> str:
        return USAGE_TEXT


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = _Parser(prog="run_prodigal.py", add_help=False)
    parser.add_argument("-h", "--help", action="help", help="Show this help and exit")
    parser.add_argument("genome", type=Path)
    parser.add_argument("summary", type=Path)
    parser.add_argument("--prefix", type=str, required=True)
    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        default=Path(f"Prodigal_{datetime.now().strftime('%Y%m%d')}"),
    )
    parser.add_argument("-t", "--threads", type=int, default=4)
    return parser.parse_args(argv)


def _open_text(path: Path):
    """Open plain-text or gzip-compressed input based on the filename suffix."""
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    """Read a FASTA file and return unique (record_id, sequence) pairs.

    Record identifiers use the first whitespace-delimited token in each header,
    matching the contig naming convention used by findGEVE.py and pyfastx.
    """
    records: List[Tuple[str, str]] = []
    seen: set[str] = set()
    current_name: Optional[str] = None
    sequence_parts: List[str] = []

    def flush_record() -> None:
        nonlocal current_name, sequence_parts
        if current_name is None:
            return
        sequence = "".join(sequence_parts).replace(" ", "").replace("\t", "").upper()
        if not sequence:
            raise ValueError(f"FASTA record '{current_name}' has an empty sequence")
        records.append((current_name, sequence))
        sequence_parts = []

    with _open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush_record()
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header at line {line_number}")
                current_name = header.split(None, 1)[0]
                if current_name in seen:
                    raise ValueError(f"Duplicate FASTA record ID: {current_name}")
                seen.add(current_name)
                continue
            if current_name is None:
                raise ValueError(
                    f"Sequence data found before the first FASTA header at line {line_number}"
                )
            sequence_parts.append("".join(line.split()))

    flush_record()
    if not records:
        raise ValueError("No FASTA records were found")
    return records


def _merge_intervals(intervals: Iterable[Interval]) -> List[Interval]:
    """Merge overlapping or directly adjacent intervals without adding flanks."""
    ordered = sorted(intervals, key=lambda interval: (interval.start, interval.end))
    merged: List[Interval] = []
    for interval in ordered:
        if merged and interval.start <= merged[-1].end + 1:
            previous = merged[-1]
            merged[-1] = Interval(previous.start, max(previous.end, interval.end))
        else:
            merged.append(interval)
    return merged


def read_geve_intervals(summary_path: Path) -> Tuple[Dict[str, List[Interval]], int]:
    """Read exact GEVE coordinates from a findGEVE summary.tsv file."""
    by_contig: Dict[str, List[Interval]] = {}
    n_rows = 0

    with open(summary_path, "rt", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("summary.tsv is empty or has no header")

        fieldnames = [name.strip() if name is not None else "" for name in reader.fieldnames]
        reader.fieldnames = fieldnames
        required = {"contig_id", "start", "end"}
        missing = sorted(required.difference(fieldnames))
        if missing:
            raise ValueError(
                "summary.tsv is missing required column(s): " + ", ".join(missing)
            )

        for line_number, row in enumerate(reader, start=2):
            if not row or not any((value or "").strip() for value in row.values()):
                continue
            contig = (row.get("contig_id") or "").strip()
            start_text = (row.get("start") or "").strip()
            end_text = (row.get("end") or "").strip()
            if not contig:
                raise ValueError(f"Missing contig_id in summary.tsv line {line_number}")
            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid start/end coordinate in summary.tsv line {line_number}: "
                    f"start='{start_text}', end='{end_text}'"
                ) from exc
            if start < 1 or end < 1:
                raise ValueError(
                    f"Coordinates must be positive in summary.tsv line {line_number}: "
                    f"{contig}:{start}-{end}"
                )
            if start > end:
                raise ValueError(
                    f"start is greater than end in summary.tsv line {line_number}: "
                    f"{contig}:{start}-{end}"
                )
            by_contig.setdefault(contig, []).append(Interval(start, end))
            n_rows += 1

    if n_rows == 0:
        raise ValueError("summary.tsv contains no GEVE coordinate rows")

    return {
        contig: _merge_intervals(intervals)
        for contig, intervals in by_contig.items()
    }, n_rows


def validate_intervals(
    records: Sequence[Tuple[str, str]],
    intervals_by_contig: Dict[str, List[Interval]],
) -> Tuple[Dict[str, int], int]:
    """Validate GEVE contigs and coordinates against the input genome."""
    contig_lengths = {name: len(sequence) for name, sequence in records}
    missing_contigs = sorted(
        set(intervals_by_contig).difference(contig_lengths), key=_natural_key
    )
    if missing_contigs:
        preview = ", ".join(missing_contigs[:10])
        suffix = "" if len(missing_contigs) <= 10 else f" ... (+{len(missing_contigs) - 10})"
        raise ValueError(
            f"{len(missing_contigs)} GEVE contig(s) from summary.tsv are absent from "
            f"the genome FASTA: {preview}{suffix}"
        )

    masked_bp = 0
    for contig, intervals in intervals_by_contig.items():
        contig_length = contig_lengths[contig]
        for interval in intervals:
            if interval.end > contig_length:
                raise ValueError(
                    f"GEVE interval exceeds contig length: {contig}:{interval.start}-"
                    f"{interval.end}, contig length={contig_length}"
                )
            masked_bp += interval.end - interval.start + 1
    return contig_lengths, masked_bp


def _mask_sequence(sequence: str, intervals: Sequence[Interval]) -> bytes:
    """Replace exact 1-based inclusive GEVE spans with N while preserving length."""
    masked = bytearray(sequence.encode("ascii", errors="replace"))
    for interval in intervals:
        length = interval.end - interval.start + 1
        masked[interval.start - 1 : interval.end] = b"N" * length
    return bytes(masked)


def _overlaps_intervals(
    starts: Sequence[int],
    intervals: Sequence[Interval],
    start: int,
    end: int,
) -> bool:
    """Return True if a 1-based inclusive ORF overlaps any GEVE base."""
    if not intervals:
        return False
    index = bisect_right(starts, end) - 1
    return index >= 0 and intervals[index].end >= start


def _safe_translation_table(gene) -> Optional[int]:
    value = getattr(gene, "translation_table", None)
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _predict_orfs_on_contig(
    args: Tuple[str, str, List[Interval]],
) -> PredictionResult:
    """Run Pyrodigal-GV on one contig after exact GEVE masking."""
    contig, sequence, intervals = args
    try:
        masked_sequence = _mask_sequence(sequence, intervals)
        finder = pyrodigal_gv.ViralGeneFinder(
            meta=True,
            mask=True,
            min_mask=1,
        )
        genes = finder.find_genes(masked_sequence)
    except Exception as exc:
        return PredictionResult(
            contig=contig,
            orfs=[],
            n_raw=0,
            n_overlap_removed=0,
            error=f"pyrodigal-gv failed on {contig}: {exc}",
        )

    interval_starts = [interval.start for interval in intervals]
    retained: List[OrfRecord] = []
    n_raw = 0
    n_overlap_removed = 0

    for gene in genes:
        n_raw += 1
        try:
            start = int(gene.begin)
            end = int(gene.end)
            strand = int(gene.strand)
        except (AttributeError, TypeError, ValueError) as exc:
            return PredictionResult(
                contig=contig,
                orfs=[],
                n_raw=n_raw,
                n_overlap_removed=n_overlap_removed,
                error=f"Invalid pyrodigal-gv gene object on {contig}: {exc}",
            )

        start, end = min(start, end), max(start, end)
        if _overlaps_intervals(interval_starts, intervals, start, end):
            n_overlap_removed += 1
            continue

        try:
            protein = str(gene.translate()).rstrip("*")
        except Exception as exc:
            return PredictionResult(
                contig=contig,
                orfs=[],
                n_raw=n_raw,
                n_overlap_removed=n_overlap_removed,
                error=f"Protein translation failed on {contig}:{start}-{end}: {exc}",
            )
        if not protein:
            continue

        retained.append(
            OrfRecord(
                contig=contig,
                start=start,
                end=end,
                strand=strand,
                protein=protein,
                partial_begin=bool(getattr(gene, "partial_begin", False)),
                partial_end=bool(getattr(gene, "partial_end", False)),
                translation_table=_safe_translation_table(gene),
            )
        )

    return PredictionResult(
        contig=contig,
        orfs=retained,
        n_raw=n_raw,
        n_overlap_removed=n_overlap_removed,
        error=None,
    )


def predict_orfs(
    records: Sequence[Tuple[str, str]],
    intervals_by_contig: Dict[str, List[Interval]],
    threads: int,
) -> Tuple[List[OrfRecord], int, int]:
    """Predict ORFs on every genome contig, excluding exact GEVE intervals."""
    work_items = [
        (contig, sequence, intervals_by_contig.get(contig, []))
        for contig, sequence in records
    ]

    executor = None
    if threads > 1 and len(work_items) > 1:
        executor = ProcessPoolExecutor(max_workers=threads)
        results_iter: Iterator[PredictionResult] = executor.map(
            _predict_orfs_on_contig,
            work_items,
            chunksize=1,
        )
    else:
        results_iter = (_predict_orfs_on_contig(item) for item in work_items)

    all_orfs: List[OrfRecord] = []
    n_raw = 0
    n_overlap_removed = 0
    failed_contigs: List[str] = []
    try:
        for result in results_iter:
            if result.error is not None:
                failed_contigs.append(result.error)
                continue
            all_orfs.extend(result.orfs)
            n_raw += result.n_raw
            n_overlap_removed += result.n_overlap_removed
            _LOG.debug(
                f"Prediction {result.contig}: raw={result.n_raw:,}, "
                f"retained={len(result.orfs):,}, "
                f"overlap_removed={result.n_overlap_removed:,}"
            )
    finally:
        if executor is not None:
            executor.shutdown(wait=True)

    if failed_contigs:
        for message in failed_contigs:
            _LOG.error(message)
        raise RuntimeError(
            f"Pyrodigal-GV failed on {len(failed_contigs)} contig(s); no output was written"
        )

    all_orfs.sort(
        key=lambda orf: (
            _natural_key(orf.contig),
            orf.start,
            orf.end,
            -orf.strand,
        )
    )
    return all_orfs, n_raw, n_overlap_removed


def _partial_code(orf: OrfRecord) -> str:
    return f"{int(orf.partial_begin)}{int(orf.partial_end)}"


def _gff_escape(value: str) -> str:
    """Percent-encode characters that are unsafe in a GFF3 attribute value."""
    replacements = {
        "%": "%25",
        ";": "%3B",
        "=": "%3D",
        "&": "%26",
        ",": "%2C",
        "\t": "%09",
        "\n": "%0A",
        "\r": "%0D",
    }
    return "".join(replacements.get(character, character) for character in value)


def assign_orf_ids(orfs: Sequence[OrfRecord], prefix: str) -> None:
    for index, orf in enumerate(orfs, start=1):
        orf.orf_id = f"{prefix}_ORF_{index:07d}"


def write_protein_fasta(orfs: Sequence[OrfRecord], path: Path) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for orf in orfs:
            strand = "+" if orf.strand >= 0 else "-"
            fields = [
                f"contig={orf.contig}",
                f"start={orf.start}",
                f"end={orf.end}",
                f"strand={strand}",
                f"length={len(orf.protein)}",
                f"partial={_partial_code(orf)}",
            ]
            if orf.translation_table is not None:
                fields.append(f"transl_table={orf.translation_table}")
            handle.write(f">{orf.orf_id} {' '.join(fields)}\n")
            for offset in range(0, len(orf.protein), 80):
                handle.write(orf.protein[offset : offset + 80] + "\n")


def write_gff3(
    orfs: Sequence[OrfRecord],
    contig_lengths: Dict[str, int],
    path: Path,
) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        for contig in sorted(contig_lengths, key=_natural_key):
            handle.write(f"##sequence-region {contig} 1 {contig_lengths[contig]}\n")

        for orf in orfs:
            strand = "+" if orf.strand >= 0 else "-"
            attributes = [
                f"ID={_gff_escape(orf.orf_id)}",
                f"Name={_gff_escape(orf.orf_id)}",
                f"partial={_partial_code(orf)}",
                f"protein_length={len(orf.protein)}",
            ]
            if orf.translation_table is not None:
                attributes.append(f"transl_table={orf.translation_table}")
            handle.write(
                f"{orf.contig}\trun_prodigal\tCDS\t{orf.start}\t{orf.end}\t.\t"
                f"{strand}\t0\t{';'.join(attributes)}\n"
            )


def _get_tool_versions() -> Dict[str, str]:
    versions: Dict[str, str] = {"python": sys.version.split()[0]}
    for package, label in [
        ("pyrodigal-gv", "pyrodigal-gv"),
        ("pyrodigal", "pyrodigal"),
    ]:
        try:
            versions[label] = _pkg_version(package)
        except PackageNotFoundError:
            versions[label] = "unknown"
    return versions


def _format_elapsed(seconds: float) -> str:
    if seconds >= 3600:
        return f"{seconds / 3600:.2f} h"
    if seconds >= 60:
        return f"{seconds / 60:.2f} min"
    return f"{seconds:.1f} s"


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    if not args.prefix or not _PREFIX_RE.fullmatch(args.prefix):
        print(
            "ERROR: --prefix must be non-empty and cannot contain whitespace, "
            "path separators, ';', or '='.",
            file=sys.stderr,
        )
        return 2
    if args.threads < 1:
        print("ERROR: --threads must be at least 1.", file=sys.stderr)
        return 2

    args.outdir.mkdir(parents=True, exist_ok=True)
    log_path = args.outdir / "run.log"
    setup_logging(log_path)
    start_time = time.time()

    _LOG.info(
        f"run_prodigal started | prefix='{args.prefix}' | threads={args.threads} | "
        f"genome={args.genome} | summary={args.summary}"
    )
    _LOG.info(f"Command line | {' '.join(sys.argv)}")

    if not args.genome.is_file():
        _LOG.error(f"Genome file not found: {args.genome}")
        return 2
    if not args.summary.is_file():
        _LOG.error(f"Summary file not found: {args.summary}")
        return 2
    if pyrodigal_gv is None:
        _LOG.error(
            "pyrodigal-gv is not installed. Install it with: "
            "conda install -c bioconda pyrodigal-gv"
        )
        return 2

    versions = _get_tool_versions()
    _LOG.info("Tool versions | " + " | ".join(f"{key}={value}" for key, value in versions.items()))

    try:
        records = read_fasta(args.genome)
        intervals_by_contig, n_summary_rows = read_geve_intervals(args.summary)
        contig_lengths, masked_bp = validate_intervals(records, intervals_by_contig)
    except (OSError, UnicodeError, ValueError) as exc:
        _LOG.error(str(exc))
        return 2

    n_merged_intervals = sum(len(intervals) for intervals in intervals_by_contig.values())
    total_bp = sum(contig_lengths.values())
    _LOG.info(
        f"Genome: {len(records):,} contig(s), {total_bp:,} bp | "
        f"GEVE summary: {n_summary_rows:,} row(s), "
        f"{n_merged_intervals:,} merged exact interval(s) on "
        f"{len(intervals_by_contig):,} contig(s), {masked_bp:,} bp masked"
    )

    try:
        orfs, n_raw, n_overlap_removed = predict_orfs(
            records,
            intervals_by_contig,
            args.threads,
        )
    except RuntimeError as exc:
        _LOG.error(str(exc))
        return 2

    assign_orf_ids(orfs, args.prefix)

    # Final independent invariant check before any output is written.
    residual_overlaps: List[OrfRecord] = []
    for orf in orfs:
        intervals = intervals_by_contig.get(orf.contig, [])
        starts = [interval.start for interval in intervals]
        if _overlaps_intervals(starts, intervals, orf.start, orf.end):
            residual_overlaps.append(orf)
    if residual_overlaps:
        preview = ", ".join(
            f"{orf.contig}:{orf.start}-{orf.end}" for orf in residual_overlaps[:5]
        )
        _LOG.error(
            f"Internal safety check failed: {len(residual_overlaps)} retained ORF(s) "
            f"still overlap GEVE intervals: {preview}"
        )
        return 2

    pep_path = args.outdir / f"{args.prefix}.nonGEVE.pep"
    gff_path = args.outdir / f"{args.prefix}.nonGEVE.gff3"

    try:
        write_protein_fasta(orfs, pep_path)
        write_gff3(orfs, contig_lengths, gff_path)
    except OSError as exc:
        _LOG.error(f"Failed to write output files: {exc}")
        for path in (pep_path, gff_path):
            try:
                path.unlink(missing_ok=True)
            except OSError:
                pass
        return 2

    _LOG.info(
        f"ORF prediction: {n_raw:,} raw ORF(s), {n_overlap_removed:,} removed "
        f"by the explicit GEVE-overlap safety filter, {len(orfs):,} retained"
    )
    if n_overlap_removed:
        _LOG.warning(
            f"The overlap safety filter removed {n_overlap_removed:,} ORF(s) returned "
            f"despite sequence masking; none were written to the final outputs."
        )
    if not orfs:
        _LOG.warning("No ORFs were retained outside the GEVE regions.")

    _LOG.output(f"Protein FASTA -> {pep_path}")
    _LOG.output(f"GFF3 annotation -> {gff_path}")
    _LOG.output(f"Run log -> {log_path}")
    _LOG.info(f"run_prodigal completed in {_format_elapsed(time.time() - start_time)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
