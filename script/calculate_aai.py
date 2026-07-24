#!/usr/bin/env python3
"""
calculate_aai_v2.py - Calculate all-vs-all AAI among GEVE nucleotide or protein sequences.
Author: Dede Kurniawan (dedekurniawan@genomics.cn)
"""

from __future__ import annotations

import argparse
import logging
import math
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

try:
    import numpy as np
except ImportError:  # pragma: no cover
    np = None

try:
    import pandas as pd
except ImportError:  # pragma: no cover
    pd = None

try:
    import pyfastx
except ImportError:  # pragma: no cover
    pyfastx = None

try:
    import pyrodigal_gv
except ImportError:  # pragma: no cover
    pyrodigal_gv = None

HELP_TEXT = """\
calculate_aai_v2.py - Calculate all-vs-all AAI among GEVE nucleotide or protein sequences.

Usage: calculate_aai_v2.py <input.fa> --prefix <prefix> [OPTIONS]

Mandatory:
  input                 GEVE nucleotide FASTA or protein FASTA from findGEVE
  --prefix              Output prefix

Options:
  -o, --outdir          Output directory                              [default: ./AAI_Result_<YYYYMMDD>]
  -t, --threads         CPU threads for pyrodigal-gv and LAST          [default: 4]
  --evalue              LAST E-value cutoff                           [default: 1e-3]
  --min-query-cover     Minimum query coverage (%)                     [default: 50]
  --min-subject-cover   Minimum subject coverage (%)                   [default: 50]
  --min-aa-length       Minimum predicted/provided protein length      [default: 30]
  --max-heatmap         Maximum number of GEVEs for heatmap PDF        [default: 50]
  -h, --help            Show this help and exit
"""

USAGE_TEXT = "Usage: calculate_aai_v2.py <input.fa> --prefix <prefix> [OPTIONS]\n"

DEFAULTS = dict(
    threads=4,
    evalue=1e-3,
    min_query_cover=50.0,
    min_subject_cover=50.0,
    min_aa_length=30,
    max_heatmap=50,
)

_NATKEY_RE = re.compile(r"(\d+)")
_ORF_SUFFIX_RE = re.compile(r"(.+?)_orf\d+$", re.IGNORECASE)
_VALID_AA_RE = re.compile(r"^[A-Z*.-]+$")
_DNA_ALPHABET = set("ACGTUNRYSWKMBDHVacgtunryswkmbdhv")

_LOG = logging.getLogger("calculate_aai")

OUTPUT = 25
logging.addLevelName(OUTPUT, "OUTPUT")


def _output(self, message, *args, **kwargs):
    if self.isEnabledFor(OUTPUT):
        self._log(OUTPUT, message, args, **kwargs)


logging.Logger.output = _output


def _natural_key(s: str):
    return [int(t) if t.isdigit() else t.lower() for t in _NATKEY_RE.split(str(s))]


def setup_logging(log_path: Optional[Path] = None) -> None:
    _LOG.setLevel(logging.DEBUG)
    _LOG.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")
    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(fmt)
    sh.setLevel(logging.INFO)
    _LOG.addHandler(sh)
    if log_path is not None:
        fh = logging.FileHandler(log_path, mode="w")
        fh.setFormatter(fmt)
        fh.setLevel(logging.DEBUG)
        _LOG.addHandler(fh)


@dataclass
class ProteinRecord:
    protein_id: str
    geve_id: str
    sequence: str
    length: int


@dataclass
class GeveInfo:
    geve_id: str
    length: Optional[int] = None
    gc: Optional[float] = None
    n_proteins: int = 0


@dataclass
class BestHit:
    query: str
    subject: str
    query_geve: str
    subject_geve: str
    pident: float
    bitscore: float
    evalue: float
    aln_len: int
    qcover: float
    scover: float


@dataclass
class AaiResult:
    query_geve: str
    target_geve: str
    aai: Optional[float]
    num_rbh: int
    query_proteins: int
    target_proteins: int
    af_query: float
    af_target: float
    mean_af: float
    min_af: float
    query_len: Optional[int]
    target_len: Optional[int]
    query_gc: Optional[float]
    target_gc: Optional[float]
    similarity_class: str


# Stage 1: command-line parsing and checks
class _Parser(argparse.ArgumentParser):
    def format_help(self) -> str:
        return HELP_TEXT

    def format_usage(self) -> str:
        return USAGE_TEXT


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = _Parser(prog="calculate_aai_v2.py", add_help=False)
    p.add_argument("-h", "--help", action="help", help="Show this help and exit")
    p.add_argument("input", type=Path)
    p.add_argument("--prefix", type=str, required=True)
    p.add_argument("-o", "--outdir", type=Path,
                   default=Path(f"AAI_Result_{datetime.now().strftime('%Y%m%d')}"))
    p.add_argument("-t", "--threads", type=int, default=DEFAULTS["threads"])
    p.add_argument("--evalue", type=float, default=DEFAULTS["evalue"])
    p.add_argument("--min-query-cover", type=float, default=DEFAULTS["min_query_cover"])
    p.add_argument("--min-subject-cover", type=float, default=DEFAULTS["min_subject_cover"])
    p.add_argument("--min-aa-length", type=int, default=DEFAULTS["min_aa_length"])
    p.add_argument("--max-heatmap", type=int, default=DEFAULTS["max_heatmap"])
    return p.parse_args(argv)


def require_executable(name: str) -> None:
    if shutil.which(name) is None:
        _LOG.error(f"Required executable not found in PATH: {name}")
        sys.exit(1)


def require_python_package(obj, package_name: str) -> None:
    if obj is None:
        _LOG.error(f"Required Python package is not installed: {package_name}")
        sys.exit(1)


# Stage 2: input detection and FASTA parsing
def gc_of_seq(seq: str) -> float:
    arr = np.frombuffer(seq.upper().encode("ascii", "ignore"), dtype=np.uint8)
    if arr.size == 0:
        return float("nan")
    gc = int(((arr == ord("G")) | (arr == ord("C"))).sum())
    at = int(((arr == ord("A")) | (arr == ord("T"))).sum())
    valid = gc + at
    return float(100.0 * gc / valid) if valid > 0 else float("nan")


def is_probably_nucleotide(seq: str) -> bool:
    clean = re.sub(r"\s+", "", seq)
    clean = clean.replace("-", "").replace(".", "")
    if not clean:
        return False
    impossible_dna = set(clean) - _DNA_ALPHABET
    return len(impossible_dna) == 0


def detect_input_type(fasta_path: Path) -> str:
    try:
        fa = pyfastx.Fasta(str(fasta_path), build_index=True, uppercase=True)
    except Exception as exc:
        _LOG.error(f"Could not read input FASTA: {fasta_path} ({exc})")
        sys.exit(1)

    checked = 0
    nuc_like = 0
    for rec in fa:
        checked += 1
        if is_probably_nucleotide(str(rec.seq)):
            nuc_like += 1
        if checked >= 20:
            break
    if checked == 0:
        _LOG.error("Input FASTA contains no sequences.")
        sys.exit(1)
    return "nucleotide" if nuc_like == checked else "protein"


def geve_id_from_protein_header(header: str) -> str:
    base = header.split()[0]
    m = _ORF_SUFFIX_RE.match(base)
    if m:
        return m.group(1)
    return base


def sanitize_protein(seq: str) -> str:
    seq = re.sub(r"\s+", "", seq).upper().replace("*", "")
    seq = seq.replace(".", "").replace("-", "")
    return seq


def read_protein_input(
    fasta_path: Path,
    min_aa_length: int,
) -> Tuple[Dict[str, List[ProteinRecord]], Dict[str, GeveInfo], Dict[str, int]]:
    proteins_by_geve: Dict[str, List[ProteinRecord]] = defaultdict(list)
    geve_info: Dict[str, GeveInfo] = {}
    prot_lengths: Dict[str, int] = {}

    fa = pyfastx.Fasta(str(fasta_path), build_index=True, uppercase=True)
    n_read = 0
    n_kept = 0
    for rec in fa:
        n_read += 1
        pid = rec.name.split()[0]
        geve_id = geve_id_from_protein_header(rec.name)
        seq = sanitize_protein(str(rec.seq))
        if len(seq) < min_aa_length:
            continue
        if not _VALID_AA_RE.match(seq):
            _LOG.warning(f"Skipping protein {pid}: unsupported residue characters detected")
            continue
        pr = ProteinRecord(pid, geve_id, seq, len(seq))
        proteins_by_geve[geve_id].append(pr)
        prot_lengths[pid] = len(seq)
        geve_info.setdefault(geve_id, GeveInfo(geve_id=geve_id))
        n_kept += 1

    for gid, prots in proteins_by_geve.items():
        geve_info[gid].n_proteins = len(prots)

    _LOG.info(
        f"Protein FASTA input: {n_kept:,}/{n_read:,} protein(s) retained "
        f"across {len(proteins_by_geve):,} GEVE(s) "
        f"(min-aa-length={min_aa_length})"
    )
    _LOG.warning("Protein input has no nucleotide GEVE length or GC information; query_len/target_len/query_gc/target_gc will be NA.")
    return dict(proteins_by_geve), geve_info, prot_lengths


# Stage 3: pyrodigal-gv ORF prediction for nucleotide GEVE FASTA
def _predict_orfs_one_geve(args: Tuple[str, str, int]) -> Tuple[str, List[Tuple[str, str]], Optional[str]]:
    geve_id, seq, min_aa_length = args
    try:
        gf = pyrodigal_gv.ViralGeneFinder(meta=True)
        genes = gf.find_genes(seq.encode("ascii", "ignore"))
    except Exception as exc:
        return geve_id, [], f"pyrodigal-gv failed on {geve_id}: {exc}"

    out: List[Tuple[str, str]] = []
    for i, gene in enumerate(genes, start=1):
        prot = gene.translate().rstrip("*")
        prot = sanitize_protein(prot)
        if len(prot) < min_aa_length:
            continue
        out.append((f"{geve_id}_orf{i:05d}", prot))
    return geve_id, out, None


def predict_proteins_from_nucleotide(
    fasta_path: Path,
    min_aa_length: int,
    threads: int,
) -> Tuple[Dict[str, List[ProteinRecord]], Dict[str, GeveInfo], Dict[str, int]]:
    fa = pyfastx.Fasta(str(fasta_path), build_index=True, uppercase=True)
    work_items: List[Tuple[str, str, int]] = []
    geve_info: Dict[str, GeveInfo] = {}

    for rec in fa:
        gid = rec.name.split()[0]
        seq = str(rec.seq).upper()
        geve_info[gid] = GeveInfo(geve_id=gid, length=len(seq), gc=gc_of_seq(seq))
        work_items.append((gid, seq, min_aa_length))

    if not work_items:
        _LOG.error("Input FASTA contains no nucleotide GEVE sequences.")
        sys.exit(1)

    proteins_by_geve: Dict[str, List[ProteinRecord]] = defaultdict(list)
    prot_lengths: Dict[str, int] = {}
    n_warn = 0

    executor = None
    if threads > 1 and len(work_items) > 1:
        executor = ProcessPoolExecutor(max_workers=threads)
        results_iter = executor.map(_predict_orfs_one_geve, work_items, chunksize=1)
    else:
        results_iter = (_predict_orfs_one_geve(wi) for wi in work_items)

    try:
        for gid, records, err in results_iter:
            if err:
                n_warn += 1
                _LOG.warning(err)
                continue
            for pid, seq in records:
                pr = ProteinRecord(pid, gid, seq, len(seq))
                proteins_by_geve[gid].append(pr)
                prot_lengths[pid] = len(seq)
    finally:
        if executor is not None:
            executor.shutdown(wait=True)

    for gid in geve_info:
        geve_info[gid].n_proteins = len(proteins_by_geve.get(gid, []))

    total_prots = sum(len(v) for v in proteins_by_geve.values())
    _LOG.info(
        f"pyrodigal-gv ORF prediction: {total_prots:,} protein(s) retained "
        f"across {len(work_items):,} GEVE sequence(s) "
        f"(min-aa-length={min_aa_length}; warnings={n_warn})"
    )
    return dict(proteins_by_geve), geve_info, prot_lengths


# Stage 4: temporary LAST input construction
def write_combined_protein_fasta(
    proteins_by_geve: Dict[str, List[ProteinRecord]],
    out_path: Path,
) -> Tuple[Dict[str, str], Dict[str, int]]:
    protein_to_geve: Dict[str, str] = {}
    protein_lengths: Dict[str, int] = {}
    seen: set = set()
    with open(out_path, "w") as fh:
        for gid in sorted(proteins_by_geve, key=_natural_key):
            for p in proteins_by_geve[gid]:
                pid = p.protein_id
                if pid in seen:
                    # LAST and RBH parsing require unique protein IDs.
                    suffix = 2
                    new_pid = f"{pid}__dup{suffix}"
                    while new_pid in seen:
                        suffix += 1
                        new_pid = f"{pid}__dup{suffix}"
                    _LOG.warning(f"Duplicate protein ID {pid}; renamed to {new_pid}")
                    pid = new_pid
                seen.add(pid)
                protein_to_geve[pid] = gid
                protein_lengths[pid] = p.length
                fh.write(f">{pid}\n")
                seq = p.sequence
                for i in range(0, len(seq), 80):
                    fh.write(seq[i:i + 80] + "\n")
    return protein_to_geve, protein_lengths


def run_lastdb(db_prefix: Path, proteins_faa: Path) -> None:
    cmd = ["lastdb", "-p", str(db_prefix), str(proteins_faa)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd, output=proc.stdout, stderr=proc.stderr)


# Stage 5: LAST parsing and reciprocal best hits
def _is_better_hit(new: BestHit, old: Optional[BestHit]) -> bool:
    if old is None:
        return True
    if new.bitscore != old.bitscore:
        return new.bitscore > old.bitscore
    if new.pident != old.pident:
        return new.pident > old.pident
    if new.aln_len != old.aln_len:
        return new.aln_len > old.aln_len
    return new.evalue < old.evalue


def run_lastal_collect_best_hits(
    db_prefix: Path,
    query_faa: Path,
    protein_to_geve: Dict[str, str],
    protein_lengths: Dict[str, int],
    threads: int,
    evalue_cutoff: float,
    min_query_cover: float,
    min_subject_cover: float,
) -> Dict[Tuple[str, str], BestHit]:
    cmd = [
        "lastal", "-P", str(max(1, int(threads))), "-m", "500",
        str(db_prefix), str(query_faa), "-f", "BlastTab",
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert proc.stdout is not None

    best_hits: Dict[Tuple[str, str], BestHit] = {}
    n_lines = 0
    n_kept = 0

    for raw in proc.stdout:
        if not raw or raw.startswith("#"):
            continue
        n_lines += 1
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < 12:
            continue
        try:
            qid = fields[0]
            sid = fields[1]
            pident = float(fields[2])
            aln_len = int(float(fields[3]))
            evalue = float(fields[10])
            bitscore = float(fields[11])
        except (ValueError, IndexError):
            continue

        qgeve = protein_to_geve.get(qid)
        sgeve = protein_to_geve.get(sid)
        if qgeve is None or sgeve is None or qgeve == sgeve:
            continue
        if evalue > evalue_cutoff:
            continue
        qlen = protein_lengths.get(qid, 0)
        slen = protein_lengths.get(sid, 0)
        if qlen <= 0 or slen <= 0:
            continue
        qcover = 100.0 * float(aln_len) / float(qlen)
        scover = 100.0 * float(aln_len) / float(slen)
        if qcover < min_query_cover or scover < min_subject_cover:
            continue

        hit = BestHit(
            query=qid, subject=sid,
            query_geve=qgeve, subject_geve=sgeve,
            pident=pident, bitscore=bitscore, evalue=evalue,
            aln_len=aln_len, qcover=qcover, scover=scover,
        )
        key = (qid, sgeve)
        old = best_hits.get(key)
        if _is_better_hit(hit, old):
            best_hits[key] = hit
            n_kept += 1

    stderr = proc.stderr.read() if proc.stderr is not None else ""
    returncode = proc.wait()
    if returncode != 0:
        raise subprocess.CalledProcessError(returncode, cmd, output="", stderr=stderr)

    _LOG.info(
        f"LAST search parsed: {n_lines:,} alignment line(s); "
        f"{len(best_hits):,} best query-protein/target-GEVE hit(s) retained"
    )
    return best_hits


def classify_similarity(aai: Optional[float], min_af: float, num_rbh: int) -> str:
    if aai is None or num_rbh == 0:
        return "no_reciprocal_hits"
    if aai >= 80.0 and min_af >= 50.0:
        return "high_similarity"
    if aai >= 80.0 and min_af < 50.0:
        return "high_aai_low_af"
    if aai >= 50.0:
        return "moderate_similarity"
    return "low_similarity"


def compute_aai_results(
    geve_ids: List[str],
    geve_info: Dict[str, GeveInfo],
    best_hits: Dict[Tuple[str, str], BestHit],
) -> List[AaiResult]:
    pair_to_identities: Dict[Tuple[str, str], List[float]] = defaultdict(list)
    seen_pairs: set = set()

    for (qid, target_geve), hit in best_hits.items():
        qgeve = hit.query_geve
        sgeve = hit.subject_geve
        sid = hit.subject
        if qgeve == sgeve:
            continue
        reciprocal = best_hits.get((sid, qgeve))
        if reciprocal is None or reciprocal.subject != qid:
            continue
        protein_pair_key = tuple(sorted((qid, sid)))
        geve_pair_key = tuple(sorted((qgeve, sgeve), key=_natural_key))
        combined_key = (geve_pair_key, protein_pair_key)
        if combined_key in seen_pairs:
            continue
        seen_pairs.add(combined_key)
        mean_identity = float(np.mean([hit.pident, reciprocal.pident]))
        pair_to_identities[geve_pair_key].append(mean_identity)

    results: List[AaiResult] = []
    geve_ids = sorted(geve_ids, key=_natural_key)
    for i, g1 in enumerate(geve_ids):
        for g2 in geve_ids[i + 1:]:
            vals = pair_to_identities.get((g1, g2), [])
            n1 = int(geve_info[g1].n_proteins)
            n2 = int(geve_info[g2].n_proteins)
            num_rbh = len(vals)
            aai = float(np.mean(vals)) if vals else None
            af1 = 100.0 * num_rbh / n1 if n1 else 0.0
            af2 = 100.0 * num_rbh / n2 if n2 else 0.0
            mean_af = (af1 + af2) / 2.0
            min_af = min(af1, af2)
            results.append(AaiResult(
                query_geve=g1, target_geve=g2,
                aai=aai, num_rbh=num_rbh,
                query_proteins=n1, target_proteins=n2,
                af_query=af1, af_target=af2,
                mean_af=mean_af, min_af=min_af,
                query_len=geve_info[g1].length,
                target_len=geve_info[g2].length,
                query_gc=geve_info[g1].gc,
                target_gc=geve_info[g2].gc,
                similarity_class=classify_similarity(aai, min_af, num_rbh),
            ))
    return results


# Stage 6: output writing
def _fmt_float(x: Optional[float], ndigits: int = 3) -> str:
    if x is None:
        return "NA"
    try:
        if math.isnan(float(x)):
            return "NA"
    except Exception:
        return "NA"
    return f"{float(x):.{ndigits}f}"


def _fmt_int(x: Optional[int]) -> str:
    return "NA" if x is None else str(int(x))


def write_pairs_tsv(results: List[AaiResult], out_path: Path) -> None:
    columns = [
        "query_geve", "target_geve", "aai", "num_rbh",
        "query_proteins", "target_proteins",
        "af_query", "af_target", "mean_af", "min_af",
        "query_len", "target_len", "query_gc", "target_gc",
        "similarity_class",
    ]
    with open(out_path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in results:
            row = [
                r.query_geve,
                r.target_geve,
                _fmt_float(r.aai, 3),
                str(r.num_rbh),
                str(r.query_proteins),
                str(r.target_proteins),
                _fmt_float(r.af_query, 3),
                _fmt_float(r.af_target, 3),
                _fmt_float(r.mean_af, 3),
                _fmt_float(r.min_af, 3),
                _fmt_int(r.query_len),
                _fmt_int(r.target_len),
                _fmt_float(r.query_gc, 3),
                _fmt_float(r.target_gc, 3),
                r.similarity_class,
            ]
            fh.write("\t".join(row) + "\n")


def build_aai_matrix(geve_ids: List[str], results: List[AaiResult]) -> pd.DataFrame:
    geve_ids = sorted(geve_ids, key=_natural_key)
    mat = pd.DataFrame(np.full((len(geve_ids), len(geve_ids)), np.nan),
                       index=geve_ids, columns=geve_ids)
    for gid in geve_ids:
        mat.loc[gid, gid] = 100.0
    for r in results:
        value = float(r.aai) if r.aai is not None else 0.0
        mat.loc[r.query_geve, r.target_geve] = value
        mat.loc[r.target_geve, r.query_geve] = value
    return mat.fillna(0.0)


def write_heatmap_pdf(
    aai_matrix: pd.DataFrame,
    out_path: Path,
    max_heatmap: int,
) -> Optional[Path]:
    n = aai_matrix.shape[0]
    if n == 0:
        _LOG.warning("Heatmap skipped: empty AAI matrix")
        return None
    if n > max_heatmap:
        _LOG.info(
            f"Heatmap skipped: {n:,} GEVEs > --max-heatmap {max_heatmap:,}. "
            "Pairwise TSV was still written."
        )
        return None

    try:
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap
        from scipy.cluster.hierarchy import leaves_list, linkage
        from scipy.spatial.distance import squareform
    except Exception as exc:
        _LOG.warning(f"Heatmap skipped: scipy/matplotlib import failed ({exc})")
        return None

    values = aai_matrix.values.astype(float)
    if n > 1:
        dist = 100.0 - values
        np.fill_diagonal(dist, 0.0)
        dist = np.maximum(dist, 0.0)
        condensed = squareform(dist, checks=False)
        Z = linkage(condensed, method="average")
        order = leaves_list(Z)
    else:
        order = np.array([0])

    ordered = aai_matrix.iloc[order, order]
    fig_size = max(6.0, min(18.0, 0.28 * n + 3.0))
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    cmap = LinearSegmentedColormap.from_list(
        "white_to_geve_red",
        ["#ffffff", "#f76e68"],
        N=256,
    )
    im = ax.imshow(ordered.values, vmin=0, vmax=100, aspect="auto", cmap=cmap)
    ax.set_xticks(np.arange(n))
    ax.set_yticks(np.arange(n))
    ax.set_xticklabels(ordered.columns, rotation=90, fontsize=max(4, min(8, int(250 / max(n, 1)))))
    ax.set_yticklabels(ordered.index, fontsize=max(4, min(8, int(250 / max(n, 1)))))
    ax.set_title("GEVE all-vs-all AAI (%)")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("AAI (%)")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


# Main
def main(argv: Optional[List[str]] = None) -> int:
    t0 = time.time()
    args = parse_args(argv)

    if args.input.is_dir():
        sys.stderr.write("ERROR: folder input is not supported. Please provide one nucleotide or protein FASTA file.\n")
        return 1
    if not args.input.is_file():
        sys.stderr.write(f"ERROR: input FASTA not found: {args.input}\n")
        return 1

    args.outdir.mkdir(parents=True, exist_ok=True)
    log_path = args.outdir / "run.log"
    setup_logging(log_path)

    _LOG.info("calculate_aai started")
    _LOG.info(f"Input FASTA: {args.input}")
    _LOG.info(f"Output directory: {args.outdir}")

    require_python_package(np, "numpy")
    require_python_package(pd, "pandas")
    require_python_package(pyfastx, "pyfastx")
    require_executable("lastdb")
    require_executable("lastal")

    input_type = detect_input_type(args.input)
    _LOG.info(f"Detected input type: {input_type}")

    if input_type == "nucleotide":
        require_python_package(pyrodigal_gv, "pyrodigal-gv")
        proteins_by_geve, geve_info, _ = predict_proteins_from_nucleotide(
            args.input, args.min_aa_length, max(1, int(args.threads))
        )
    else:
        proteins_by_geve, geve_info, _ = read_protein_input(
            args.input, args.min_aa_length
        )

    # Drop GEVE records with no retained proteins.
    no_prot = [gid for gid, info in geve_info.items() if info.n_proteins == 0]
    if no_prot:
        _LOG.warning(f"{len(no_prot):,} GEVE(s) have no retained proteins and will be excluded from AAI.")
        for gid in no_prot:
            geve_info.pop(gid, None)
            proteins_by_geve.pop(gid, None)

    geve_ids = sorted(proteins_by_geve.keys(), key=_natural_key)
    if len(geve_ids) < 2:
        _LOG.error("At least two GEVE records with retained proteins are required for AAI.")
        return 1

    total_prots = sum(len(v) for v in proteins_by_geve.values())
    _LOG.info(f"AAI dataset: {len(geve_ids):,} GEVE(s), {total_prots:,} protein(s)")

    pairs_tsv = args.outdir / f"{args.prefix}.aai_pairs.tsv"
    heatmap_pdf = args.outdir / f"{args.prefix}.aai_heatmap.pdf"

    with tempfile.TemporaryDirectory(prefix="calculate_aai_") as tmpdir:
        tmp = Path(tmpdir)
        combined_faa = tmp / "all_geve_proteins.faa"
        db_prefix = tmp / "all_geve_proteins.lastdb"

        protein_to_geve, protein_lengths = write_combined_protein_fasta(proteins_by_geve, combined_faa)
        _LOG.info("Formatting LAST protein database")
        try:
            run_lastdb(db_prefix, combined_faa)
        except subprocess.CalledProcessError as exc:
            _LOG.error(f"lastdb failed: {(exc.stderr or exc.output or '').strip()[:500]}")
            return 1

        _LOG.info("Running LAST all-vs-all protein search")
        try:
            best_hits = run_lastal_collect_best_hits(
                db_prefix=db_prefix,
                query_faa=combined_faa,
                protein_to_geve=protein_to_geve,
                protein_lengths=protein_lengths,
                threads=max(1, int(args.threads)),
                evalue_cutoff=float(args.evalue),
                min_query_cover=float(args.min_query_cover),
                min_subject_cover=float(args.min_subject_cover),
            )
        except subprocess.CalledProcessError as exc:
            _LOG.error(f"lastal failed: {(exc.stderr or exc.output or '').strip()[:500]}")
            return 1

    _LOG.info("Computing reciprocal-best-hit AAI values")
    results = compute_aai_results(geve_ids, geve_info, best_hits)
    write_pairs_tsv(results, pairs_tsv)
    _LOG.output(f"AAI pair table -> {pairs_tsv}")

    aai_matrix = build_aai_matrix(geve_ids, results)
    hp = write_heatmap_pdf(aai_matrix, heatmap_pdf, int(args.max_heatmap))
    if hp is not None:
        _LOG.output(f"AAI heatmap    -> {hp}")

    _LOG.output(f"Run log        -> {log_path}")
    elapsed = time.time() - t0
    if elapsed >= 3600:
        timing = f"{elapsed / 3600:.2f} h"
    elif elapsed >= 60:
        timing = f"{elapsed / 60:.2f} min"
    else:
        timing = f"{elapsed:.1f} s"
    _LOG.info(f"calculate_aai completed in {timing}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
