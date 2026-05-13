#!/usr/bin/env python3
"""
findGEVE.py - Identify Giant Endogenous Viral Elements (GEVEs) in eukaryotic genome assemblies.
Author: Dede Kurniawan (dedekurniawan@genomics.cn)
"""

from __future__ import annotations

import argparse
import logging
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version as _pkg_version
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
import numpy as np
import pandas as pd
import pyfastx
import pyhmmer
import pyrodigal_gv

HELP_TEXT = """\
findGEVE.py - Identify Giant Endogenous Viral Elements (GEVEs) in eukaryotic genome assemblies.

Usage: findGEVE.py -db <directory> --prefix <prefix> <genome.fa> [OPTIONS]

Mandatory:
  -db, --db            HMM database directory (must contain NCLDV_markers.hmm
                       and gvog.complete.hmm; Pfam-A.hmm is optional)
  --prefix             Output prefix for GEVE IDs and file names
  genome               Input genome assembly FASTA (gzip is acceptable)

Optionals:
  -o, --outdir         Output directory                              [default: ./Result_<YYYYMMDD>]
  -t, --threads        CPU threads for ORF prediction and HMM search [default: 4]
  -e, --evalue         E-value cutoff for HMM searches               [default: 1e-5]
  --blastn-jobs        Parallel TIR-detection workers                [default: --threads]
  -m, --min-hallmark-type
                       Minimum number of distinct hallmark types
                       required in the final retained GEVE
                       (seeding always uses >= 1)                    [default: 2]
  --min-contig         Minimum contig length to scan for GEVE 
                       detection                                     [default: 50_000]
  --cluster-merge-gap  Maximum gap (bp) between same-contig clusters
                       eligible for merging                          [default: 100_000]
  -h, --help           Show this help and exit
"""

USAGE_TEXT = "Usage: findGEVE.py -db <DB directory> --prefix <prefix> <genome.fa> [OPTIONS]\n"

# tunable parameters.
DEFAULTS = dict(
    min_contig          = 50_000,
    min_geve_length     = 50_000,
    min_hallmarks_seed  = 1,
    min_hallmarks       = 2,
    seed_window         = 300_000,
    cluster_merge_gap   = 100_000,
    max_cluster_span    = 2_000_000,
    host_territory_fraction = 0.7,
    rolling_window      = 15,
    tir_flank_start     = 100_000,
    tir_flank_step      = 100_000,
    tir_flank_max       = 200_000,
    tir_min_insert      = 30_000,  
    tir_max_insert      = 1_500_000,
    tir_min_len         = 10,
    tir_max_len         = 10_000,
    tir_min_id          = 65.0, 
    tir_bracket_fraction= 1.0,
    tir_min_dinuc_entropy = 2.0,
    tir_max_kmer_fraction = 0.70,
    tir_max_tandem_fraction = 0.70,
    tir_tandem_max_period = 12,
    blastn_word_size    = 7,        
    blastn_reward       = 1,         
    blastn_penalty      = -1,        
    blastn_gapopen      = 2,         
    blastn_gapextend    = 1,         
    blastn_evalue       = 10.0,      
    blastn_max_targets  = 10_000,   
    tsd_min             = 4,
    tsd_max             = 12,
    tsd_max_slide       = 2,
    tsd_search_window   = 60,     
    extend_tirless      = True,
    extend_threshold    = -1.0,
    extend_start_threshold = 0.0,
    extend_max_bp       = 200_000,
    extend_max_drops    = 5,
    n_max_fraction      = 0.05,
    evalue              = 1e-5,
    hallmark_score_cutoffs = {
        "A32":   80.0,
        "D5":    80.0,
        "SFII": 100.0,
        "mcp":   80.0,
        "mRNAc": 80.0,
        "PolB": 200.0,
        "RNAPL": 200.0,
        "RNAPS": 200.0,
        "RNR":   80.0,
        "VLTF3": 80.0,
    },
)

_NATKEY_RE = re.compile(r"(\d+)")

def _natural_key(s: str):
    return [int(t) if t.isdigit() else t.lower() for t in _NATKEY_RE.split(str(s))]

_LOG = logging.getLogger("findGEVE")

OUTPUT = 25
logging.addLevelName(OUTPUT, "OUTPUT")

def _output(self, message, *args, **kwargs):
    if self.isEnabledFor(OUTPUT):
        self._log(OUTPUT, message, args, **kwargs)

logging.Logger.output = _output

def setup_logging(log_path: Optional[Path] = None) -> None:
    _LOG.setLevel(logging.DEBUG)
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
class Orf:
    orf_id: str
    contig: str
    start: int          
    end: int             
    strand: int          
    protein: str
    hallmark: Optional[str] = None
    hallmark_bitscore: float = 0.0
    hallmark_evalue: float = float("inf")
    gvog: Optional[str] = None
    gvog_bitscore: float = 0.0
    gvog_evalue: float = float("inf")
    best_pfam_acc: Optional[str] = None
    best_pfam_name: Optional[str] = None
    best_pfam_bitscore: float = 0.0
    best_pfam_evalue: float = float("inf")
    virbit: float = 0.0
    pfambit: float = 0.0
    net_score: float = 0.0

@dataclass
class TirPair:
    left_start: int   
    left_end: int
    right_start: int
    right_end: int
    tir_length: int
    insert_size: int
    tir_identity: float
    score: int
    matches: int
    total: int
    gaps: int
    tir_evalue: float = float("nan")

@dataclass
class Tsd:
    sequence_left: str
    sequence_right: str
    length: int
    mismatches: int
    identity: float
    left_shift: int
    right_shift: int

# Stage 1: ORF prediction
def _predict_orfs_on_contig(
    args: Tuple[str, str]
) -> Tuple[str, List[Tuple], Optional[str]]:
    contig_name, seq = args
    try:
        gf = pyrodigal_gv.ViralGeneFinder(meta=True)
        genes = gf.find_genes(seq.encode("ascii"))
    except Exception as exc:
        return contig_name, [], f"pyrodigal-gv failed on {contig_name}: {exc}"
    out = []
    for i, gene in enumerate(genes, start=1):
        prot = gene.translate().rstrip("*")
        if not prot:
            continue
        out.append((
            f"orf{i:05d}",
            int(gene.begin), int(gene.end), int(gene.strand),
            prot,
        ))
    return contig_name, out, None

def predict_orfs(
    genome_path: Path,
    min_contig: int,
    threads: int,
) -> Tuple[Dict[str, Orf], Dict[str, int]]:
    """Run pyrodigal-gv meta on all contigs >= min_contig in parallel."""
    fa = pyfastx.Fasta(str(genome_path), build_index=True, uppercase=True)
    work_items: List[Tuple[str, str]] = []
    contig_lengths: Dict[str, int] = {}
    n_kept = n_skipped = 0
    for rec in fa:
        contig_lengths[rec.name] = len(rec.seq)
        if len(rec.seq) < min_contig:
            n_skipped += 1
            continue
        n_kept += 1
        work_items.append((rec.name, str(rec.seq)))

    if not work_items:
        _LOG.error("No contigs passed the minimum-length filter.")
        sys.exit(1)

    executor = None
    if threads > 1 and len(work_items) > 1:
        executor = ProcessPoolExecutor(max_workers=threads)
        results_iter = executor.map(_predict_orfs_on_contig, work_items, chunksize=1)
    else:
        results_iter = (_predict_orfs_on_contig(wi) for wi in work_items)

    orfs_by_id: Dict[str, Orf] = {}
    total_orfs = 0
    try:
        for contig_name, records, err in results_iter:
            if err:
                _LOG.warning(f"{err}; skipping contig")
                continue
            for orf_suffix, start, end, strand, prot in records:
                orf_id = f"{contig_name}__{orf_suffix}"
                orfs_by_id[orf_id] = Orf(
                    orf_id=orf_id, contig=contig_name,
                    start=start, end=end, strand=strand,
                    protein=prot,
                )
                total_orfs += 1
    finally:
        if executor is not None:
            executor.shutdown(wait=True)

    _LOG.info(
        f"ORF prediction: {total_orfs:,} ORFs on {n_kept:,} contig(s) "
        f"(>= {min_contig:,} bp; {n_skipped:,} skipped)"
    )
    if total_orfs == 0:
        _LOG.error("No ORFs predicted. Check input FASTA.")
        sys.exit(1)
    return orfs_by_id, contig_lengths

# Stage 2: HMM scans
_HMMER_MAX_TARGET_LEN = 100_000

def _digital_seqs(
    orfs: Iterable[Orf],
    alphabet: pyhmmer.easel.Alphabet,
) -> List[pyhmmer.easel.DigitalSequence]:
    out = []
    for o in orfs:
        if len(o.protein) > _HMMER_MAX_TARGET_LEN:
            _LOG.warning(
                f"Skipping ORF {o.orf_id} on {o.contig}: protein length "
                f"{len(o.protein):,} aa exceeds HMMER pipeline limit "
                f"({_HMMER_MAX_TARGET_LEN:,} aa)"
            )
            continue
        try:
            t = pyhmmer.easel.TextSequence(
                name=o.orf_id.encode(), sequence=o.protein
            )
            out.append(t.digitize(alphabet))
        except Exception as exc:
            _LOG.warning(f"Skipping ORF {o.orf_id} during digitization: {exc}")
    return out

def _hmm_query_name(top_hits) -> str:
    try:
        n = top_hits.query.name
    except AttributeError:
        n = top_hits.query_name
    return n.decode() if isinstance(n, bytes) else n

def scan_hallmarks(
    orfs_by_id: Dict[str, Orf],
    hmm_path: Path,
    evalue: float,
    threads: int,
    score_cutoffs: Optional[Dict[str, float]] = None,
) -> Dict[str, List[str]]:
    """Scan all proteins against NCLDV_markers.hmm.
    """
    alphabet = pyhmmer.easel.Alphabet.amino()
    seqs = _digital_seqs(orfs_by_id.values(), alphabet)

    with pyhmmer.plan7.HMMFile(str(hmm_path)) as hf:
        hmms = list(hf)
    if not hmms:
        _LOG.error(f"No HMM profiles in {hmm_path}")
        sys.exit(1)

    score_cutoffs = score_cutoffs or {}
    contig2hits: Dict[str, List[str]] = defaultdict(list)
    n_hits = 0
    for top_hits in pyhmmer.hmmsearch(hmms, seqs, cpus=threads, E=evalue):
        hmm_name = _hmm_query_name(top_hits)
        cutoff = score_cutoffs.get(hmm_name, 0.0)
        for hit in top_hits:
            if not hit.included:
                continue
            score = float(hit.score)
            if score < cutoff:
                continue
            target = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
            o = orfs_by_id.get(target)
            if o is None:
                continue
            if score > o.hallmark_bitscore:
                o.hallmark = hmm_name
                o.hallmark_bitscore = score
                o.hallmark_evalue = float(hit.evalue)
                o.virbit = score          # provisional; overwritten by GVOG scan
            contig2hits[o.contig].append(hmm_name)
            n_hits += 1

    _LOG.info(
        f"Hallmark scan: {n_hits:,} hit(s) on {len(contig2hits):,} contig(s)"
    )
    return dict(contig2hits)

def scan_gvog(
    orfs_to_scan: List[Orf],
    hmm_path: Path,
    evalue: float,
    threads: int,
) -> None:
    """Scan a subset of proteins against gvog.complete.hmm Sets Orf.virbit."""
    if not orfs_to_scan:
        return
    alphabet = pyhmmer.easel.Alphabet.amino()
    seqs = _digital_seqs(orfs_to_scan, alphabet)
    by_id = {o.orf_id: o for o in orfs_to_scan}

    n_hits = 0
    with pyhmmer.plan7.HMMFile(str(hmm_path)) as hf:
        for top_hits in pyhmmer.hmmsearch(hf, seqs, cpus=threads, E=evalue):
            hmm_name = _hmm_query_name(top_hits)
            for hit in top_hits:
                if not hit.included:
                    continue
                target = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
                o = by_id.get(target)
                if o is None:
                    continue
                score = float(hit.score)
                if score > o.gvog_bitscore:
                    o.gvog = hmm_name
                    o.gvog_bitscore = score
                    o.gvog_evalue = float(hit.evalue)
                    o.virbit = score
                n_hits += 1

    _LOG.info(f"GVOG scan: {n_hits:,} hit(s)")

def scan_pfam(
    orfs_to_scan: List[Orf],
    hmm_path: Path,
    evalue: float,
    threads: int,
) -> None:
    """Scan a subset of proteins against Pfam-A.hmm. Sets Orf.pfambit."""
    if not orfs_to_scan:
        return
    alphabet = pyhmmer.easel.Alphabet.amino()
    seqs = _digital_seqs(orfs_to_scan, alphabet)
    by_id = {o.orf_id: o for o in orfs_to_scan}

    n_hits = 0
    with pyhmmer.plan7.HMMFile(str(hmm_path)) as hf:
        for top_hits in pyhmmer.hmmsearch(hf, seqs, cpus=threads, E=evalue):
            hmm_name = _hmm_query_name(top_hits)
            try:
                acc_raw = top_hits.query.accession
                acc = (acc_raw.decode() if isinstance(acc_raw, bytes) else acc_raw) or hmm_name
                acc = acc.split(".")[0]
            except AttributeError:
                acc = hmm_name
            for hit in top_hits:
                if not hit.included:
                    continue
                target = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
                o = by_id.get(target)
                if o is None:
                    continue
                score = float(hit.score)
                if score > o.best_pfam_bitscore:
                    o.best_pfam_bitscore = score
                    o.best_pfam_evalue = float(hit.evalue)
                    o.best_pfam_acc = acc
                    o.best_pfam_name = hmm_name
                    o.pfambit = score
                n_hits += 1

    _LOG.info(f"Pfam-A scan: {n_hits:,} hit(s)")

# Stage 3: seeding
def _gap_is_host_territory(
    prev_end: int,
    curr_start: int,
    contig_orfs: List[Orf],
    rolling: Dict[str, float],
    host_fraction: float = 0.7,
) -> bool:
    """Gap is host territory when no hallmark/GVOG hits and a majority of gap ORFs
    have negative rolling viral score (fraction > host_fraction)."""
    gap_orfs = [o for o in contig_orfs
                if o.start > prev_end and o.end < curr_start]
    if not gap_orfs:
        return False
    if any(o.hallmark is not None or o.gvog is not None for o in gap_orfs):
        return False
    n_neg = sum(1 for o in gap_orfs if rolling.get(o.orf_id, 0.0) < 0)
    return (n_neg / len(gap_orfs)) > host_fraction

def find_seed_clusters(
    orfs_by_contig: Dict[str, List[Orf]],
    window_size: int,
    min_hallmarks: int,
    cluster_merge_gap: int,
    max_cluster_span: int,
    rolling_by_orf_per_contig: Optional[Dict[str, Dict[str, float]]] = None,
    host_fraction: float = 0.7,
) -> List[dict]:
    """Identify candidate GEVE clusters by sliding window over hallmark ORFs."""
    half = window_size // 2
    raw: List[dict] = []
    rolling_by_orf_per_contig = rolling_by_orf_per_contig or {}

    for contig, orfs in orfs_by_contig.items():
        hallmark_orfs = [o for o in orfs if o.hallmark is not None]
        if len(hallmark_orfs) < min_hallmarks:
            continue
        hallmark_orfs.sort(key=lambda x: x.start)
        starts = np.array([o.start for o in hallmark_orfs])
        ends   = np.array([o.end   for o in hallmark_orfs])
        for anchor in hallmark_orfs:
            mid = (anchor.start + anchor.end) // 2
            wstart, wend = mid - half, mid + half
            mask = (starts <= wend) & (ends >= wstart)
            members = [hallmark_orfs[i] for i in np.nonzero(mask)[0]]
            fams = {m.hallmark for m in members}
            if len(fams) < min_hallmarks:
                continue
            raw.append(dict(
                contig=contig,
                cluster_start=min(m.start for m in members),
                cluster_end=max(m.end for m in members),
                orf_ids=[m.orf_id for m in members],
                hallmarks=sorted(fams),
            ))

    raw.sort(key=lambda c: (c["contig"], c["cluster_start"]))
    overlap_merged: List[dict] = []
    for c in raw:
        if (overlap_merged
                and overlap_merged[-1]["contig"] == c["contig"]
                and c["cluster_start"] <= overlap_merged[-1]["cluster_end"]):
            prev = overlap_merged[-1]
            prev["cluster_end"] = max(prev["cluster_end"], c["cluster_end"])
            prev["orf_ids"]   = sorted(set(prev["orf_ids"])   | set(c["orf_ids"]))
            prev["hallmarks"] = sorted(set(prev["hallmarks"]) | set(c["hallmarks"]))
        else:
            overlap_merged.append(dict(c))

    merged: List[dict] = []
    n_refused = 0
    for c in overlap_merged:
        if merged and merged[-1]["contig"] == c["contig"]:
            prev = merged[-1]
            gap = c["cluster_start"] - prev["cluster_end"] - 1
            new_span = c["cluster_end"] - prev["cluster_start"] + 1
            if 0 <= gap <= cluster_merge_gap and new_span <= max_cluster_span:
                contig_orfs = orfs_by_contig.get(c["contig"], [])
                rolling = rolling_by_orf_per_contig.get(c["contig"], {})
                if not _gap_is_host_territory(
                        prev["cluster_end"], c["cluster_start"], contig_orfs, rolling,
                        host_fraction):
                    prev["cluster_end"] = c["cluster_end"]
                    prev["orf_ids"]   = sorted(set(prev["orf_ids"])   | set(c["orf_ids"]))
                    prev["hallmarks"] = sorted(set(prev["hallmarks"]) | set(c["hallmarks"]))
                    continue
                n_refused += 1
        merged.append(dict(c))

    _LOG.info(
        f"Seeding: {len(merged)} candidate region(s) "
        f"(>= {min_hallmarks} distinct hallmarks within {window_size:,} bp; "
        f"merged within {cluster_merge_gap:,} bp, max span {max_cluster_span:,} bp; "
        f"{n_refused} merge(s) refused as host territory)"
    )
    return merged

# Stage 4: rolling viral score
def compute_rolling_scores(
    orfs_by_contig: Dict[str, List[Orf]],
    rolling_window: int,
    candidate_contigs: set,
) -> Dict[str, Dict[str, float]]:
    """Centered rolling mean of net_score per ORF, vectorized via cumulative sums."""
    half = rolling_window // 2
    result: Dict[str, Dict[str, float]] = {}
    for contig in candidate_contigs:
        orfs = sorted(orfs_by_contig.get(contig, []), key=lambda x: x.start)
        if not orfs:
            continue
        scores = np.fromiter((o.net_score for o in orfs), dtype=np.float64, count=len(orfs))
        n = scores.size
        if n == 0:
            continue
        idx = np.arange(n)
        lo = np.maximum(0, idx - half)
        hi = np.minimum(n, idx + half + 1)
        cumsum = np.concatenate(([0.0], np.cumsum(scores)))
        window_sums = cumsum[hi] - cumsum[lo]
        window_sizes = (hi - lo).astype(np.float64)
        rolling = np.where(window_sizes >= 3, window_sums / np.maximum(window_sizes, 1), 0.0)
        result[contig] = {orfs[i].orf_id: float(rolling[i]) for i in range(n)}
    return result

# Stage 5: TIR detection (blastn self-vs-self, -strand minus)
_BLASTN_OUTFMT = "6 qstart qend sstart send length nident pident gaps evalue bitscore"
_BLASTN_COLUMNS = [
    "qstart", "qend", "sstart", "send", "length", "nident", "pident",
    "gaps", "evalue", "bitscore",
]

def parse_blastn_tabular(tab_path: Path) -> List[TirPair]:
    """Parse blastn outfmt-6 output into TirPair objects."""
    if not tab_path.exists():
        return []
    text = tab_path.read_text()
    if not text.strip():
        return []

    pairs: List[TirPair] = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < len(_BLASTN_COLUMNS):
            continue
        try:
            qstart   = int(fields[0])
            qend     = int(fields[1])
            sstart   = int(fields[2])
            send     = int(fields[3])
            aln_len  = int(fields[4])
            nident   = int(fields[5])
            pident   = float(fields[6])
            gaps     = int(fields[7])
            evalue   = float(fields[8])
            bitscore = float(fields[9])
        except (ValueError, IndexError):
            continue

        left_local_start  = min(qstart, qend)
        left_local_end    = max(qstart, qend)
        right_local_start = min(sstart, send)
        right_local_end   = max(sstart, send)

        if left_local_end >= right_local_start:
            continue
        if left_local_start > right_local_start:
            continue

        tir_length  = left_local_end - left_local_start + 1
        insert_size = right_local_start - left_local_end - 1

        pairs.append(TirPair(
            left_start=left_local_start,
            left_end=left_local_end,
            right_start=right_local_start,
            right_end=right_local_end,
            tir_length=tir_length,
            insert_size=insert_size,
            tir_identity=pident,
            score=int(round(bitscore)),
            matches=nident,
            total=aln_len,
            gaps=gaps,
            tir_evalue=evalue,
        ))
    return pairs

def run_blastn_self(
    region_fa_path: Path,
    tab_out_path: Path,
    cfg: dict,
    threads: int = 1,
) -> None:
    """Run blastn -query R -subject R -strand minus to find inverted repeats."""
    cmd = [
        "blastn",
        "-query",         str(region_fa_path),
        "-subject",       str(region_fa_path),
        "-strand",        "minus",
        "-task",          "blastn",
        "-word_size",     str(cfg.get("blastn_word_size", 7)),
        "-reward",        str(cfg.get("blastn_reward", 1)),
        "-penalty",       str(cfg.get("blastn_penalty", -1)),
        "-gapopen",       str(cfg.get("blastn_gapopen", 2)),
        "-gapextend",     str(cfg.get("blastn_gapextend", 1)),
        "-evalue",        str(cfg.get("blastn_evalue", 10.0)),
        "-dust",          "no",
        "-soft_masking",  "false",
        "-max_target_seqs", str(cfg.get("blastn_max_targets", 10000)),
        "-num_threads",   str(max(1, int(threads))),
        "-outfmt",        _BLASTN_OUTFMT,
        "-out",           str(tab_out_path),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(
            proc.returncode, cmd, output=proc.stdout, stderr=proc.stderr
        )

def _dinucleotide_entropy(seq: str) -> float:
    """Shannon entropy (bits) of overlapping ACGT dinucleotides."""
    if len(seq) < 2:
        return 0.0
    seq = seq.upper()
    counts: Counter = Counter()
    for i in range(len(seq) - 1):
        di = seq[i:i + 2]
        if "N" in di:
            continue
        counts[di] += 1
    total = sum(counts.values())
    if total == 0:
        return 0.0
    entropy = 0.0
    for c in counts.values():
        p = c / total
        entropy -= p * np.log2(p)
    return entropy

def _max_kmer_fraction(seq: str, k: int) -> float:
    """Maximum fraction of any single k-mer across all phase offsets."""
    arr = seq.upper().encode("ascii")
    n = len(arr)
    if n < k:
        return 0.0
    best = 0.0
    for phase in range(k):
        tiles = [arr[i:i + k] for i in range(phase, n - k + 1, k)]
        valid = [t for t in tiles if b"N" not in t]
        if not valid:
            continue
        counts = Counter(valid)
        max_count = max(counts.values())
        frac = max_count / len(valid)
        if frac > best:
            best = frac
    return best

def _max_tandem_period_fraction(seq: str, max_period: int) -> float:
    """Maximum fraction of positions i where seq[i] == seq[i+p] across periods p=1..max_period."""
    arr = np.frombuffer(seq.upper().encode("ascii"), dtype=np.uint8)
    n = arr.size
    if n < 2 or max_period < 1:
        return 0.0
    n_byte = ord("N")
    upper_p = min(max_period, n - 1)
    best = 0.0
    for p in range(1, upper_p + 1):
        a = arr[:-p]
        b = arr[p:]
        valid = (a != n_byte) & (b != n_byte)
        total = int(valid.sum())
        if total == 0:
            continue
        matches = int((valid & (a == b)).sum())
        frac = matches / total
        if frac > best:
            best = frac
    return best

def _tir_is_low_complexity(
    seq: str,
    min_entropy: float,
    max_kmer_fraction: float,
    max_tandem_fraction: float,
    max_tandem_period: int,
) -> Tuple[bool, str]:
    """Reject TIRs that are simple repeats or low-complexity sequences.
    """
    if not seq:
        return True, "empty TIR sequence"
    entropy = _dinucleotide_entropy(seq)
    if entropy < min_entropy:
        return True, f"dinuc_entropy={entropy:.2f}<{min_entropy}"
    for k in (1, 2, 3, 4):
        frac = _max_kmer_fraction(seq, k)
        if frac > max_kmer_fraction:
            return True, f"{k}-mer_fraction={frac:.2f}>{max_kmer_fraction}"
    tandem_frac = _max_tandem_period_fraction(seq, max_tandem_period)
    if tandem_frac > max_tandem_fraction:
        return True, f"tandem_fraction={tandem_frac:.2f}>{max_tandem_fraction}"
    return False, ""

def _count_bracketed(
    tir: TirPair,
    hallmark_intervals: List[Tuple[int, int]],
) -> int:
    """Count how many hallmark ORFs fall fully inside a TIR pair (genome-absolute)."""
    return sum(
        1 for hs, he in hallmark_intervals
        if tir.left_start <= hs and he <= tir.right_end
    )

def select_best_tir(
    pairs: List[TirPair],
    region_offset: int,
    hallmark_intervals: List[Tuple[int, int]],
    cfg: dict,
    region_seq: Optional[str] = None,
) -> Tuple[Optional[TirPair], dict]:
    """Filter and rank TIR pairs. Returns (best_tir, diagnostics_dict)."""
    bracket_fraction = cfg.get("tir_bracket_fraction", 0.5)
    require_bracket  = bool(hallmark_intervals)
    n_hallmarks_total = len(hallmark_intervals)
    min_required = (
        max(1, int(np.ceil(n_hallmarks_total * bracket_fraction)))
        if require_bracket else 0
    )

    min_entropy        = cfg.get("tir_min_dinuc_entropy", 2.0)
    max_kmer_frac      = cfg.get("tir_max_kmer_fraction", 0.70)
    max_tandem_frac    = cfg.get("tir_max_tandem_fraction", 0.70)
    max_tandem_period  = cfg.get("tir_tandem_max_period", 12)

    diag = dict(
        n_raw=len(pairs),
        n_pass_insert=0, n_pass_len=0, n_pass_id=0,
        n_pass_complexity=0, n_pass_bracket=0,
        best_near_miss="",
    )

    valid: List[Tuple[int, TirPair]] = []  
    near_miss_score = -1.0

    for t in pairs:
        gl_start = t.left_start  + region_offset - 1
        gl_end   = t.left_end    + region_offset - 1
        gr_start = t.right_start + region_offset - 1
        gr_end   = t.right_end   + region_offset - 1

        ok_insert = cfg["tir_min_insert"] <= t.insert_size <= cfg["tir_max_insert"]
        ok_len    = cfg["tir_min_len"]    <= t.tir_length  <= cfg["tir_max_len"]
        ok_id     = t.tir_identity >= cfg["tir_min_id"]

        ok_complexity = True
        complexity_reason = "ok"
        if region_seq is not None:
            left_seq  = region_seq[t.left_start  - 1: t.left_end]
            right_seq = region_seq[t.right_start - 1: t.right_end]
            is_lc_l, reason_l = _tir_is_low_complexity(
                left_seq, min_entropy, max_kmer_frac, max_tandem_frac, max_tandem_period
            )
            is_lc_r, reason_r = _tir_is_low_complexity(
                right_seq, min_entropy, max_kmer_frac, max_tandem_frac, max_tandem_period
            )
            if is_lc_l or is_lc_r:
                ok_complexity = False
                complexity_reason = f"L:{reason_l or 'ok'} R:{reason_r or 'ok'}"

        abs_tir = TirPair(
            left_start=gl_start, left_end=gl_end,
            right_start=gr_start, right_end=gr_end,
            tir_length=t.tir_length, insert_size=t.insert_size,
            tir_identity=t.tir_identity, score=t.score,
            matches=t.matches, total=t.total, gaps=t.gaps,
            tir_evalue=t.tir_evalue,
        )

        if require_bracket:
            n_bracketed = _count_bracketed(abs_tir, hallmark_intervals)
            ok_bracket  = n_bracketed >= min_required
        else:
            n_bracketed = 0
            ok_bracket  = True

        if ok_insert:
            diag["n_pass_insert"] += 1
        if ok_insert and ok_len:
            diag["n_pass_len"] += 1
        if ok_insert and ok_len and ok_id:
            diag["n_pass_id"] += 1
        if ok_insert and ok_len and ok_id and ok_complexity:
            diag["n_pass_complexity"] += 1
        if ok_insert and ok_len and ok_id and ok_complexity and ok_bracket:
            diag["n_pass_bracket"] += 1

        score = (sum([ok_insert, ok_len, ok_id, ok_complexity, ok_bracket]) * 1e6
                 + t.tir_identity * t.tir_length)
        if score > near_miss_score:
            near_miss_score = score
            diag["best_near_miss"] = (
                f"insert={t.insert_size:,} tir_len={t.tir_length} "
                f"id={t.tir_identity:.1f}% "
                f"bracketed={n_bracketed}/{n_hallmarks_total} "
                f"(need>={min_required}) "
                f"complexity={complexity_reason} "
                f"passed=[insert:{ok_insert},len:{ok_len},id:{ok_id},"
                f"complexity:{ok_complexity},bracket:{ok_bracket}]"
            )
        if not (ok_insert and ok_len and ok_id and ok_complexity and ok_bracket):
            continue

        valid.append((n_bracketed, abs_tir))

    if not valid:
        return None, diag

    valid.sort(
        key=lambda x: (x[0], x[1].tir_identity, x[1].insert_size, x[1].tir_length),
        reverse=True,
    )
    best = valid[0][1]
    _LOG.debug(
        f"TIR selected: insert={best.insert_size:,} bp, "
        f"tir_len={best.tir_length} bp, id={best.tir_identity:.1f}%, "
        f"bracketed={valid[0][0]}/{n_hallmarks_total} hallmarks"
    )
    return best, diag

# Stage 6: TSD detection
def find_tsd(
    left_flank: str,
    right_flank: str,
    k_min: int,
    k_max: int,
    max_slide: int,
) -> Optional[Tsd]:
    """Search for Target Site Duplications flanking the TIR boundaries.
    """
    left  = left_flank.upper()
    right = right_flank.upper()
    best: Optional[Tsd] = None
    for k in range(k_max, k_min - 1, -1):
        if k > len(left) or k > len(right):
            continue
        max_mm = 0 if k <= 5 else (1 if k <= 8 else 2)
        for sl in range(max_slide + 1):
            if k + sl > len(left):
                break
            for sr in range(max_slide + 1):
                if k + sr > len(right):
                    break
                lk = left[len(left) - k - sl: len(left) - sl] if sl > 0 else left[-k:]
                rk = right[sr: sr + k]
                if "N" in lk or "N" in rk:
                    continue
                mm = sum(1 for a, b in zip(lk, rk) if a != b)
                if mm <= max_mm:
                    cand = Tsd(
                        sequence_left=lk, sequence_right=rk,
                        length=k, mismatches=mm,
                        identity=100.0 * (k - mm) / k,
                        left_shift=sl, right_shift=sr,
                    )
                    if best is None or (cand.length, cand.identity) > (best.length, best.identity):
                        best = cand
        if best is not None and best.length == k:
            return best
    return best

# Stage 7: GC composition
def gc_of_seq(seq: str) -> float:
    """Mean GC% (ignoring N) for a sequence."""
    arr = np.frombuffer(seq.upper().encode("ascii"), dtype=np.uint8)
    gc = int(((arr == ord("G")) | (arr == ord("C"))).sum())
    at = int(((arr == ord("A")) | (arr == ord("T"))).sum())
    valid = gc + at
    return float(100.0 * gc / valid) if valid > 0 else float("nan")

# Stage 8: per-cluster worker (TIR + TSD + GC)
def extend_tirless_boundaries(
    contig_orfs: List[Orf],
    cluster_start: int,
    cluster_end: int,
    rolling_by_orf: Dict[str, float],
    cfg: dict,
) -> Tuple[int, int, dict]:
    """Walk the rolling viral score outward from a seed cluster to find GEVE boundaries."""
    threshold       = cfg.get("extend_threshold", -0.5)
    start_threshold = cfg.get("extend_start_threshold", 0.0)
    max_bp    = cfg.get("extend_max_bp", 200_000)
    max_drops = cfg.get("extend_max_drops", 2)
    max_span  = cfg.get("max_cluster_span", 2_000_000)

    diag = dict(
        applied=False,
        n_left_added=0, n_right_added=0,
        extended_left_bp=0, extended_right_bp=0,
        threshold=threshold, start_threshold=start_threshold, max_bp=max_bp,
        stop_left="not_started", stop_right="not_started",
    )

    sorted_orfs = sorted(contig_orfs, key=lambda x: x.start)
    if not sorted_orfs:
        return cluster_start, cluster_end, diag

    left_idx: Optional[int] = None
    right_idx: Optional[int] = None
    for i, o in enumerate(sorted_orfs):
        if o.start >= cluster_start and o.end <= cluster_end:
            if left_idx is None:
                left_idx = i
            right_idx = i
    if left_idx is None or right_idx is None:
        return cluster_start, cluster_end, diag

    edge_left_score  = rolling_by_orf.get(sorted_orfs[left_idx].orf_id,  0.0)
    edge_right_score = rolling_by_orf.get(sorted_orfs[right_idx].orf_id, 0.0)

    new_start = cluster_start
    if edge_left_score <= start_threshold:
        diag["stop_left"] = "edge_below_start_threshold"
    else:
        drops = 0
        i = left_idx - 1
        diag["stop_left"] = "contig_edge"
        while i >= 0:
            o = sorted_orfs[i]
            if cluster_start - o.start > max_bp:
                diag["stop_left"] = "extend_max_bp"
                break
            if (cluster_end - o.start + 1) > max_span:
                diag["stop_left"] = "max_cluster_span"
                break
            score = rolling_by_orf.get(o.orf_id, 0.0)
            if score > threshold:
                new_start = o.start
                drops = 0
                diag["n_left_added"] += 1
            else:
                drops += 1
                if drops >= max_drops:
                    diag["stop_left"] = "below_threshold"
                    break
            i -= 1

    new_end = cluster_end
    if edge_right_score <= start_threshold:
        diag["stop_right"] = "edge_below_start_threshold"
    else:
        drops = 0
        i = right_idx + 1
        diag["stop_right"] = "contig_edge"
        while i < len(sorted_orfs):
            o = sorted_orfs[i]
            if o.end - cluster_end > max_bp:
                diag["stop_right"] = "extend_max_bp"
                break
            if (o.end - new_start + 1) > max_span:
                diag["stop_right"] = "max_cluster_span"
                break
            score = rolling_by_orf.get(o.orf_id, 0.0)
            if score > threshold:
                new_end = o.end
                drops = 0
                diag["n_right_added"] += 1
            else:
                drops += 1
                if drops >= max_drops:
                    diag["stop_right"] = "below_threshold"
                    break
            i += 1

    diag["extended_left_bp"]  = cluster_start - new_start
    diag["extended_right_bp"] = new_end - cluster_end
    diag["applied"] = (new_start < cluster_start) or (new_end > cluster_end)
    return new_start, new_end, diag

def prescan_and_merge_clusters(
    clusters: List[dict],
    orfs_by_contig: Dict[str, List[Orf]],
    rolling_by_orf_per_contig: Dict[str, Dict[str, float]],
    cfg: dict,
) -> List[dict]:
    """Run viral-score pre-scan on each seed cluster, then merge same-contig
    clusters whose extended boundaries overlap or are within cluster_merge_gap.
    """
    extend    = cfg.get("extend_tirless", True)
    merge_gap = cfg.get("cluster_merge_gap", 50_000)
    max_span  = cfg.get("max_cluster_span", 2_000_000)
    host_fraction = cfg.get("host_territory_fraction", 0.7)

    extended: List[dict] = []
    for cl in clusters:
        contig      = cl["contig"]
        contig_orfs = orfs_by_contig.get(contig, [])
        rolling     = rolling_by_orf_per_contig.get(contig, {})
        if extend:
            pre_s, pre_e, diag = extend_tirless_boundaries(
                contig_orfs, cl["cluster_start"], cl["cluster_end"], rolling, cfg,
            )
            if diag.get("applied"):
                _LOG.info(
                    f"Pre-scan ({contig}:{cl['cluster_start']:,}-{cl['cluster_end']:,}): "
                    f"seed -> {pre_s:,}-{pre_e:,} "
                    f"(+{diag['extended_left_bp']:,} bp / "
                    f"+{diag['n_left_added']} ORFs left; "
                    f"+{diag['extended_right_bp']:,} bp / "
                    f"+{diag['n_right_added']} ORFs right; "
                    f"stop_left={diag['stop_left']}, stop_right={diag['stop_right']})"
                )
        else:
            pre_s, pre_e, diag = cl["cluster_start"], cl["cluster_end"], None
        out = dict(cl)
        out["pre_start"]    = pre_s
        out["pre_end"]      = pre_e
        out["prescan_diag"] = diag
        extended.append(out)

    extended.sort(key=lambda c: (c["contig"], c["pre_start"]))
    merged: List[dict] = []
    n_merged = 0
    n_refused = 0
    for c in extended:
        if merged and merged[-1]["contig"] == c["contig"]:
            prev = merged[-1]
            gap = c["pre_start"] - prev["pre_end"] - 1
            new_span = max(prev["pre_end"], c["pre_end"]) - prev["pre_start"] + 1
            if gap <= merge_gap and new_span <= max_span:
                if gap <= 0 or not _gap_is_host_territory(
                        prev["pre_end"], c["pre_start"],
                        orfs_by_contig.get(c["contig"], []),
                        rolling_by_orf_per_contig.get(c["contig"], {}),
                        host_fraction):
                    prev["pre_end"]       = max(prev["pre_end"], c["pre_end"])
                    prev["cluster_start"] = min(prev["cluster_start"], c["cluster_start"])
                    prev["cluster_end"]   = max(prev["cluster_end"], c["cluster_end"])
                    prev["orf_ids"]   = sorted(set(prev["orf_ids"])   | set(c["orf_ids"]))
                    prev["hallmarks"] = sorted(set(prev["hallmarks"]) | set(c["hallmarks"]))
                    prev["prescan_diag"] = dict(prev.get("prescan_diag") or {})
                    prev["prescan_diag"]["applied"] = True
                    n_merged += 1
                    continue
                else:
                    n_refused += 1
        merged.append(dict(c))

    _LOG.info(
        f"Post-extension merge: {len(extended)} -> {len(merged)} cluster(s) "
        f"({n_merged} merged within {merge_gap:,} bp; "
        f"{n_refused} refused as host territory)"
    )
    return merged

def _process_cluster(task: dict) -> dict:
    """Process one hallmark cluster end-to-end inside a worker process."""
    ci             = task["cluster_index"]
    clabel         = task.get("cluster_label", ci)
    contig         = task["contig"]
    cstart         = task["cluster_start"]
    cend           = task["cluster_end"]
    clen           = task["contig_length"]
    cfg            = task["cfg"]
    genome_path    = task["genome_path"]
    contig_orfs    = task["contig_orfs"]
    blastn_threads = int(task.get("blastn_threads", 1))

    log_msgs: List[Tuple[str, str]] = []
    fa = pyfastx.Fasta(str(genome_path), build_index=True, uppercase=True)

    pre_start    = task["pre_start"]
    pre_end      = task["pre_end"]
    prescan_diag = task.get("prescan_diag")

    hallmark_intervals: List[Tuple[int, int]] = [
        (o.start, o.end)
        for o in contig_orfs
        if o.hallmark is not None and o.start <= pre_end and o.end >= pre_start
    ]

    best_tir: Optional[TirPair] = None
    tir_note  = ""
    flank_used = 0
    flanks_tried: List[int] = []

    flank_start = cfg["tir_flank_start"]
    flank_max   = cfg["tir_flank_max"]
    step        = cfg["tir_flank_step"]

    with tempfile.TemporaryDirectory(prefix="findGEVE_") as tmp:
        tmp = Path(tmp)
        rstart_max = max(1, pre_start - flank_max)
        rend_max   = min(clen, pre_end + flank_max)
        region_seq_max = fa.fetch(contig, (rstart_max, rend_max))

        region_fa = tmp / f"cluster_{ci:04d}.fa"
        tab_out   = tmp / f"cluster_{ci:04d}.tsv"
        with open(region_fa, "w") as fh:
            fh.write(f">cluster_{ci:04d}\n")
            for j in range(0, len(region_seq_max), 80):
                fh.write(region_seq_max[j:j + 80] + "\n")

        try:
            run_blastn_self(region_fa, tab_out, cfg, threads=blastn_threads)
        except FileNotFoundError:
            return dict(
                status="fatal", cluster_index=ci, contig=contig,
                message="blastn not found in PATH inside worker",
            )
        except subprocess.CalledProcessError as exc:
            err = (exc.stderr or "").strip()[:120]
            log_msgs.append(("info",
                f"Cluster {clabel} ({contig}:{cstart:,}-{cend:,}): "
                f"blastn failed: {err}; retaining TIR-less "
                f"(boundary={pre_start:,}-{pre_end:,})"
            ))
            all_pairs: List[TirPair] = []
        else:
            all_pairs = parse_blastn_tabular(tab_out)

        last_note = ""
        flank = flank_start
        while flank <= flank_max:
            flanks_tried.append(flank)
            rstart = max(1, pre_start - flank)
            rend   = min(clen, pre_end + flank)
            visible_pairs = [
                p for p in all_pairs
                if (p.left_start  + rstart_max - 1) >= rstart
                and (p.right_end + rstart_max - 1) <= rend
            ]
            best, diag = select_best_tir(
                visible_pairs, region_offset=rstart_max,
                hallmark_intervals=hallmark_intervals, cfg=cfg,
                region_seq=region_seq_max,
            )
            if best is not None:
                best_tir   = best
                flank_used = flank
                break
            if diag["n_raw"] == 0:
                last_note = f"no inverted-repeat pairs in region (flank={flank:,} bp)"
            else:
                last_note = (
                    f"all TIR candidates filtered at flank={flank:,} bp "
                    f"({diag['n_raw']} raw, "
                    f"{diag['n_pass_insert']} insert, "
                    f"{diag['n_pass_len']} len, "
                    f"{diag['n_pass_id']} id, "
                    f"{diag['n_pass_complexity']} complexity, "
                    f"{diag['n_pass_bracket']} bracket); "
                    f"best near-miss: {diag['best_near_miss']}"
                )
            if rstart == 1 and rend == clen:
                break
            flank += step

        if best_tir is None:
            tir_note = (last_note or "no TIR detected") + \
                       f"; flanks tried: {','.join(f'{f:,}' for f in flanks_tried)}"
            log_msgs.append(("info",
                f"Cluster {clabel} ({contig}:{cstart:,}-{cend:,}): "
                f"{tir_note}; retaining TIR-less "
                f"(boundary={pre_start:,}-{pre_end:,})"
            ))

    if best_tir is not None:
        geve_start      = best_tir.left_start
        geve_end        = best_tir.right_end
        boundary_method = "TIR"
    else:
        geve_start = pre_start
        geve_end   = pre_end
        boundary_method = (
            "viral_score_boundary"
            if (prescan_diag is not None and prescan_diag.get("applied"))
            else "seed_cluster"
        )

    geve_length = geve_end - geve_start + 1

    if geve_length < cfg["min_geve_length"]:
        return dict(
            status="skip", cluster_index=ci, contig=contig, log_msgs=log_msgs,
            message=(
                f"Cluster {clabel} ({contig}:{cstart:,}-{cend:,}): span "
                f"{geve_length:,} bp < {cfg['min_geve_length']:,} bp threshold; discarded"
            ),
        )

    geve_orfs = [o for o in contig_orfs
                 if o.start >= geve_start and o.end <= geve_end]
    geve_orfs.sort(key=lambda x: x.start)
    hallmarks_in_geve = [o.hallmark for o in geve_orfs if o.hallmark]
    hallmark_types_in_geve = sorted(set(hallmarks_in_geve))
    if not hallmark_types_in_geve:
        return dict(
            status="skip", cluster_index=ci, contig=contig, log_msgs=log_msgs,
            message=(
                f"Cluster {clabel} ({contig}:{cstart:,}-{cend:,}): "
                f"boundary excludes all hallmark ORFs; "
                f"boundary_method={boundary_method}; discarded"
            ),
        )

    # TSD detection
    tsd: Optional[Tsd] = None
    if best_tir is not None:
        win = cfg["tsd_search_window"]
        l_flank_end   = best_tir.left_start - 1
        l_flank_start = max(1, l_flank_end - win)
        r_flank_start = best_tir.right_end + 1
        r_flank_end   = min(clen, r_flank_start + win)
        left_flank  = fa.fetch(contig, (l_flank_start, l_flank_end)) if l_flank_end >= l_flank_start else ""
        right_flank = fa.fetch(contig, (r_flank_start, r_flank_end)) if r_flank_end >= r_flank_start else ""
        tsd = find_tsd(
            left_flank, right_flank,
            cfg["tsd_min"], cfg["tsd_max"], cfg["tsd_max_slide"],
        )

    # GC content of the GEVE
    geve_seq = fa.fetch(contig, (geve_start, geve_end))
    n_count = geve_seq.upper().count("N")
    n_fraction = n_count / len(geve_seq) if geve_seq else 0.0
    n_max = cfg.get("n_max_fraction", 0.05)
    if n_fraction > n_max:
        return dict(
            status="skip", cluster_index=ci, contig=contig, log_msgs=log_msgs,
            message=(
                f"Cluster {clabel} ({contig}:{cstart:,}-{cend:,}): "
                f"N fraction {n_fraction:.2%} > {n_max:.0%} "
                f"(likely contig-junction artifact); discarded"
            ),
        )
    gc_geve  = gc_of_seq(geve_seq)

    geve_obj = dict(
        geve_id="TBD",
        contig=contig, contig_length=clen,
        geve_start=geve_start, geve_end=geve_end, geve_length=geve_length,
        tir=best_tir, tsd=tsd,
        orfs=geve_orfs,
        n_hallmarks=len(hallmarks_in_geve),
        hallmarks_present=hallmark_types_in_geve,
        gc_geve=gc_geve,
        has_tir=(best_tir is not None),
        flank_used=flank_used,
        boundary_method=boundary_method,
    )
    return dict(
        status="ok", cluster_index=ci, contig=contig,
        log_msgs=log_msgs, geve=geve_obj,
    )

def coding_density(orfs: List[Orf], span_start: int, span_end: int) -> float:
    """Return coding_density (percent) for ORFs in a span."""
    span = span_end - span_start + 1
    if span <= 0 or not orfs:
        return 0.0
    intervals = []
    for o in orfs:
        s = max(o.start, span_start)
        e = min(o.end, span_end)
        if e >= s:
            intervals.append((s, e))
    intervals.sort()
    merged: List[Tuple[int, int]] = []
    for s, e in intervals:
        if merged and s <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    coding_bases = sum(e - s + 1 for s, e in merged)
    pct = min(100.0 * coding_bases / span, 100.0)
    return round(pct, 3)

def _build_tirless_merge(
    group: List[dict],
    fa,
    contig_orfs: List[Orf],
) -> dict:
    """Build a merged GEVE without TIR/TSD (fallback when re-detection fails)."""
    group = sorted(group, key=lambda x: x["geve_start"])
    contig    = group[0]["contig"]
    new_start = group[0]["geve_start"]
    new_end   = max(g["geve_end"] for g in group)
    geve_orfs = sorted(
        [o for o in contig_orfs
         if o.start >= new_start and o.end <= new_end],
        key=lambda x: x.start,
    )
    hms = [o.hallmark for o in geve_orfs if o.hallmark]
    return dict(
        geve_id="TBD",
        contig=contig, contig_length=group[0]["contig_length"],
        geve_start=new_start, geve_end=new_end,
        geve_length=new_end - new_start + 1,
        tir=None, tsd=None, orfs=geve_orfs,
        n_hallmarks=len(hms),
        hallmarks_present=sorted(set(hms)),
        gc_geve=gc_of_seq(fa.fetch(contig, (new_start, new_end))),
        has_tir=False, flank_used=0,
        boundary_method="viral_score_boundary",
    )

def _redetect_on_merged_span(
    group: List[dict],
    fa,
    contig_orfs: List[Orf],
    contig_length: int,
    cfg: dict,
    genome_path: Path,
    blastn_threads: int,
) -> Tuple[dict, List[Tuple[str, str]]]:
    """Re-run TIR/TSD detection on the combined span of merged GEVEs.
    """
    group = sorted(group, key=lambda x: x["geve_start"])
    contig    = group[0]["contig"]
    new_start = group[0]["geve_start"]
    new_end   = max(g["geve_end"] for g in group)

    task = dict(
        cluster_index=0,
        cluster_label="merge",
        contig=contig,
        cluster_start=new_start, cluster_end=new_end,
        pre_start=new_start,     pre_end=new_end,
        prescan_diag={"applied": True},
        contig_length=contig_length,
        cfg=cfg,
        genome_path=str(genome_path),
        contig_orfs=contig_orfs,
        blastn_threads=blastn_threads,
    )
    result = _process_cluster(task)
    log_msgs = result.get("log_msgs", []) or []
    if result.get("status") == "ok":
        return result["geve"], log_msgs
    return _build_tirless_merge(group, fa, contig_orfs), log_msgs


def _merge_adjacent_geves(
    geves: List[dict],
    fa,
    orfs_by_contig: Dict[str, List[Orf]],
    rolling_by_orf_per_contig: Dict[str, Dict[str, float]],
    cfg: dict,
    contig_lengths: Dict[str, int],
    genome_path: Path,
    blastn_threads: int,
) -> List[dict]:
    """Iteratively bridge same-contig GEVEs within cluster_merge_gap and
    re-run TIR/TSD on each merged span.
    """
    if not geves:
        return geves
    merge_gap = cfg.get("cluster_merge_gap", 200_000)
    max_span  = cfg.get("max_cluster_span", 2_000_000)
    host_frac = cfg.get("host_territory_fraction", 0.7)
    n_initial = len(geves)

    iteration = 0
    total_refused = 0
    while True:
        iteration += 1
        by_contig: Dict[str, List[dict]] = {}
        for g in geves:
            by_contig.setdefault(g["contig"], []).append(g)

        next_geves: List[dict] = []
        any_merge = False
        n_refused_pass = 0

        for contig, lst in by_contig.items():
            lst.sort(key=lambda x: x["geve_start"])
            contig_orfs = orfs_by_contig.get(contig, [])
            rolling     = rolling_by_orf_per_contig.get(contig, {})

            group: List[dict] = [lst[0]]
            group_end = lst[0]["geve_end"]

            def _flush(group_local):
                """Emit the current group: re-detect TIR if 2+ pieces, else pass through."""
                nonlocal any_merge
                if len(group_local) > 1:
                    merged, lmsgs = _redetect_on_merged_span(
                        group_local, fa, contig_orfs, contig_lengths[contig],
                        cfg, genome_path, blastn_threads,
                    )
                    for level, msg in lmsgs:
                        (_LOG.warning if level == "warning" else _LOG.info)(msg)
                    tir_str = (
                        f"TIR={merged['tir'].tir_length} bp @ "
                        f"{merged['tir'].tir_identity:.1f}%"
                        if merged.get("tir") is not None else "TIR=none"
                    )
                    _LOG.info(
                        f"Adjacent-GEVE merge iter {iteration}: {len(group_local)} pieces on "
                        f"{contig} bridged into {contig}:{merged['geve_start']:,}-"
                        f"{merged['geve_end']:,} ({merged['geve_length']:,} bp, "
                        f"boundary={merged.get('boundary_method', 'NA')}, {tir_str})"
                    )
                    next_geves.append(merged)
                    any_merge = True
                else:
                    next_geves.append(group_local[0])

            for nxt in lst[1:]:
                gap = nxt["geve_start"] - group_end - 1
                proposed_end  = max(group_end, nxt["geve_end"])
                proposed_span = proposed_end - group[0]["geve_start"] + 1
                is_overlap = gap < 0
                eligible = is_overlap or (
                    (gap <= merge_gap) and (proposed_span <= max_span)
                )
                if eligible and (is_overlap or gap <= 0 or not _gap_is_host_territory(
                        group_end, nxt["geve_start"], contig_orfs, rolling, host_frac)):
                    group.append(nxt)
                    group_end = proposed_end
                    continue
                if eligible:
                    n_refused_pass += 1
                _flush(group)
                group = [nxt]
                group_end = nxt["geve_end"]

            _flush(group)

        total_refused += n_refused_pass
        _LOG.info(
            f"Adjacent-GEVE merge iter {iteration}: "
            f"{len(geves)} -> {len(next_geves)} GEVE(s) "
            f"({n_refused_pass} refused as host territory)"
        )

        if not any_merge:
            break
        geves = next_geves

    if iteration > 1 or len(next_geves) != n_initial:
        _LOG.info(
            f"Adjacent-GEVE merge: converged after {iteration} iteration(s); "
            f"{n_initial} -> {len(next_geves)} GEVE(s) "
            f"({total_refused} total refused as host territory)"
        )
    return next_geves


def _resolve_overlapping_geves(
    geves: List[dict],
    fa,
    orfs_by_contig: Dict[str, List[Orf]],
) -> List[dict]:
    """Merge GEVEs whose spans overlap on the same contig.
    """
    if not geves:
        return geves
    by_contig: Dict[str, List[dict]] = {}
    for g in geves:
        by_contig.setdefault(g["contig"], []).append(g)

    resolved: List[dict] = []
    merge_events: List[Tuple[str, int, int, int]] = []
    for contig, lst in by_contig.items():
        lst.sort(key=lambda x: x["geve_start"])
        cur = dict(lst[0])
        merged_count = 1
        for nxt in lst[1:]:
            if nxt["geve_start"] <= cur["geve_end"]:
                new_start = min(cur["geve_start"], nxt["geve_start"])
                new_end   = max(cur["geve_end"],   nxt["geve_end"])
                cur["geve_start"]     = new_start
                cur["geve_end"]       = new_end
                cur["geve_length"]    = new_end - new_start + 1
                cur["tir"]            = None
                cur["tsd"]            = None
                cur["has_tir"]        = False
                cur["flank_used"]     = 0
                cur["boundary_method"] = "viral_score_boundary"
                contig_orfs = orfs_by_contig.get(contig, [])
                cur["orfs"] = sorted(
                    [o for o in contig_orfs
                     if o.start >= new_start and o.end <= new_end],
                    key=lambda x: x.start,
                )
                hms = [o.hallmark for o in cur["orfs"] if o.hallmark]
                cur["n_hallmarks"]       = len(hms)
                cur["hallmarks_present"] = sorted(set(hms))
                cur["gc_geve"] = gc_of_seq(fa.fetch(contig, (new_start, new_end)))
                merged_count += 1
                continue
            if merged_count > 1:
                merge_events.append((contig, cur["geve_start"], cur["geve_end"], merged_count))
            resolved.append(cur)
            cur = dict(nxt)
            merged_count = 1
        if merged_count > 1:
            merge_events.append((contig, cur["geve_start"], cur["geve_end"], merged_count))
        resolved.append(cur)

    if merge_events:
        for contig, gs, ge, n in merge_events:
            _LOG.info(
                f"Overlap resolution: {n} overlapping GEVE(s) on {contig} merged into "
                f"{contig}:{gs:,}-{ge:,} ({ge - gs + 1:,} bp, boundary=viral_score_boundary)"
            )
        _LOG.info(
            f"Overlap resolution: {len(geves)} -> {len(resolved)} non-overlapping GEVE(s)"
        )
    return resolved

# Output writers
def _get_tool_versions() -> Dict[str, str]:
    versions: Dict[str, str] = {}
    for pkg in ("pyrodigal-gv", "pyhmmer", "pyfastx", "numpy", "pandas"):
        try:
            versions[pkg] = _pkg_version(pkg)
        except PackageNotFoundError:
            versions[pkg] = "unknown"
    try:
        out = subprocess.run(
            ["blastn", "-version"], capture_output=True, text=True, check=False,
        )
        first = (out.stdout or out.stderr).strip().splitlines()
        versions["blastn"] = first[0] if first else "unknown"
    except FileNotFoundError:
        versions["blastn"] = "not_found"
    versions["python"] = sys.version.split()[0]
    return versions

def _tsd_fields(tsd: Optional[Tsd]) -> dict:
    if tsd is None:
        return dict(
            tsd_len="NA", tsd_left="NA", tsd_right="NA",
            tsd_mismatch="NA", tsd_conservation="NODETECT",
        )
    return dict(
        tsd_len=tsd.length,
        tsd_left=tsd.sequence_left,
        tsd_right=tsd.sequence_right,
        tsd_mismatch=tsd.mismatches,
        tsd_conservation="PERFECT" if tsd.mismatches == 0 else "IMPERFECT",
    )

def _tir_fields(tir: Optional[TirPair]) -> dict:
    if tir is None:
        return dict(
            tir_length="NA", tir_score="NA",
            tir_identity_pct="NA", tir_gaps="NA",
        )
    return dict(
        tir_length=tir.tir_length,
        tir_score=tir.score,
        tir_identity_pct=round(tir.tir_identity, 2),
        tir_gaps=tir.gaps,
    )

def load_gvog_annotations(db: Path) -> Dict[str, str]:
    """Load GVOG name lookup from gvog.complete.annot.tsv if present.
    """
    path = db / "gvog.complete.annot.tsv"
    if not path.exists():
        _LOG.warning(
            f"GVOG annotation file not found: {path}; "
            f"gvog_name will be NA in {{prefix}}.geve.func.tsv"
        )
        return {}
    try:
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except Exception as exc:
        _LOG.warning(f"Failed to load GVOG annotations from {path}: {exc}")
        return {}
    if "GVOG" not in df.columns or "NCVOG_descs" not in df.columns:
        _LOG.warning(
            f"Expected columns 'GVOG' and 'NCVOG_descs' missing in "
            f"{path.name}; gvog_name will be NA"
        )
        return {}
    name_map: Dict[str, str] = {}
    for _, row in df.iterrows():
        gid = (row.get("GVOG") or "").strip()
        if not gid:
            continue
        descs = (row.get("NCVOG_descs") or "").strip()
        if not descs:
            continue

        first = descs.split(" | ")[0].strip()
        if first:
            name_map[gid] = first
    _LOG.info(
        f"Loaded {len(name_map):,} GVOG name annotations from {path.name}"
    )
    return name_map


def write_protein_func_tsv(
    geves: List[dict],
    path: Path,
    gvog_name_map: Dict[str, str],
) -> None:
    """Write per-protein functional annotation table for retained GEVEs.
    """
    def _fmt_score(v: float) -> str:
        return f"{v:.1f}" if v and v > 0 else "NA"

    def _fmt_evalue(v: float) -> str:
        if v is None or v == float("inf") or not np.isfinite(v):
            return "NA"
        return f"{v:.2e}"

    def _fmt_str(s: Optional[str]) -> str:
        return s if s else "NA"

    rows = []
    for g in geves:
        geve_id = g["geve_id"]
        for orf_idx, o in enumerate(g["orfs"], start=1):
            if o.hallmark is None and o.gvog is None and o.best_pfam_acc is None:
                continue
            label     = f"orf{orf_idx:05d}"
            gvog_name = gvog_name_map.get(o.gvog, "") if o.gvog else ""
            rows.append(dict(
                geve_id           = geve_id,
                protein_id        = label,
                gvog_bitscore     = _fmt_score(o.gvog_bitscore),
                gvog_evalue       = _fmt_evalue(o.gvog_evalue),
                gvog_id           = _fmt_str(o.gvog),
                gvog_name         = _fmt_str(gvog_name),
                pfam_bitscore     = _fmt_score(o.best_pfam_bitscore),
                pfam_evalue       = _fmt_evalue(o.best_pfam_evalue),
                pfam_id           = _fmt_str(o.best_pfam_acc),
                pfam_name         = _fmt_str(o.best_pfam_name),
            ))
    df = pd.DataFrame(rows)
    df.to_csv(path, sep="\t", index=False)


def write_summary_tsv(geves: List[dict], path: Path) -> None:
    rows = []
    for g in geves:
        orfs = g["orfs"]
        cod = coding_density(orfs, g["geve_start"], g["geve_end"])
        row = dict(
            contig_id      = g["contig"],
            geve_name      = g["geve_id"],
            start          = g["geve_start"],
            end            = g["geve_end"],
            geve_length    = g["geve_length"],
            gc             = round(g["gc_geve"], 2),
            total_cds      = len(orfs),
            NCLDV_hits     = sum(1 for o in orfs if o.hallmark or o.gvog),
            coding_density = cod,
            n_hallmarks    = g["n_hallmarks"],
            hallmarks      = ",".join(g["hallmarks_present"]),
            has_tir        = "yes" if g["has_tir"] else "no",
        )
        row.update(_tir_fields(g["tir"]))
        row.update(_tsd_fields(g["tsd"]))
        rows.append(row)
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(
            by="geve_name",
            key=lambda col: col.map(_natural_key),
        ).reset_index(drop=True)
    df.to_csv(path, sep="\t", index=False)

def write_markerout(geves: List[dict], path: Path) -> None:
    with open(path, "w") as fh:
        fh.write("contig\tgeve_name\tfeature\tstart\tend\tstrand\te_value\tscore\n")
        for g in geves:
            geve_name = g["geve_id"]
            contig = g["contig"]
            for o in g["orfs"]:
                if o.hallmark:
                    strand = "+" if o.strand >= 0 else "-"
                    fh.write(
                        f"{contig}\t{geve_name}\t{o.hallmark}\t{o.start}\t{o.end}\t"
                        f"{strand}\t{o.hallmark_evalue:.2e}\t{o.hallmark_bitscore:.1f}\n"
                    )
            tir = g["tir"]
            if tir is not None:
                fh.write(
                    f"{contig}\t{geve_name}\tTIR_left\t{tir.left_start}\t{tir.left_end}\t"
                    f"+\tNA\t{tir.score}\n"
                )
                fh.write(
                    f"{contig}\t{geve_name}\tTIR_right\t{tir.right_start}\t{tir.right_end}\t"
                    f"-\tNA\t{tir.score}\n"
                )
            tsd = g["tsd"]
            if tsd is not None and tir is not None:
                ltsd_end   = tir.left_start - 1 - tsd.left_shift
                ltsd_start = ltsd_end - tsd.length + 1
                rtsd_start = tir.right_end + 1 + tsd.right_shift
                rtsd_end   = rtsd_start + tsd.length - 1
                fh.write(
                    f"{contig}\t{geve_name}\tTSD_5p\t{ltsd_start}\t{ltsd_end}\t"
                    f"+\tNA\t{tsd.identity:.1f}\n"
                )
                fh.write(
                    f"{contig}\t{geve_name}\tTSD_3p\t{rtsd_start}\t{rtsd_end}\t"
                    f"+\tNA\t{tsd.identity:.1f}\n"
                )

def write_multifasta(
    geves: List[dict],
    fa_index: pyfastx.Fasta,
    path: Path,
) -> None:
    with open(path, "w") as fh:
        for g in geves:
            seq = fa_index.fetch(g["contig"], (g["geve_start"], g["geve_end"]))
            tir = g["tir"]
            tsd = g["tsd"]
            header = (
                f">{g['geve_id']} "
                f"contig={g['contig']} "
                f"start={g['geve_start']} "
                f"end={g['geve_end']} "
                f"length={g['geve_length']} "
                f"hallmarks={','.join(g['hallmarks_present'])} "
                f"gc={g['gc_geve']:.2f}%"
            )
            if tir is not None:
                header += (
                    f" tirL={tir.left_start}..{tir.left_end}"
                    f" tirR={tir.right_start}..{tir.right_end}"
                    f" tir_id={tir.tir_identity:.1f}%"
                )
            if tsd is not None:
                header += f" tsd={tsd.sequence_left}|{tsd.sequence_right}"
            fh.write(header + "\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")

def write_protein_fasta(geves: List[dict], path: Path) -> None:
    with open(path, "w") as fh:
        for g in geves:
            geve_name = g["geve_id"]
            for orf_idx, o in enumerate(g["orfs"], start=1):
                label = f"orf{orf_idx:05d}"
                if o.hallmark:
                    header = f">{geve_name}_{label} {o.hallmark} length={len(o.protein)}"
                else:
                    header = f">{geve_name}_{label} length={len(o.protein)}"
                fh.write(header + "\n")
                seq = o.protein
                for i in range(0, len(seq), 80):
                    fh.write(seq[i:i + 80] + "\n")

_REVCOMP_TABLE = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

def _revcomp(seq: str) -> str:
    return seq.translate(_REVCOMP_TABLE)[::-1]

def write_cds_fasta(geves: List[dict], fa_index: pyfastx.Fasta, path: Path) -> None:
    with open(path, "w") as fh:
        for g in geves:
            geve_name = g["geve_id"]
            for orf_idx, o in enumerate(g["orfs"], start=1):
                label = f"orf{orf_idx:05d}"
                seq = fa_index.fetch(o.contig, (o.start, o.end))
                if o.strand < 0:
                    seq = _revcomp(seq)
                if o.hallmark:
                    header = f">{geve_name}_{label} {o.hallmark} length={len(seq)}"
                else:
                    header = f">{geve_name}_{label} length={len(seq)}"
                fh.write(header + "\n")
                for i in range(0, len(seq), 80):
                    fh.write(seq[i:i + 80] + "\n")

def write_hallmark_peps(geves: List[dict], outdir: Path, prefix: str) -> List[Path]:
    """For each hallmark type present across GEVEs, write a `<prefix>.<hallmark>.pep`
    file containing one entry per GEVE (longest copy if duplicated)."""
    by_hallmark: Dict[str, Dict[str, Orf]] = defaultdict(dict)
    for g in geves:
        geve_id = g["geve_id"]
        for o in g["orfs"]:
            if not o.hallmark:
                continue
            current = by_hallmark[o.hallmark].get(geve_id)
            if current is None or len(o.protein) > len(current.protein):
                by_hallmark[o.hallmark][geve_id] = o

    hallmark_dir = outdir / "hallmark"
    hallmark_dir.mkdir(parents=True, exist_ok=True)

    written: List[Path] = []
    for hallmark, geve_map in by_hallmark.items():
        out_path = hallmark_dir / f"{prefix}.{hallmark.lower()}.pep"
        items = sorted(geve_map.items(), key=lambda kv: _natural_key(kv[0]))
        with open(out_path, "w") as fh:
            for geve_id, o in items:
                fh.write(f">{geve_id}_{hallmark}\n")
                seq = o.protein
                for i in range(0, len(seq), 80):
                    fh.write(seq[i:i + 80] + "\n")
        written.append(out_path)
    return written

def write_gff3(geves: List[dict], path: Path) -> None:
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        for g in geves:
            geve_id = g["geve_id"]
            contig  = g["contig"]
            gstart  = g["geve_start"]
            gend    = g["geve_end"]
            fh.write(
                f"{contig}\tfindGEVE\tmobile_genetic_element\t{gstart}\t{gend}\t.\t+\t.\t"
                f"ID={geve_id};Name={geve_id}\n"
            )
            tir = g["tir"]
            if tir is not None:
                fh.write(
                    f"{contig}\tfindGEVE\tterminal_inverted_repeat\t"
                    f"{tir.left_start}\t{tir.left_end}\t.\t+\t.\t"
                    f"ID={geve_id}.TIR_left;Parent={geve_id}\n"
                )
                fh.write(
                    f"{contig}\tfindGEVE\tterminal_inverted_repeat\t"
                    f"{tir.right_start}\t{tir.right_end}\t.\t-\t.\t"
                    f"ID={geve_id}.TIR_right;Parent={geve_id}\n"
                )
            tsd = g["tsd"]
            if tsd is not None and tir is not None:
                ltsd_end   = tir.left_start - 1 - tsd.left_shift
                ltsd_start = ltsd_end - tsd.length + 1
                rtsd_start = tir.right_end + 1 + tsd.right_shift
                rtsd_end   = rtsd_start + tsd.length - 1
                fh.write(
                    f"{contig}\tfindGEVE\ttarget_site_duplication\t"
                    f"{ltsd_start}\t{ltsd_end}\t.\t+\t.\t"
                    f"ID={geve_id}.TSD_5p;Parent={geve_id}\n"
                )
                fh.write(
                    f"{contig}\tfindGEVE\ttarget_site_duplication\t"
                    f"{rtsd_start}\t{rtsd_end}\t.\t+\t.\t"
                    f"ID={geve_id}.TSD_3p;Parent={geve_id}\n"
                )
            for orf_idx, o in enumerate(g["orfs"], start=1):
                label = f"orf{orf_idx:05d}"
                strand = "+" if o.strand >= 0 else "-"
                if o.hallmark:
                    score_field = f"{o.hallmark_bitscore:.1f}"
                    attrs = (
                        f"ID={geve_id}.{label};Parent={geve_id};"
                        f"Name={label};hallmark={o.hallmark}"
                    )
                else:
                    score_field = "."
                    attrs = f"ID={geve_id}.{label};Parent={geve_id};Name={label}"
                fh.write(
                    f"{contig}\tfindGEVE\tCDS\t{o.start}\t{o.end}\t"
                    f"{score_field}\t{strand}\t0\t{attrs}\n"
                )

def log_run_summary(
    geves: List[dict],
    genome_path: Path,
) -> None:
    n = len(geves)
    n_tir = sum(1 for g in geves if g["has_tir"])
    _LOG.info("Result Summary")
    _LOG.info(f"Input genome: {genome_path}")
    _LOG.info(f"GEVE candidates: {n}")
    if n > 0:
        _LOG.info(f"  With TIR:      {n_tir} ({100 * n_tir / n:.1f}%)")
        _LOG.info(f"  Without TIR:   {n - n_tir} ({100 * (n - n_tir) / n:.1f}%)")
        tir_gevs = [g for g in geves if g["has_tir"]]
        if tir_gevs:
            tsd_counts = dict(PERFECT=0, IMPERFECT=0, NODETECT=0)
            for g in tir_gevs:
                if g["tsd"] is None:
                    tsd_counts["NODETECT"] += 1
                elif g["tsd"].mismatches == 0:
                    tsd_counts["PERFECT"] += 1
                else:
                    tsd_counts["IMPERFECT"] += 1
            _LOG.info("TSD (TIR-bearing GEVEs):")
            for k, v in tsd_counts.items():
                _LOG.info(f"  {k}: {v} ({100 * v / len(tir_gevs):.1f}%)")
        lengths = [g["geve_length"] for g in geves]
        _LOG.info("GEVE length (bp):")
        _LOG.info(f"  Min:    {min(lengths):,}")
        _LOG.info(f"  Max:    {max(lengths):,}")
        _LOG.info(f"  Mean:   {int(np.mean(lengths)):,}")
        _LOG.info(f"  Median: {int(np.median(lengths)):,}")
        hm_counter: Counter = Counter()
        for g in geves:
            for h in g["hallmarks_present"]:
                hm_counter[h] += 1
        top = hm_counter.most_common(6)
        _LOG.info("Most frequent hallmarks: " +
                  ", ".join(f"{h}({c})" for h, c in top))


# Main
class _Parser(argparse.ArgumentParser):
    def format_help(self) -> str:
        return HELP_TEXT

    def format_usage(self) -> str:
        return USAGE_TEXT


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = _Parser(prog="findGEVE.py", add_help=False)

    p.add_argument("-h", "--help", action="help",
                   help="Show this help and exit")

    p.add_argument("genome", type=Path)
    p.add_argument("-db", "--db", type=Path, required=True)
    p.add_argument("--prefix", type=str, required=True)
    p.add_argument("-o", "--outdir", type=Path,
                   default=Path(f"Result_{datetime.now().strftime('%Y%m%d')}"))
    p.add_argument("-t", "--threads", type=int, default=4)
    p.add_argument("-e", "--evalue", type=float, default=DEFAULTS["evalue"])
    p.add_argument("--blastn-jobs", type=int, default=None)
    p.add_argument("-m", "--min-hallmark-type", type=int, default=DEFAULTS["min_hallmarks"])
    p.add_argument("--blastn-threads",       type=int,   default=1)
    p.add_argument("--min-contig",           type=int,   default=DEFAULTS["min_contig"])
    p.add_argument("--cluster-merge-gap",    type=int,   default=DEFAULTS["cluster_merge_gap"])

    return p.parse_args(argv)

def _write_empty_outputs(outdir: Path, prefix: str, genome_path: Path) -> None:
    (outdir / f"{prefix}.geve.fna").write_text("")
    (outdir / f"{prefix}.geve.pep").write_text("")
    (outdir / f"{prefix}.geve.cds").write_text("")
    (outdir / f"{prefix}.markerout").write_text(
        "contig\tgeve_name\tfeature\tstart\tend\tstrand\te_value\tscore\n"
    )
    (outdir / f"{prefix}.geve.gff3").write_text("##gff-version 3\n")
    pd.DataFrame().to_csv(outdir / f"{prefix}.summary.tsv", sep="\t", index=False)
    pd.DataFrame().to_csv(outdir / f"{prefix}.geve.func.tsv", sep="\t", index=False)
    log_run_summary([], genome_path)

def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    cfg = dict(DEFAULTS)
    cfg["min_contig"]           = args.min_contig
    cfg["min_hallmarks"]        = args.min_hallmark_type
    cfg["cluster_merge_gap"]    = args.cluster_merge_gap
    cfg["evalue"]               = args.evalue

    args.outdir.mkdir(parents=True, exist_ok=True)
    log_path = args.outdir / "run.log"
    setup_logging(log_path)

    blastn_threads = max(1, int(args.blastn_threads))
    if args.blastn_jobs is not None:
        blastn_jobs = max(1, int(args.blastn_jobs))
    else:
        blastn_jobs = max(1, args.threads // blastn_threads)

    t0 = time.time()
    _LOG.info(
        f"findGEVE started | prefix='{args.prefix}' | "
        f"threads={args.threads} | blastn_jobs={blastn_jobs} | "
        f"blastn_threads={blastn_threads} | "
        f"genome={args.genome}"
    )
    _LOG.info(
        f"Parameters | min_contig={cfg['min_contig']:,} bp | "
        f"min_geve_length={cfg['min_geve_length']:,} bp | "
        f"min_hallmarks_seed={cfg['min_hallmarks_seed']} | "
        f"min_hallmarks_final={cfg['min_hallmarks']} | "
        f"seed_window={cfg['seed_window']:,} bp | "
        f"cluster_merge_gap={cfg['cluster_merge_gap']:,} bp | "
        f"max_cluster_span={cfg['max_cluster_span']:,} bp | "
        f"tir_flank={cfg['tir_flank_start']:,}->{cfg['tir_flank_max']:,} "
        f"step {cfg['tir_flank_step']:,} bp | "
        f"tir_bracket_fraction={cfg['tir_bracket_fraction']} | "
        f"tir_min_len={cfg['tir_min_len']} bp | "
        f"host_territory_fraction={cfg['host_territory_fraction']}"
    )
    if cfg["extend_tirless"]:
        _LOG.info(
            f"TIR-less extension | enabled | "
            f"start_threshold={cfg['extend_start_threshold']} | "
            f"continue_threshold={cfg['extend_threshold']} | "
            f"max={cfg['extend_max_bp']:,} bp/side | "
            f"max_drops={cfg['extend_max_drops']}"
        )
    else:
        _LOG.info("TIR-less extension | disabled (--no-extend-tirless)")

    if not args.genome.exists():
        _LOG.error(f"Genome file not found: {args.genome}")
        return 2

    db = args.db
    hallmark_hmm = db / "NCLDV_markers.hmm"
    gvog_hmm     = db / "gvog.complete.hmm"
    pfam_hmm     = db / "Pfam-A.hmm"

    for f, label in [(hallmark_hmm, "NCLDV_markers.hmm"), (gvog_hmm, "gvog.complete.hmm")]:
        if not f.exists():
            _LOG.error(f"Required HMM not found: {f} ({label})")
            return 2

    run_pfam = pfam_hmm.exists()
    if not run_pfam:
        _LOG.info("Pfam-A.hmm not found in -db directory; skipping Pfam annotation.")

    if not shutil.which("blastn"):
        _LOG.error(
            "blastn not found in PATH. "
            "Install BLAST+: conda install -c bioconda blast"
        )
        return 2

    versions = _get_tool_versions()
    _LOG.info("Tool versions | " + " | ".join(f"{k}={v}" for k, v in versions.items()))
    _LOG.info(f"Command line | {' '.join(sys.argv)}")

    # Stage 1: ORF prediction
    orfs_by_id, contig_lengths = predict_orfs(args.genome, cfg["min_contig"], args.threads)
    fa = pyfastx.Fasta(str(args.genome), build_index=True, uppercase=True)

    # Build a contig -> ORF list once and reuse throughout
    orfs_by_contig: Dict[str, List[Orf]] = {}
    for o in orfs_by_id.values():
        orfs_by_contig.setdefault(o.contig, []).append(o)

    # Stage 2a: hallmark scan
    contig2hallmark_hits = scan_hallmarks(
        orfs_by_id, hallmark_hmm,
        cfg["evalue"], args.threads,
        score_cutoffs=cfg.get("hallmark_score_cutoffs"),
    )

    hallmark_contigs = set(contig2hallmark_hits.keys())
    if not hallmark_contigs:
        _LOG.warning("No NCLDV hallmark hits detected. No GEVE candidates.")
        _write_empty_outputs(args.outdir, args.prefix, args.genome)
        return 0

    _LOG.info(f"Hallmark-positive contigs: {len(hallmark_contigs):,}")

    # Stage 2b/c: GVOG and (optionally) Pfam scans, scoped to hallmark-positive contigs
    gvog_targets = [o for o in orfs_by_id.values() if o.contig in hallmark_contigs]
    scan_gvog(gvog_targets, gvog_hmm, cfg["evalue"], args.threads)
    if run_pfam:
        scan_pfam(gvog_targets, pfam_hmm, cfg["evalue"], args.threads)

    for o in orfs_by_id.values():
        o.net_score = o.virbit - max(0.0, o.pfambit - o.virbit)

    # Rolling viral score (computed before seeding for the gap merge check)
    rolling_by_orf_per_contig = compute_rolling_scores(
        orfs_by_contig, cfg["rolling_window"], hallmark_contigs,
    )

    # Stage 3: seeding uses min_hallmarks_seed (always 1) — every hallmark anchors a seed
    clusters = find_seed_clusters(
        orfs_by_contig,
        cfg["seed_window"],
        cfg["min_hallmarks_seed"],
        cfg["cluster_merge_gap"],
        cfg["max_cluster_span"],
        rolling_by_orf_per_contig=rolling_by_orf_per_contig,
        host_fraction=cfg["host_territory_fraction"],
    )
    if not clusters:
        _LOG.warning("No clusters passed hallmark-density criterion. No GEVE candidates.")
        _write_empty_outputs(args.outdir, args.prefix, args.genome)
        return 0

    # Stage 3b: viral-score pre-scan + post-extension same-contig merge
    clusters = prescan_and_merge_clusters(
        clusters, orfs_by_contig, rolling_by_orf_per_contig, cfg,
    )

    tasks = []
    for ci, cl in enumerate(clusters, start=1):
        contig = cl["contig"]
        tasks.append(dict(
            cluster_index=ci,
            contig=contig,
            cluster_start=cl["cluster_start"],
            cluster_end=cl["cluster_end"],
            pre_start=cl["pre_start"],
            pre_end=cl["pre_end"],
            prescan_diag=cl.get("prescan_diag"),
            contig_length=contig_lengths[contig],
            cfg=cfg,
            genome_path=str(args.genome),
            contig_orfs=orfs_by_contig.get(contig, []),
            blastn_threads=blastn_threads,
        ))

    if (args.blastn_jobs is None
            and args.blastn_threads == 1
            and len(tasks) < args.threads):
        per_job_threads = max(1, args.threads // max(1, len(tasks)))
        if per_job_threads > 1:
            blastn_threads = per_job_threads
            blastn_jobs = max(1, len(tasks))
            for t in tasks:
                t["blastn_threads"] = blastn_threads
            _LOG.info(
                f"Auto-balance: {len(tasks)} cluster(s) < {args.threads} threads; "
                f"using blastn_jobs={blastn_jobs}, blastn_threads={blastn_threads}"
            )

    n_workers = max(1, min(blastn_jobs, len(tasks)))
    _LOG.info(
        f"TIR/TSD/composition: dispatching {len(tasks)} cluster(s) "
        f"across {n_workers} parallel worker(s) "
        f"(blastn_jobs={blastn_jobs}, blastn_threads={blastn_threads})"
    )

    executor = None
    if n_workers > 1:
        executor = ProcessPoolExecutor(max_workers=n_workers)
        results_iter = executor.map(_process_cluster, tasks, chunksize=1)
    else:
        results_iter = (_process_cluster(t) for t in tasks)

    raw_geves: List[dict] = []
    fatal_msg: Optional[str] = None
    n_done = 0
    try:
        for res in results_iter:
            n_done += 1
            ci = res.get("cluster_index", n_done)
            for level, msg in res.get("log_msgs", []) or []:
                (_LOG.warning if level == "warning" else _LOG.info)(msg)
            status = res.get("status")
            if status == "fatal":
                fatal_msg = res.get("message", "unknown fatal error in worker")
                _LOG.error(fatal_msg)
                break
            if status == "skip":
                _LOG.info(res.get("message", f"Cluster {ci}: discarded"))
                continue
            if status == "ok":
                g = res["geve"]
                raw_geves.append(g)
                tir_str = (
                    f"TIR={g['tir'].tir_length} bp @ {g['tir'].tir_identity:.1f}% "
                    f"(flank={g['flank_used']:,} bp)"
                    if g["tir"] is not None else "TIR=none"
                )
                _LOG.info(
                    f"Cluster {ci}/{len(tasks)} accepted | "
                    f"{g['contig']}:{g['geve_start']:,}-{g['geve_end']:,} | "
                    f"length={g['geve_length']:,} bp | "
                    f"{tir_str} | "
                    f"boundary={g.get('boundary_method', 'NA')} | "
                    f"hallmarks={g['n_hallmarks']} ({','.join(g['hallmarks_present'])}) | "
                    f"TSD={'yes' if g['tsd'] else 'no'}"
                )
    finally:
        if executor is not None:
            executor.shutdown(wait=True)

    if fatal_msg is not None:
        return 2

    # Stage 7.5a: iteratively bridge same-contig GEVEs and re-run TIR/TSD
    raw_geves = _merge_adjacent_geves(
        raw_geves, fa, orfs_by_contig, rolling_by_orf_per_contig, cfg,
        contig_lengths, args.genome, blastn_threads,
    )

    # Stage 7.5b: resolve overlapping GEVE spans (drop spurious TIRs, merge regions)
    raw_geves = _resolve_overlapping_geves(raw_geves, fa, orfs_by_contig)

    # Stage 7.6: final QC — require >= min_hallmarks distinct types on the merged GEVE
    before = len(raw_geves)
    raw_geves = [
        g for g in raw_geves
        if len(g["hallmarks_present"]) >= cfg["min_hallmarks"]
    ]
    n_dropped = before - len(raw_geves)
    if n_dropped:
        _LOG.info(
            f"Final hallmark filter: dropped {n_dropped} GEVE(s) with "
            f"< {cfg['min_hallmarks']} distinct hallmark types"
        )

    # Sort and assign final IDs by natural contig order, then position
    raw_geves.sort(key=lambda g: (_natural_key(g["contig"]), g["geve_start"]))
    for i, g in enumerate(raw_geves, start=1):
        g["geve_id"] = f"{args.prefix}_GEVE_{i:03d}"

    n_tir          = sum(1 for g in raw_geves if g["has_tir"])
    n_score_bound  = sum(1 for g in raw_geves if g.get("boundary_method") == "viral_score_boundary")
    n_seed_only    = sum(1 for g in raw_geves if g.get("boundary_method") == "seed_cluster")
    _LOG.info(
        f"Final: {len(raw_geves)} GEVE candidate(s) "
        f"({n_tir} TIR-bounded, "
        f"{n_score_bound} viral-score boundary, "
        f"{n_seed_only} seed-cluster only)"
    )

    # Stage 8: write outputs
    out = args.outdir
    fasta_path    = out / f"{args.prefix}.geve.fna"
    pep_path      = out / f"{args.prefix}.geve.pep"
    cds_path      = out / f"{args.prefix}.geve.cds"
    marker_path   = out / f"{args.prefix}.markerout"
    summary_path  = out / f"{args.prefix}.summary.tsv"
    gff3_path     = out / f"{args.prefix}.geve.gff3"
    func_path     = out / f"{args.prefix}.geve.func.tsv"

    gvog_name_map = load_gvog_annotations(db)

    write_multifasta(raw_geves, fa, fasta_path)
    write_protein_fasta(raw_geves, pep_path)
    write_cds_fasta(raw_geves, fa, cds_path)
    write_markerout(raw_geves, marker_path)
    write_summary_tsv(raw_geves, summary_path)
    write_gff3(raw_geves, gff3_path)
    write_protein_func_tsv(raw_geves, func_path, gvog_name_map)
    hallmark_pep_paths = write_hallmark_peps(raw_geves, out, args.prefix)

    _LOG.output(f"GEVE sequences  -> {fasta_path}")
    _LOG.output(f"GEVE proteins   -> {pep_path}")
    _LOG.output(f"GEVE CDS        -> {cds_path}")
    _LOG.output(f"Marker annot    -> {marker_path}")
    _LOG.output(f"Summary table   -> {summary_path}")
    _LOG.output(f"GFF3 annotation -> {gff3_path}")
    _LOG.output(f"Functional annot-> {func_path}")
    for hp in hallmark_pep_paths:
        _LOG.output(f"Hallmark pep    -> {hp}")
    _LOG.output(f"Run log         -> {log_path}")

    elapsed = time.time() - t0
    if elapsed >= 3600:
        timing = f"{elapsed / 3600:.2f} h"
    elif elapsed >= 60:
        timing = f"{elapsed / 60:.2f} min"
    else:
        timing = f"{elapsed:.1f} s"
    _LOG.info(f"findGEVE completed in {timing}")
    log_run_summary(raw_geves, args.genome)
    return 0

if __name__ == "__main__":
    sys.exit(main())
