#!/usr/bin/env python3
"""
findGEVE_review.py - Apply manual review decisions to findGEVE outputs.
Author: Dede Kurniawan (dedekurniawan@genomics.cn)
"""

from __future__ import annotations

import argparse
import gzip
import logging
import math
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
import pandas as pd

try:
    from openpyxl import Workbook, load_workbook
    from openpyxl.styles import Alignment, Border, Font, PatternFill, Protection, Side
    from openpyxl.worksheet.datavalidation import DataValidation
except ImportError as exc:
    raise SystemExit(
        "Error: openpyxl is required for review.xlsx support. "
        "Install with: pip install openpyxl"
    ) from exc

HELP_TEXT = """\
findGEVE_review.py - Manual review helper for findGEVE results.

Usage:
  findGEVE_review.py make-template <prefix.summary.tsv>
  findGEVE_review.py apply \
    --review <prefix.review.xlsx> \
    --summary <prefix.summary.tsv> \
    --markerout <prefix.markerout> \
    --bed <prefix.geve.bed> \
    --genome genome.fa [OPTIONS]

Review actions:
  unchanged   Keep the original GEVE call.
  remove      Remove the GEVE from reviewed outputs.
  change      Use review_start and review_end as a curated candidate interval.
"""

_LOG = logging.getLogger("findGEVE_review")
OUTPUT = 25
logging.addLevelName(OUTPUT, "OUTPUT")

def _output(self, message, *args, **kwargs):
    if self.isEnabledFor(OUTPUT):
        self._log(OUTPUT, message, args, **kwargs)

logging.Logger.output = _output

ACTIONS = ("unchanged", "remove", "change")
REVIEW_COLUMNS = [
    "geve_name",
    "action",
    "contig",
    "original_start",
    "original_end",
    "review_start",
    "review_end",
]

FEATURE_IGNORE = {
    "GEVE",
    "flank_left",
    "flank_right",
    "TIR_left",
    "TIR_right",
    "TSD_5p",
    "TSD_3p",
}

BLASTN_COLUMNS = [
    "qstart", "qend", "sstart", "send", "length", "nident",
    "pident", "gaps", "evalue", "bitscore",
]
BLASTN_OUTFMT = "6 " + " ".join(BLASTN_COLUMNS)
NATKEY_SPLIT = re.compile(r"(\d+)")
REVCOMP_TABLE = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

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

def setup_logging(log_path: Optional[Path] = None) -> None:
    _LOG.setLevel(logging.DEBUG)
    _LOG.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(fmt)
    sh.setLevel(logging.INFO)
    _LOG.addHandler(sh)
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_path, mode="w")
        fh.setFormatter(fmt)
        fh.setLevel(logging.DEBUG)
        _LOG.addHandler(fh)

def _natural_key(s: str):
    return [int(t) if t.isdigit() else t.lower() for t in NATKEY_SPLIT.split(str(s))]

def _require_columns(df: pd.DataFrame, cols: Iterable[str], label: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise SystemExit(f"Error: {label} is missing required column(s): {', '.join(missing)}")

def _read_table(path: Path, label: str) -> pd.DataFrame:
    if path is None or not path.is_file():
        raise SystemExit(f"Error: {label} file not found: {path}")
    try:
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except Exception as exc:
        raise SystemExit(f"Error: failed to read {label}: {path}: {exc}") from exc

def _safe_int(value, default: Optional[int] = None) -> Optional[int]:
    if value is None:
        return default
    if isinstance(value, float) and math.isnan(value):
        return default
    s = str(value).strip()
    if not s or s.upper() in {"NA", "NAN", "NONE"}:
        return default
    try:
        return int(float(s.replace(",", "")))
    except ValueError:
        return default

def _safe_float(value, default: float = float("nan")) -> float:
    if value is None:
        return default
    s = str(value).strip()
    if not s or s.upper() in {"NA", "NAN", "NONE"}:
        return default
    try:
        return float(s.replace(",", ""))
    except ValueError:
        return default

def _normalize_action(value) -> str:
    return str(value or "").strip().lower()

def infer_prefix(summary_path: Path, summary: Optional[pd.DataFrame] = None) -> str:
    name = summary_path.name
    for suffix in (".summary.tsv", ".tsv", ".txt"):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break
    if name:
        return name
    if summary is not None and "geve_name" in summary.columns and not summary.empty:
        first = str(summary["geve_name"].iloc[0])
        return first.split("_GEVE_", 1)[0] if "_GEVE_" in first else first.split("_", 1)[0]
    return "findGEVE"

def default_outdir(base: Optional[Path]) -> Path:
    root = base if base is not None else Path.cwd()
    date_tag = datetime.now().strftime("%Y-%m-%d")
    out = root / f"review_{date_tag}"
    if not out.exists():
        return out
    idx = 1
    while True:
        cand = root / f"review_{date_tag}_{idx:02d}"
        if not cand.exists():
            return cand
        idx += 1

def read_fasta(path: Optional[Path]) -> Dict[str, str]:
    if path is None:
        return {}
    if not path.is_file():
        raise SystemExit(f"Error: genome FASTA not found: {path}")
    opener = gzip.open if str(path).endswith(".gz") else open
    seqs: Dict[str, List[str]] = {}
    current: Optional[str] = None
    try:
        with opener(path, "rt") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    current = line[1:].split()[0]
                    seqs.setdefault(current, [])
                elif current is not None:
                    seqs[current].append(line.upper())
    except Exception as exc:
        raise SystemExit(f"Error: failed to read genome FASTA {path}: {exc}") from exc
    return {k: "".join(v) for k, v in seqs.items()}

def fetch_seq(seqs: Dict[str, str], contig: str, start: int, end: int) -> str:
    seq = seqs.get(str(contig))
    if seq is None:
        return ""
    start = max(1, int(start))
    end = min(len(seq), int(end))
    if end < start:
        return ""
    return seq[start - 1:end]

def gc_of_seq(seq: str) -> float:
    s = seq.upper()
    gc = s.count("G") + s.count("C")
    at = s.count("A") + s.count("T")
    return float(100.0 * gc / (gc + at)) if (gc + at) else float("nan")

def wrap_fasta(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))

def revcomp(seq: str) -> str:
    return seq.translate(REVCOMP_TABLE)[::-1]

def translate_cds(seq: str) -> str:
    seq = seq.upper().replace("U", "T")
    aa = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        if any(b not in "ACGT" for b in codon):
            aa.append("X")
        else:
            aa.append(CODON_TABLE.get(codon, "X"))
    return "".join(aa).rstrip("*")

def viz_flank_size(geve_length: int) -> int:
    geve_length = max(1, int(geve_length))
    return int(min(200_000, max(10_000, round(geve_length * 0.10))))

def make_template(summary_path: Path, overwrite: bool = False) -> Path:
    summary = _read_table(summary_path, "summary")
    _require_columns(summary, ["geve_name", "start", "end"], "summary")
    contig_col = "contig_id" if "contig_id" in summary.columns else "contig"
    if contig_col not in summary.columns:
        raise SystemExit("Error: summary must contain contig_id or contig column")
    prefix = infer_prefix(summary_path, summary)
    out_path = summary_path.with_name(f"{prefix}.review.xlsx")
    if out_path.exists() and not overwrite:
        raise SystemExit(f"Error: output exists, use --overwrite to replace: {out_path}")

    wb = Workbook()
    ws = wb.active
    ws.title = "review"

    header_fill = PatternFill("solid", fgColor="1F4E78")
    locked_fill = PatternFill("solid", fgColor="D9EAF7")
    editable_fill = PatternFill("solid", fgColor="FFF2CC")
    thin = Side(style="thin", color="B7B7B7")
    border = Border(left=thin, right=thin, top=thin, bottom=thin)

    ws.append(REVIEW_COLUMNS)
    for cell in ws[1]:
        cell.font = Font(bold=True, color="FFFFFF")
        cell.fill = header_fill
        cell.alignment = Alignment(horizontal="center")
        cell.border = border

    for _, row in summary.sort_values("geve_name", key=lambda c: c.map(_natural_key)).iterrows():
        ws.append([
            row.get("geve_name", ""),
            "unchanged",
            row.get(contig_col, ""),
            _safe_int(row.get("start", "")),
            _safe_int(row.get("end", "")),
            "",
            "",
        ])

    dv = DataValidation(type="list", formula1='"unchanged,remove,change"', allow_blank=False)
    dv.error = "Choose one of: unchanged, remove, change"
    dv.errorTitle = "Invalid action"
    dv.prompt = "Choose unchanged, remove, or change"
    dv.promptTitle = "GEVE review action"
    ws.add_data_validation(dv)
    if ws.max_row >= 2:
        dv.add(f"B2:B{ws.max_row}")

    for row in ws.iter_rows(min_row=2, max_row=ws.max_row):
        for idx, cell in enumerate(row, start=1):
            cell.border = border
            if idx in (2, 6, 7):
                cell.fill = editable_fill
                cell.protection = Protection(locked=False)
            else:
                cell.fill = locked_fill
                cell.protection = Protection(locked=True)

    widths = {"A": 28, "B": 14, "C": 28, "D": 16, "E": 16, "F": 16, "G": 16}
    for col, width in widths.items():
        ws.column_dimensions[col].width = width
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = f"A1:G{ws.max_row}"
    ws.protection.sheet = True
    ws.protection.enable()

    note = wb.create_sheet("README")
    lines = [
        "findGEVE review.xlsx instructions",
        "Allowed action values: unchanged, remove, change",
        "unchanged: keep the original GEVE call. Leave review_start/review_end empty.",
        "remove: drop the GEVE from reviewed outputs. Leave review_start/review_end empty.",
        "change: fill review_start and review_end. These coordinates become the curated candidate interval.",
        "Use the Plotly HTML coordinate-review plot to click and copy boundary coordinates.",
    ]
    for line in lines:
        note.append([line])
    note["A1"].font = Font(bold=True, size=14)
    note.column_dimensions["A"].width = 120

    wb.save(out_path)
    _LOG.info(f"Wrote review template: {out_path}")
    return out_path

def read_review_xlsx(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise SystemExit(f"Error: review file not found: {path}")
    try:
        wb = load_workbook(path, data_only=True)
    except Exception as exc:
        raise SystemExit(f"Error: failed to open review workbook {path}: {exc}") from exc
    if "review" not in wb.sheetnames:
        raise SystemExit("Error: review workbook must contain a sheet named 'review'")
    ws = wb["review"]
    rows = list(ws.iter_rows(values_only=True))
    if not rows:
        raise SystemExit("Error: review sheet is empty")
    header = [str(x).strip() if x is not None else "" for x in rows[0]]
    df = pd.DataFrame(rows[1:], columns=header).dropna(how="all")
    _require_columns(df, REVIEW_COLUMNS, "review.xlsx")
    return df[REVIEW_COLUMNS].copy()

def validate_review(review: pd.DataFrame, summary: pd.DataFrame) -> Tuple[pd.DataFrame, List[str], List[str]]:
    _require_columns(summary, ["geve_name", "start", "end"], "summary")
    contig_col = "contig_id" if "contig_id" in summary.columns else "contig"
    if contig_col not in summary.columns:
        raise SystemExit("Error: summary must contain contig_id or contig column")
    known = {str(r["geve_name"]): r for _, r in summary.iterrows()}

    errors: List[str] = []
    warnings: List[str] = []
    clean_rows = []
    seen = set()

    for i, row in review.iterrows():
        excel_row = i + 2
        geve_name = str(row.get("geve_name", "")).strip()
        action = _normalize_action(row.get("action", ""))
        if not geve_name:
            errors.append(f"row {excel_row}: geve_name is empty")
            continue
        if geve_name in seen:
            errors.append(f"row {excel_row}: duplicate geve_name: {geve_name}")
            continue
        seen.add(geve_name)
        if geve_name not in known:
            errors.append(f"row {excel_row}: geve_name not present in summary: {geve_name}")
            continue
        if action not in ACTIONS:
            errors.append(
                f"row {excel_row} ({geve_name}): invalid action {row.get('action')!r}; "
                "must be unchanged, remove, or change"
            )
            continue

        original = known[geve_name]
        expected_contig = str(original.get(contig_col, "")).strip()
        contig = str(row.get("contig", "")).strip()
        if contig and contig != expected_contig:
            errors.append(
                f"row {excel_row} ({geve_name}): contig {contig!r} does not match summary {expected_contig!r}"
            )

        orig_start = _safe_int(original.get("start"))
        orig_end = _safe_int(original.get("end"))
        if orig_start is None or orig_end is None:
            errors.append(f"row {excel_row} ({geve_name}): original summary start/end is not numeric")
            continue

        sheet_orig_start = _safe_int(row.get("original_start"))
        sheet_orig_end = _safe_int(row.get("original_end"))
        if sheet_orig_start is not None and sheet_orig_start != orig_start:
            warnings.append(f"row {excel_row} ({geve_name}): original_start differs from summary; summary value was used")
        if sheet_orig_end is not None and sheet_orig_end != orig_end:
            warnings.append(f"row {excel_row} ({geve_name}): original_end differs from summary; summary value was used")

        review_start = _safe_int(row.get("review_start"))
        review_end = _safe_int(row.get("review_end"))
        if action == "change":
            if review_start is None or review_end is None:
                errors.append(f"row {excel_row} ({geve_name}): change requires review_start and review_end")
            elif review_start >= review_end:
                errors.append(f"row {excel_row} ({geve_name}): review_start must be smaller than review_end")
        elif review_start is not None or review_end is not None:
            warnings.append(f"row {excel_row} ({geve_name}): review_start/review_end ignored for action={action}")

        clean_rows.append(dict(
            geve_name=geve_name,
            action=action,
            contig=expected_contig,
            original_start=orig_start,
            original_end=orig_end,
            review_start=review_start,
            review_end=review_end,
        ))

    missing = sorted(set(known) - seen, key=_natural_key)
    for geve_name in missing:
        warnings.append(f"{geve_name}: missing from review.xlsx; treated as unchanged")
        r = known[geve_name]
        clean_rows.append(dict(
            geve_name=geve_name,
            action="unchanged",
            contig=str(r.get(contig_col, "")),
            original_start=_safe_int(r.get("start")),
            original_end=_safe_int(r.get("end")),
            review_start=None,
            review_end=None,
        ))

    clean = pd.DataFrame(clean_rows)
    if not clean.empty:
        clean = clean.sort_values("geve_name", key=lambda c: c.map(_natural_key)).reset_index(drop=True)
    return clean, errors, warnings

def parse_blastn_tabular(tab_path: Path) -> List[TirPair]:
    if not tab_path.exists() or not tab_path.read_text().strip():
        return []
    pairs: List[TirPair] = []
    for line in tab_path.read_text().splitlines():
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 10:
            continue
        try:
            qstart, qend = int(fields[0]), int(fields[1])
            sstart, send = int(fields[2]), int(fields[3])
            aln_len, nident = int(fields[4]), int(fields[5])
            pident, gaps = float(fields[6]), int(fields[7])
            evalue, bitscore = float(fields[8]), float(fields[9])
        except ValueError:
            continue
        left_start, left_end = min(qstart, qend), max(qstart, qend)
        right_start, right_end = min(sstart, send), max(sstart, send)
        if left_end >= right_start:
            continue
        tir_length = left_end - left_start + 1
        insert_size = right_start - left_end - 1
        pairs.append(TirPair(left_start, left_end, right_start, right_end,
                             tir_length, insert_size, pident, int(round(bitscore)),
                             nident, aln_len, gaps, evalue))
    return pairs

def run_blastn_self(region_seq: str, cfg: dict, threads: int = 1) -> List[TirPair]:
    blastn = shutil.which("blastn")
    if blastn is None:
        raise RuntimeError("blastn executable was not found in PATH")
    with tempfile.TemporaryDirectory(prefix="findGEVE_review_tir_") as tmp:
        tmpdir = Path(tmp)
        fa_path = tmpdir / "candidate.fa"
        tab_path = tmpdir / "candidate.blastn.tsv"
        fa_path.write_text(">candidate\n" + wrap_fasta(region_seq) + "\n")
        cmd = [
            blastn, "-query", str(fa_path), "-subject", str(fa_path),
            "-strand", "minus", "-task", "blastn",
            "-word_size", str(cfg["blastn_word_size"]),
            "-reward", str(cfg["blastn_reward"]),
            "-penalty", str(cfg["blastn_penalty"]),
            "-gapopen", str(cfg["blastn_gapopen"]),
            "-gapextend", str(cfg["blastn_gapextend"]),
            "-evalue", str(cfg["blastn_evalue"]),
            "-dust", "no", "-soft_masking", "false",
            "-max_target_seqs", str(cfg["blastn_max_targets"]),
            "-num_threads", str(max(1, int(threads))),
            "-outfmt", BLASTN_OUTFMT,
            "-out", str(tab_path),
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(proc.stderr.strip() or "blastn failed")
        return parse_blastn_tabular(tab_path)

def _count_bracketed(tir: TirPair, intervals: List[Tuple[int, int]]) -> int:
    return sum(1 for s, e in intervals if tir.left_start <= s and e <= tir.right_end)

def select_best_tir(pairs: List[TirPair], region_offset: int, candidate_start: int,
                    candidate_end: int, hallmark_intervals: List[Tuple[int, int]], cfg: dict) -> Optional[TirPair]:
    valid = []
    for t in pairs:
        abs_t = TirPair(
            left_start=t.left_start + region_offset - 1,
            left_end=t.left_end + region_offset - 1,
            right_start=t.right_start + region_offset - 1,
            right_end=t.right_end + region_offset - 1,
            tir_length=t.tir_length,
            insert_size=t.insert_size,
            tir_identity=t.tir_identity,
            score=t.score,
            matches=t.matches,
            total=t.total,
            gaps=t.gaps,
            tir_evalue=t.tir_evalue,
        )
        if not (cfg["tir_min_len"] <= abs_t.tir_length <= cfg["tir_max_len"]):
            continue
        if abs_t.tir_identity < cfg["tir_min_id"]:
            continue
        if not (cfg["tir_min_insert"] <= abs_t.insert_size <= cfg["tir_max_insert"]):
            continue
        edge_distance = abs(abs_t.left_start - candidate_start) + abs(candidate_end - abs_t.right_end)
        if edge_distance > cfg["tir_edge_slop"]:
            continue
        bracketed = _count_bracketed(abs_t, hallmark_intervals)
        if hallmark_intervals and bracketed < max(1, math.ceil(len(hallmark_intervals) * 0.5)):
            continue
        valid.append((bracketed, -edge_distance, abs_t.tir_identity, abs_t.tir_length, abs_t.score, abs_t))
    if not valid:
        return None
    valid.sort(reverse=True)
    return valid[0][-1]

def find_tsd(left_flank: str, right_flank: str, k_min: int, k_max: int, max_slide: int) -> Optional[Tsd]:
    left = left_flank.upper()
    right = right_flank.upper()
    best: Optional[Tsd] = None
    for k in range(k_max, k_min - 1, -1):
        if k > len(left) or k > len(right):
            continue
        max_mm = 0 if k <= 5 else (1 if k <= 8 else 2)
        for sl in range(max_slide + 1):
            if k + sl > len(left):
                break
            lk = left[len(left) - k - sl: len(left) - sl] if sl > 0 else left[-k:]
            for sr in range(max_slide + 1):
                if k + sr > len(right):
                    break
                rk = right[sr:sr + k]
                if "N" in lk or "N" in rk:
                    continue
                mm = sum(1 for a, b in zip(lk, rk) if a != b)
                if mm <= max_mm:
                    cand = Tsd(lk, rk, k, mm, 100.0 * (k - mm) / k, sl, sr)
                    if best is None or (cand.length, cand.identity) > (best.length, best.identity):
                        best = cand
        if best is not None and best.length == k:
            return best
    return best


def hallmark_intervals(marker: pd.DataFrame, geve_name: str, start: int, end: int) -> List[Tuple[int, int]]:
    if marker.empty:
        return []
    q = marker[(marker["geve_name"].astype(str) == geve_name) & (marker["feature"].astype(str) == "hallmark")]
    out = []
    for _, r in q.iterrows():
        s, e = _safe_int(r.get("start")), _safe_int(r.get("end"))
        if s is not None and e is not None and s >= start and e <= end:
            out.append((s, e))
    return out

def redetect_tir_tsd(seqs: Dict[str, str], marker: pd.DataFrame, geve_name: str,
                     contig: str, start: int, end: int, cfg: dict,
                     threads: int) -> Tuple[Optional[TirPair], Optional[Tsd], List[str]]:
    messages: List[str] = []
    if not seqs:
        messages.append("genome FASTA not provided; TIR/TSD detection skipped")
        return None, None, messages
    contig_seq = seqs.get(str(contig), "")
    if not contig_seq:
        messages.append(f"sequence unavailable for {contig}; TIR/TSD detection skipped")
        return None, None, messages
    search_start = max(1, start - cfg["tir_flank"])
    search_end = min(len(contig_seq), end + cfg["tir_flank"])
    region_seq = fetch_seq(seqs, contig, search_start, search_end)
    if not region_seq:
        messages.append(f"sequence unavailable for {contig}:{search_start}-{search_end}; TIR/TSD detection skipped")
        return None, None, messages
    try:
        raw_pairs = run_blastn_self(region_seq, cfg, threads=threads)
    except Exception as exc:
        messages.append(f"TIR detection skipped/failed: {exc}")
        return None, None, messages
    hms = hallmark_intervals(marker, geve_name, start, end)
    tir = select_best_tir(raw_pairs, search_start, start, end, hms, cfg)
    if tir is None:
        messages.append(f"no TIR passed filters among {len(raw_pairs)} raw inverted-repeat pairs")
        return None, None, messages
    left_flank = fetch_seq(seqs, contig, max(1, tir.left_start - cfg["tsd_flank"]), tir.left_start - 1)
    right_flank = fetch_seq(seqs, contig, tir.right_end + 1, tir.right_end + cfg["tsd_flank"])
    tsd = find_tsd(left_flank, right_flank, cfg["tsd_min_len"], cfg["tsd_max_len"], cfg["tsd_max_slide"])
    messages.append(
        f"TIR detected: {tir.left_start}-{tir.left_end} / {tir.right_start}-{tir.right_end}, "
        f"identity={tir.tir_identity:.2f}%"
    )
    if tsd is None:
        messages.append("TSD not detected")
    else:
        messages.append(f"TSD detected: {tsd.sequence_left}|{tsd.sequence_right}, len={tsd.length}")
    return tir, tsd, messages

def get_original_tir_tsd(marker: pd.DataFrame, geve_name: str) -> Tuple[Optional[TirPair], Optional[Tsd]]:
    q = marker[marker["geve_name"].astype(str) == geve_name]
    left = q[q["feature"].astype(str) == "TIR_left"]
    right = q[q["feature"].astype(str) == "TIR_right"]
    tir = None
    if not left.empty and not right.empty:
        l, r = left.iloc[0], right.iloc[0]
        ls, le = _safe_int(l.get("start")), _safe_int(l.get("end"))
        rs, re_ = _safe_int(r.get("start")), _safe_int(r.get("end"))
        if None not in (ls, le, rs, re_):
            score = _safe_int(l.get("score"), 0) or 0
            length = le - ls + 1
            tir = TirPair(ls, le, rs, re_, length, rs - le - 1, float("nan"), score, 0, 0, 0)
    tsd = None
    t5 = q[q["feature"].astype(str) == "TSD_5p"]
    t3 = q[q["feature"].astype(str) == "TSD_3p"]
    if not t5.empty and not t3.empty:
        a, b = t5.iloc[0], t3.iloc[0]
        sleft = str(a.get("name", "") or "")
        sright = str(b.get("name", "") or "")
        length = len(sleft) if sleft and sleft != "." else len(sright)
        ident = _safe_float(a.get("score"), 100.0)
        mismatches = int(round(length * (100.0 - ident) / 100.0)) if length else 0
        tsd = Tsd(sleft, sright, length, mismatches, ident, 0, 0)
    return tir, tsd

def tir_fields(tir: Optional[TirPair]) -> dict:
    if tir is None:
        return dict(tir_length="NA", tir_score="NA", tir_identity_pct="NA", tir_gaps="NA")
    return dict(
        tir_length=tir.tir_length,
        tir_score=tir.score,
        tir_identity_pct="NA" if math.isnan(tir.tir_identity) else round(tir.tir_identity, 2),
        tir_gaps=tir.gaps,
    )

def tsd_fields(tsd: Optional[Tsd]) -> dict:
    if tsd is None:
        return dict(tsd_len="NA", tsd_left="NA", tsd_right="NA", tsd_mismatch="NA", tsd_conservation="NODETECT")
    return dict(
        tsd_len=tsd.length,
        tsd_left=tsd.sequence_left,
        tsd_right=tsd.sequence_right,
        tsd_mismatch=tsd.mismatches,
        tsd_conservation="PERFECT" if tsd.mismatches == 0 else "IMPERFECT",
    )

def build_reviewed_records(review: pd.DataFrame, summary: pd.DataFrame, marker: pd.DataFrame,
                           seqs: Dict[str, str], cfg: dict, threads: int) -> Tuple[List[dict], List[str]]:
    summary_by_name = {str(r["geve_name"]): r for _, r in summary.iterrows()}
    records: List[dict] = []
    messages: List[str] = []
    kept = review[review["action"] != "remove"].copy().reset_index(drop=True)
    for idx, row in kept.iterrows():
        old_name = row["geve_name"]
        new_name = re.sub(r"_GEVE_\d+$", f"_GEVE_{idx + 1:03d}", old_name)
        if new_name == old_name:
            prefix = old_name.split("_GEVE_", 1)[0] if "_GEVE_" in old_name else infer_prefix(Path("findGEVE.summary.tsv"))
            new_name = f"{prefix}_GEVE_{idx + 1:03d}"
        original = summary_by_name[old_name]
        contig = row["contig"]
        action = row["action"]
        if action == "change":
            candidate_start = int(row["review_start"])
            candidate_end = int(row["review_end"])
            tir, tsd, tir_messages = redetect_tir_tsd(seqs, marker, old_name, contig, candidate_start, candidate_end, cfg, threads)
            messages.extend([f"{old_name}: {m}" for m in tir_messages])
            final_start = tir.left_start if tir is not None else candidate_start
            final_end = tir.right_end if tir is not None else candidate_end
            boundary_method = "reviewed_tir_boundary" if tir is not None else "reviewed_manual_boundary"
        else:
            candidate_start = int(row["original_start"])
            candidate_end = int(row["original_end"])
            final_start = candidate_start
            final_end = candidate_end
            tir, tsd = get_original_tir_tsd(marker, old_name)
            boundary_method = "original_boundary"
        geve_len = final_end - final_start + 1
        seq = fetch_seq(seqs, contig, final_start, final_end)
        gc = gc_of_seq(seq) if seq else _safe_float(original.get("gc"))
        records.append(dict(
            original_geve_name=old_name,
            reviewed_geve_name=new_name,
            action=action,
            contig=contig,
            candidate_start=candidate_start,
            candidate_end=candidate_end,
            geve_start=final_start,
            geve_end=final_end,
            geve_length=geve_len,
            gc_geve=gc,
            tir=tir,
            tsd=tsd,
            has_tir=tir is not None,
            boundary_method=boundary_method,
            original_summary=original.to_dict(),
        ))
    return records, messages

def write_reviewed_summary(records: List[dict], path: Path) -> None:
    rows = []
    for r in records:
        old = dict(r["original_summary"])
        row = dict(old)
        row["original_geve_name"] = r["original_geve_name"]
        row["geve_name"] = r["reviewed_geve_name"]
        row["review_action"] = r["action"]
        row["boundary_method"] = r["boundary_method"]
        row["contig_id"] = r["contig"]
        row["start"] = r["geve_start"]
        row["end"] = r["geve_end"]
        row["geve_length"] = r["geve_length"]
        row["gc"] = "NA" if math.isnan(r["gc_geve"]) else round(r["gc_geve"], 2)
        row["has_tir"] = "yes" if r["has_tir"] else "no"
        row.update(tir_fields(r["tir"]))
        row.update(tsd_fields(r["tsd"]))
        rows.append(row)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)

def write_reviewed_fna(records: List[dict], seqs: Dict[str, str], path: Path) -> None:
    with path.open("w") as fh:
        for r in records:
            seq = fetch_seq(seqs, r["contig"], r["geve_start"], r["geve_end"])
            header = (
                f">{r['reviewed_geve_name']} contig={r['contig']} start={r['geve_start']} "
                f"end={r['geve_end']} length={r['geve_length']} boundary_method={r['boundary_method']}"
            )
            fh.write(header + "\n")
            fh.write(wrap_fasta(seq) + "\n")

def feature_rows_for_record(marker: pd.DataFrame, record: dict, include_flank: bool = False) -> pd.DataFrame:
    old = record["original_geve_name"]
    start = record["geve_start"]
    end = record["geve_end"]
    if include_flank:
        flank = viz_flank_size(record["geve_length"])
        start = max(1, start - flank)
        end = end + flank
    q = marker[(marker["geve_name"].astype(str) == old) & (~marker["feature"].astype(str).isin(FEATURE_IGNORE))].copy()
    if q.empty:
        return q
    q["start_i"] = q["start"].map(_safe_int)
    q["end_i"] = q["end"].map(_safe_int)
    q = q.dropna(subset=["start_i", "end_i"])
    q["start_i"] = q["start_i"].astype(int)
    q["end_i"] = q["end_i"].astype(int)
    q = q[(q["end_i"] >= start) & (q["start_i"] <= end)].copy()
    return q.sort_values(["start_i", "end_i"], kind="mergesort")

def write_reviewed_cds_pep(records: List[dict], marker: pd.DataFrame, seqs: Dict[str, str], cds_path: Path,
                           pep_path: Path, hallmark_dir: Path, prefix: str) -> None:
    hallmark_best: Dict[str, Dict[str, Tuple[str, str]]] = {}
    with cds_path.open("w") as cds_fh, pep_path.open("w") as pep_fh:
        for r in records:
            gid = r["reviewed_geve_name"]
            feats = feature_rows_for_record(marker, r, include_flank=False)
            for orf_idx, (_, row) in enumerate(feats.iterrows(), start=1):
                s, e = int(row["start_i"]), int(row["end_i"])
                strand = str(row.get("strand", "+") or "+")
                cds = fetch_seq(seqs, r["contig"], s, e)
                if strand == "-":
                    cds = revcomp(cds)
                pep = translate_cds(cds)
                feature = str(row.get("feature", "orf") or "orf")
                name = str(row.get("name", ".") or ".")
                label = f"orf{orf_idx:05d}"
                annot = name if feature == "hallmark" and name not in {"", "."} else ""
                cds_header = f">{gid}_{label}"
                pep_header = f">{gid}_{label}"
                if annot:
                    cds_header += f" {annot}"
                    pep_header += f" {annot}"
                cds_header += f" contig={r['contig']} start={s} end={e} strand={strand} length={len(cds)}"
                pep_header += f" length={len(pep)}"
                cds_fh.write(cds_header + "\n" + wrap_fasta(cds) + "\n")
                pep_fh.write(pep_header + "\n" + wrap_fasta(pep) + "\n")
                if feature == "hallmark" and annot:
                    current = hallmark_best.setdefault(annot, {}).get(gid)
                    if current is None or len(pep) > len(current[1]):
                        hallmark_best[annot][gid] = (f">{gid}_{annot}", pep)
    hallmark_dir.mkdir(parents=True, exist_ok=True)
    for hallmark, geve_map in sorted(hallmark_best.items(), key=lambda kv: _natural_key(kv[0])):
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", hallmark.lower())
        out = hallmark_dir / f"{prefix}.{safe}.pep"
        with out.open("w") as fh:
            for gid, (header, pep) in sorted(geve_map.items(), key=lambda kv: _natural_key(kv[0])):
                fh.write(header + "\n" + wrap_fasta(pep) + "\n")

def write_reviewed_gff3(records: List[dict], marker: pd.DataFrame, path: Path) -> None:
    with path.open("w") as fh:
        fh.write("##gff-version 3\n")
        for r in records:
            gid = r["reviewed_geve_name"]
            contig = r["contig"]
            fh.write(
                f"{contig}\tfindGEVE_review\tmobile_genetic_element\t{r['geve_start']}\t{r['geve_end']}\t.\t+\t.\t"
                f"ID={gid};Name={gid};original_geve_name={r['original_geve_name']};boundary_method={r['boundary_method']}\n"
            )
            feats = feature_rows_for_record(marker, r, include_flank=False)
            for orf_idx, (_, row) in enumerate(feats.iterrows(), start=1):
                feature = str(row.get("feature", "ORF") or "ORF")
                name = str(row.get("name", ".") or ".")
                strand = str(row.get("strand", ".") or ".")
                score = str(row.get("score", ".") or ".")
                fh.write(
                    f"{contig}\tfindGEVE_review\t{feature}\t{int(row['start_i'])}\t{int(row['end_i'])}\t{score}\t{strand}\t.\t"
                    f"ID={gid}.orf{orf_idx:05d};Parent={gid};Name={name}\n"
                )

def write_reviewed_markerout(records: List[dict], marker: pd.DataFrame, seqs: Dict[str, str], path: Path) -> None:
    with path.open("w") as fh:
        fh.write("contig\tgeve_name\tfeature\tname\tstart\tend\tstrand\te_value\tscore\n")
        for r in records:
            gid = r["reviewed_geve_name"]
            old = r["original_geve_name"]
            contig = r["contig"]
            gstart, gend = r["geve_start"], r["geve_end"]
            clen = len(seqs.get(contig, "")) if seqs else gend
            flank = viz_flank_size(r["geve_length"])
            region_start = max(1, gstart - flank)
            region_end = min(clen, gend + flank) if clen else gend + flank
            fh.write(f"{contig}\t{gid}\tGEVE\t.\t{gstart}\t{gend}\t.\tNA\t{r['geve_length']}\n")
            if region_start < gstart:
                fh.write(f"{contig}\t{gid}\tflank_left\t.\t{region_start}\t{gstart - 1}\t.\tNA\tNA\n")
            if region_end > gend:
                fh.write(f"{contig}\t{gid}\tflank_right\t.\t{gend + 1}\t{region_end}\t.\tNA\tNA\n")
            tir = r["tir"]
            if tir is not None:
                fh.write(f"{contig}\t{gid}\tTIR_left\t.\t{tir.left_start}\t{tir.left_end}\t+\tNA\t{tir.score}\n")
                fh.write(f"{contig}\t{gid}\tTIR_right\t.\t{tir.right_start}\t{tir.right_end}\t-\tNA\t{tir.score}\n")
            tsd = r["tsd"]
            if tsd is not None and tir is not None:
                ltsd_end = tir.left_start - 1 - tsd.left_shift
                ltsd_start = ltsd_end - tsd.length + 1
                rtsd_start = tir.right_end + 1 + tsd.right_shift
                rtsd_end = rtsd_start + tsd.length - 1
                fh.write(f"{contig}\t{gid}\tTSD_5p\t{tsd.sequence_left}\t{ltsd_start}\t{ltsd_end}\t+\tNA\t{tsd.identity:.1f}\n")
                fh.write(f"{contig}\t{gid}\tTSD_3p\t{tsd.sequence_right}\t{rtsd_start}\t{rtsd_end}\t+\tNA\t{tsd.identity:.1f}\n")
            q = marker[(marker["geve_name"].astype(str) == old) & (~marker["feature"].astype(str).isin(FEATURE_IGNORE))].copy()
            for _, row in q.iterrows():
                s, e = _safe_int(row.get("start")), _safe_int(row.get("end"))
                if s is None or e is None or e < region_start or s > region_end:
                    continue
                vals = [
                    contig, gid, row.get("feature", "."), row.get("name", "."),
                    str(s), str(e), row.get("strand", "."), row.get("e_value", "NA"), row.get("score", "NA"),
                ]
                fh.write("\t".join(map(str, vals)) + "\n")

def write_reviewed_bed(records: List[dict], bed: pd.DataFrame, path: Path) -> None:
    if bed.empty:
        path.write_text("contig_id\twindow_start\twindow_end\tgeve_name\trel_start\trel_end\tregion_type\tgc\trolling_score_mean\tn_orfs\tgvog_hits\tpfam_hits\n")
        return
    out_rows = []
    for r in records:
        gid = r["reviewed_geve_name"]
        old = r["original_geve_name"]
        gstart, gend = r["geve_start"], r["geve_end"]
        flank = viz_flank_size(r["geve_length"])
        region_start = max(1, gstart - flank)
        region_end = gend + flank
        q = bed[bed["geve_name"].astype(str) == old].copy()
        if q.empty:
            continue
        for col in ["window_start", "window_end"]:
            q[col] = pd.to_numeric(q[col], errors="coerce")
        q = q[(q["window_end"] >= region_start) & (q["window_start"] <= region_end)].copy()
        if q.empty:
            continue
        centers = ((q["window_start"] + q["window_end"]) / 2.0)
        q["geve_name"] = gid
        q["region_type"] = ["flank_left" if c < gstart else ("geve" if c <= gend else "flank_right") for c in centers]
        out_rows.append(q)
    if out_rows:
        pd.concat(out_rows, ignore_index=True).to_csv(path, sep="\t", index=False)
    else:
        path.write_text("contig_id\twindow_start\twindow_end\tgeve_name\trel_start\trel_end\tregion_type\tgc\trolling_score_mean\tn_orfs\tgvog_hits\tpfam_hits\n")


def run_geve_plot(marker_path: Path, bed_path: Path, outdir: Path) -> None:
    candidates = [
        Path(__file__).resolve().parent / "findGEVE_plot_v4.py",
        Path(__file__).resolve().parent / "findGEVE_plot.py",
    ]
    script = next((p for p in candidates if p.is_file()), None)
    if script is None:
        _LOG.warning("findGEVE_plot.py not found; skipping plot step")
        return
    if not bed_path.is_file() or not bed_path.read_text().strip():
        _LOG.warning("reviewed geve.bed is empty; skipping plot step")
        return
    cmd = [sys.executable, str(script), str(marker_path.resolve()), str(bed_path.resolve())]
    _LOG.info(f"Plotting GEVEs: {' '.join(cmd)}")
    try:
        proc = subprocess.run(cmd, cwd=str(outdir), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
    except OSError as exc:
        _LOG.warning(f"Failed to launch findGEVE_plot.py: {exc}")
        return
    for line in proc.stdout.splitlines():
        if line.strip():
            _LOG.info(f"[PLOT] {line}")
    for line in proc.stderr.splitlines():
        if line.strip():
            _LOG.warning(f"[PLOT] {line}")
    if proc.returncode != 0:
        _LOG.warning(f"findGEVE_plot.py exited with code {proc.returncode}; no plot produced")

def apply_review(args) -> None:
    summary = _read_table(args.summary, "summary")
    marker = _read_table(args.markerout, "markerout")
    bed = _read_table(args.bed, "geve.bed")
    review_raw = read_review_xlsx(args.review)
    review, errors, warnings = validate_review(review_raw, summary)
    prefix = args.prefix or infer_prefix(args.summary, summary)
    outdir = args.outdir or default_outdir(args.outbase)
    if outdir.exists() and any(outdir.iterdir()) and not args.overwrite:
        raise SystemExit(f"Error: output directory is not empty; use --overwrite: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    setup_logging(outdir / "review.log")

    _LOG.info("findGEVE review started")
    _LOG.info(f"Command line | {' '.join(sys.argv)}")
    _LOG.info(f"Review file | {args.review}")
    _LOG.info(f"Summary file | {args.summary}")
    _LOG.info(f"Markerout file | {args.markerout}")
    _LOG.info(f"BED file | {args.bed}")
    _LOG.info(f"Genome file | {args.genome}")
    _LOG.info(f"Output directory | {outdir}")

    if errors:
        for err in errors:
            _LOG.error(err)
        raise SystemExit(f"Error: review validation failed with {len(errors)} error(s); see {outdir / 'review.log'}")
    for warn in warnings:
        _LOG.warning(warn)

    seqs = read_fasta(args.genome) if args.genome is not None else {}
    if args.genome is not None:
        _LOG.info(f"Loaded genome FASTA: {len(seqs):,} contig(s)")
    else:
        _LOG.warning("Genome FASTA not provided; FASTA/CDS/PEP output will be empty and changed-GEVE TIR search will be skipped")

    cfg = dict(
        tir_flank=args.tir_flank,
        tir_min_len=args.tir_min_len,
        tir_max_len=args.tir_max_len,
        tir_min_id=args.tir_min_id,
        tir_min_insert=args.tir_min_insert,
        tir_max_insert=args.tir_max_insert,
        tir_edge_slop=args.tir_edge_slop,
        blastn_word_size=args.blastn_word_size,
        blastn_reward=args.blastn_reward,
        blastn_penalty=args.blastn_penalty,
        blastn_gapopen=args.blastn_gapopen,
        blastn_gapextend=args.blastn_gapextend,
        blastn_evalue=args.blastn_evalue,
        blastn_max_targets=args.blastn_max_targets,
        tsd_min_len=args.tsd_min_len,
        tsd_max_len=args.tsd_max_len,
        tsd_max_slide=args.tsd_max_slide,
        tsd_flank=args.tsd_flank,
    )
    _LOG.info(f"Review TIR flank | {args.tir_flank:,} bp")

    n_remove = int((review["action"] == "remove").sum())
    n_change = int((review["action"] == "change").sum())
    _LOG.info(f"Review rows | total={len(review):,} unchanged={int((review['action'] == 'unchanged').sum()):,} change={n_change:,} remove={n_remove:,}")

    records, tir_messages = build_reviewed_records(review, summary, marker, seqs, cfg, args.threads)
    for msg in tir_messages:
        _LOG.info(msg)

    summary_path = outdir / f"{prefix}.reviewed.summary.tsv"
    marker_path = outdir / f"{prefix}.reviewed.markerout"
    bed_path = outdir / f"{prefix}.reviewed.geve.bed"
    fna_path = outdir / f"{prefix}.reviewed.geve.fna"
    gff_path = outdir / f"{prefix}.reviewed.geve.gff3"
    cds_path = outdir / f"{prefix}.reviewed.geve.cds"
    pep_path = outdir / f"{prefix}.reviewed.geve.pep"
    hallmark_dir = outdir / "hallmark"

    write_reviewed_summary(records, summary_path)
    _LOG.output(f"Wrote {summary_path}")
    write_reviewed_markerout(records, marker, seqs, marker_path)
    _LOG.output(f"Wrote {marker_path}")
    write_reviewed_bed(records, bed, bed_path)
    _LOG.output(f"Wrote {bed_path}")
    write_reviewed_fna(records, seqs, fna_path)
    _LOG.output(f"Wrote {fna_path}")
    write_reviewed_gff3(records, marker, gff_path)
    _LOG.output(f"Wrote {gff_path}")
    write_reviewed_cds_pep(records, marker, seqs, cds_path, pep_path, hallmark_dir, prefix)
    _LOG.output(f"Wrote {cds_path}")
    _LOG.output(f"Wrote {pep_path}")
    _LOG.output(f"Wrote hallmark protein folder: {hallmark_dir}")

    if not args.no_plot:
        run_geve_plot(marker_path, bed_path, outdir)

    _LOG.info("Result Summary")
    _LOG.info(f"GEVE reviewed: {len(review):,}")
    _LOG.info(f"  Retained: {len(records):,}")
    _LOG.info(f"  Removed:  {n_remove:,}")
    _LOG.info(f"  Changed:  {n_change:,}")
    n_tir = sum(1 for r in records if r["tir"] is not None)
    _LOG.info(f"  With TIR: {n_tir:,}")
    _LOG.info(f"Review log: {outdir / 'review.log'}")

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Manual review helper for findGEVE results.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=HELP_TEXT,
    )
    sub = parser.add_subparsers(dest="command", required=True)
    p_template = sub.add_parser("make-template", help="Create <prefix>.review.xlsx")
    p_template.add_argument("summary", type=Path, help="Input <prefix>.summary.tsv")
    p_template.add_argument("--overwrite", action="store_true", help="Overwrite existing output")
    p_apply = sub.add_parser("apply", help="Apply reviewed Excel file")
    p_apply.add_argument("--review", required=True, type=Path, help="Reviewed Excel file")
    p_apply.add_argument("--summary", required=True, type=Path, help="Original <prefix>.summary.tsv")
    p_apply.add_argument("--markerout", required=True, type=Path, help="Original <prefix>.markerout")
    p_apply.add_argument("--bed", required=True, type=Path, help="Original <prefix>.geve.bed")
    p_apply.add_argument("--genome", type=Path, help="Genome FASTA; gzip is acceptable")
    p_apply.add_argument("--prefix", help="Output prefix; inferred from summary when omitted")
    p_apply.add_argument("--outdir", type=Path, help="Output review directory")
    p_apply.add_argument("--outbase", type=Path, help="Base directory for automatic review_<date> output folder")
    p_apply.add_argument("--overwrite", action="store_true", help="Allow writing into a non-empty output directory")
    p_apply.add_argument("--no-plot", action="store_true", help="Skip automatic plotting")
    p_apply.add_argument("-t", "--threads", type=int, default=1, help="blastn threads [default: 1]")
    p_apply.add_argument("--tir-flank", type=int, default=50_000, help="TIR search flank around changed interval [default: 50000]")
    p_apply.add_argument("--tir-min-len", type=int, default=20)
    p_apply.add_argument("--tir-max-len", type=int, default=5000)
    p_apply.add_argument("--tir-min-id", type=float, default=75.0)
    p_apply.add_argument("--tir-min-insert", type=int, default=1000)
    p_apply.add_argument("--tir-max-insert", type=int, default=2_000_000)
    p_apply.add_argument("--tir-edge-slop", type=int, default=50_000)
    p_apply.add_argument("--tsd-min-len", type=int, default=4)
    p_apply.add_argument("--tsd-max-len", type=int, default=12)
    p_apply.add_argument("--tsd-max-slide", type=int, default=3)
    p_apply.add_argument("--tsd-flank", type=int, default=100)
    p_apply.add_argument("--blastn-word-size", type=int, default=7)
    p_apply.add_argument("--blastn-reward", type=int, default=1)
    p_apply.add_argument("--blastn-penalty", type=int, default=-1)
    p_apply.add_argument("--blastn-gapopen", type=int, default=2)
    p_apply.add_argument("--blastn-gapextend", type=int, default=1)
    p_apply.add_argument("--blastn-evalue", type=float, default=10.0)
    p_apply.add_argument("--blastn-max-targets", type=int, default=10000)
    return parser

def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "make-template":
        setup_logging()
        make_template(args.summary, overwrite=args.overwrite)
        return 0
    if args.command == "apply":
        setup_logging()
        apply_review(args)
        return 0
    parser.error("unknown command")
    return 2

if __name__ == "__main__":
    sys.exit(main())
