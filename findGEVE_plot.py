#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib import patches
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.utils import PlotlyJSONEncoder

USAGE = """Usage: python findGEVE_plot.py [options] <prefix.markerout> <prefix.geve.bed>

Options:
  -g, --gff FILE       Gene annotation GFF3
  -r, --repeat FILE    Repeat/TE annotation GFF3

Output files:
  <prefix>.plot.pdf   Matplotlib static report
  <prefix>.plot.html  Plotly interactive evidence report
"""

PLOT_COLORS = {
    "viral_score": "#ec3028",
    "gc_content": "#fec44a",
    "gc_midline": "#d9d9d9",
    "gvog": "#1497a5",
    "pfam": "#756bb1",
    "intron": "#dd53ea", 
    "exon": "#3ba9ef",   
    "repeat": "#ffa46c",   
    "geve_highlight": "#e6e6e6",
    "track_background": "#ffffff",
    "line_track_background": "#fbfbfb",
    "hallmark_default": "#636363",
    "hallmark_edge": "black",
}

PLOT_STYLE = {
    "hallmark_star_size": 13,
    "hallmark_star_edge_width": 0.0,
    "geve_highlight_alpha": 0.55,
    "title_y": 0.955,
    "legend_y": 0.045,
}

FIG_WIDTH = 14.0
FIG_HEIGHT = 8.5

_NATKEY_SPLIT = re.compile(r"(\d+)")

def _natural_key(s: str):
    return [int(t) if t.isdigit() else t.lower() for t in _NATKEY_SPLIT.split(str(s))]

def _read_markerout(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype={
        "contig": str, "geve_name": str, "feature": str,
        "name": str, "strand": str,
    })
    if df.empty:
        return df
    for col in ["start", "end"]:
        df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    return df

def _read_bed(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype={
        "contig_id": str, "geve_name": str, "region_type": str,
    })
    if df.empty:
        return df
    for col in ["window_start", "window_end", "rel_start", "rel_end",
                "n_orfs", "gvog_hits", "pfam_hits"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in ["gc", "rolling_score_mean"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df

def _read_gff(path: Optional[Path]) -> pd.DataFrame:
    """Read a GFF/GFF3 file into a minimal interval dataframe.
    """
    if path is None:
        return pd.DataFrame(columns=["seqid", "source", "feature", "start", "end",
                                     "score", "strand", "phase", "attributes"])
    rows = []
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t", 8)
            if len(parts) < 9:
                continue
            rows.append(parts[:9])
    if not rows:
        return pd.DataFrame(columns=["seqid", "source", "feature", "start", "end",
                                     "score", "strand", "phase", "attributes"])
    df = pd.DataFrame(rows, columns=["seqid", "source", "feature", "start", "end",
                                     "score", "strand", "phase", "attributes"])
    df["seqid"] = df["seqid"].astype(str)
    df["feature"] = df["feature"].astype(str)
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["start", "end"]).copy()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    swap = df["start"] > df["end"]
    if swap.any():
        df.loc[swap, ["start", "end"]] = df.loc[swap, ["end", "start"]].to_numpy()
    return df

def _parse_gff_attributes(attributes: str) -> Dict[str, str]:
    """Parse GFF3 attributes into a dictionary.
    """
    out: Dict[str, str] = {}
    for item in str(attributes or "").strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(None, 1)
        else:
            key, value = item, ""
        out[key.strip()] = value.strip()
    return out

def _gff_attr(attributes: str, key: str) -> str:
    return _parse_gff_attributes(attributes).get(key, "")

def _derive_introns_from_exons(gff: pd.DataFrame) -> pd.DataFrame:
    """Infer intron records from exon gaps in GFF3 transcript models.
    """
    cols = ["seqid", "source", "feature", "start", "end", "score",
            "strand", "phase", "attributes"]
    if gff is None or gff.empty or "feature" not in gff:
        return pd.DataFrame(columns=cols)

    exons = gff[gff["feature"].astype(str).str.lower().eq("exon")].copy()
    if exons.empty:
        return pd.DataFrame(columns=cols)

    exons["_parent"] = exons["attributes"].map(lambda x: _gff_attr(x, "Parent"))
    exons["_id"] = exons["attributes"].map(lambda x: _gff_attr(x, "ID"))
    exons["_parent"] = exons["_parent"].astype(str).str.split(",").str[0]
    exons.loc[exons["_parent"].eq("") | exons["_parent"].isna(), "_parent"] = exons["_id"]

    intron_rows = []
    group_cols = ["seqid", "strand", "_parent"]
    for (seqid, strand, parent), grp in exons.groupby(group_cols, dropna=False, sort=False):
        if not parent or len(grp) < 2:
            continue
        intervals = sorted((int(s), int(e)) for s, e in grp[["start", "end"]].itertuples(index=False))
        merged: List[List[int]] = []
        for s, e in intervals:
            if not merged or s > merged[-1][1] + 1:
                merged.append([s, e])
            else:
                merged[-1][1] = max(merged[-1][1], e)
        for idx, (left, right) in enumerate(zip(merged, merged[1:]), start=1):
            intron_start = left[1] + 1
            intron_end = right[0] - 1
            if intron_start <= intron_end:
                intron_rows.append({
                    "seqid": str(seqid),
                    "source": "derived_from_exon",
                    "feature": "intron",
                    "start": intron_start,
                    "end": intron_end,
                    "score": ".",
                    "strand": str(strand),
                    "phase": ".",
                    "attributes": f"ID={parent}.derived_intron{idx};Parent={parent};Note=derived_from_exon_gaps",
                })

    if not intron_rows:
        return pd.DataFrame(columns=cols)
    return pd.DataFrame(intron_rows, columns=cols)

def _prepare_gene_gff(gff: Optional[pd.DataFrame]) -> Optional[pd.DataFrame]:
    """Normalize gene GFFs for plotting.
    """
    if gff is None or gff.empty:
        return gff
    if _has_feature(gff, "intron"):
        return gff
    derived = _derive_introns_from_exons(gff)
    if derived.empty:
        return gff
    return pd.concat([gff, derived], ignore_index=True)

def _window_overlap_fraction(
    windows: pd.DataFrame,
    intervals: pd.DataFrame,
    contig: str,
    feature_names: Optional[set[str]] = None,
) -> np.ndarray:
    """Return the fraction of each window covered by selected intervals.
    """
    out = np.zeros(len(windows), dtype=float)
    if intervals is None or intervals.empty or windows.empty:
        return out

    q = intervals[intervals["seqid"].astype(str) == str(contig)].copy()
    if feature_names is not None:
        q = q[q["feature"].str.lower().isin(feature_names)]
    if q.empty:
        return out

    w_start = windows["window_start"].to_numpy(dtype=float)
    w_end = windows["window_end"].to_numpy(dtype=float)
    x0 = np.nanmin(w_start)
    x1 = np.nanmax(w_end)

    q = q[(q["end"] >= x0) & (q["start"] <= x1)]
    if q.empty:
        return out

    interval_pairs = [(float(s), float(e)) for s, e in q[["start", "end"]].itertuples(index=False)]
    for i, (ws, we) in enumerate(zip(w_start, w_end)):
        overlaps = []
        for s, e in interval_pairs:
            os = max(ws, s)
            oe = min(we, e)
            if oe >= os:
                overlaps.append((os, oe))
        if not overlaps:
            continue
        overlaps.sort()
        merged = []
        for os, oe in overlaps:
            if not merged or os > merged[-1][1] + 1:
                merged.append([os, oe])
            else:
                merged[-1][1] = max(merged[-1][1], oe)
        covered = sum(oe - os + 1 for os, oe in merged)
        out[i] = covered / max(1.0, we - ws + 1)
    return out

def _has_feature(gff: Optional[pd.DataFrame], feature_name: str) -> bool:
    """True if a GFF dataframe contains at least one requested feature type."""
    if gff is None or gff.empty or "feature" not in gff:
        return False
    return gff["feature"].astype(str).str.lower().eq(feature_name.lower()).any()

_TE_SUPERFAMILIES = frozenset({
    "ltr", "line", "sine", "dna", "rc", "helitron", "mite",
    "dirs", "penelope", "maverick", "polinton", "crypton",
    "retroposon", "retrotransposon", "transposon",
})

_NON_TE_CLASSES = frozenset({
    "simple_repeat", "low_complexity", "satellite",
    "microsatellite", "minisatellite",
    "trna", "rrna", "snrna", "scrna", "srna",
    "rna", "buffer", "tandem_repeat", "unknown", 
})

def _repeat_target_class(attributes: str) -> str:
    """Extract the repeat class from RepeatMasker-style GFF attributes.
    """
    text = str(attributes or "")
    m = re.search(r"(?:^|;)\s*Target=([^;]+)", text)
    if not m:
        return ""
    target = m.group(1).strip()
    parts = target.split()
    if len(parts) >= 2:
        return parts[1]
    return ""

def _classify_te(class_token: str) -> bool:
    """True if `class_token` looks like a transposable-element class."""
    if not class_token:
        return False
    t = class_token.strip().lower().rstrip("?")
    superfam = t.split("/", 1)[0]
    if superfam in _NON_TE_CLASSES or t in _NON_TE_CLASSES:
        return False
    return superfam in _TE_SUPERFAMILIES

def _is_te_repeat_record(row: pd.Series) -> bool:
    """Return True for TE-like repeat records, excluding simple/non-TE repeats.

    Order of precedence:
      1. RepeatMasker Target= class token in attributes (most reliable).
      2. The GFF `feature` column itself if it already encodes a TE class
         (some pipelines emit e.g. feature=LTR/Gypsy directly).
    """
    cls = _repeat_target_class(str(row.get("attributes", "")))
    if cls:
        return _classify_te(cls)
    return _classify_te(str(row.get("feature", "")))

def _filter_te_repeats(repeat_gff: Optional[pd.DataFrame]) -> Optional[pd.DataFrame]:
    """Keep only transposable-element-like records from a repeat GFF dataframe."""
    if repeat_gff is None:
        return None
    if repeat_gff.empty:
        return repeat_gff
    mask = repeat_gff.apply(_is_te_repeat_record, axis=1)
    return repeat_gff.loc[mask].copy()

def list_geves(marker: pd.DataFrame, bed: pd.DataFrame) -> List[str]:
    return sorted(
        set(marker.get("geve_name", pd.Series(dtype=str)).dropna().astype(str))
        | set(bed.get("geve_name", pd.Series(dtype=str)).dropna().astype(str)),
        key=_natural_key,
    )

def infer_output_prefix(geve_names: List[str]) -> str:
    """Prefix = substring before the first underscore in a GEVE name."""
    for name in geve_names:
        if "_" in name:
            return name.split("_", 1)[0]
        if name:
            return name
    return "GEVE"

def bp_tick_formatter(x: float) -> str:
    if abs(x) >= 1_000_000:
        return f"{x/1_000_000:.2f} Mb"
    if abs(x) >= 1_000:
        return f"{x/1_000:.0f} Kb"
    return f"{int(x)} bp"

def nice_ticks(start: int, end: int, n: int = 6) -> List[int]:
    span = max(1, end - start)
    raw = span / max(1, (n - 1))
    mag = 10 ** math.floor(math.log10(raw))
    step = mag
    for mult in [1, 2, 2.5, 5, 10]:
        step = mag * mult
        if raw <= step:
            break
    tick0 = math.ceil(start / step) * step
    ticks = []
    v = tick0
    while v <= end:
        ticks.append(int(v))
        v += step
    if start not in ticks:
        ticks = [start] + ticks
    if end not in ticks:
        ticks.append(end)
    return sorted(set(ticks))

def build_hallmark_palette(names: List[str]) -> Dict[str, str]:
    fixed = {
        "A32":   "#1f77b4",
        "D5":    "#ff7f0e",
        "PolB":  "#7f7f7f",
        "RNAPL": "#9467bd",
        "RNAPS": "#8c564b",
        "RNR":   "#e377c2",
        "SFII":  "#17becf",
        "VLTF3": "#bcbd22",
        "mRNAc": "#2ca02c",
        "mcp":   "#d62728",
    }
    missing = [n for n in names if n not in fixed]
    if missing:
        cmap = plt.get_cmap("tab20")
        for i, name in enumerate(missing):
            fixed[name] = mcolors.to_hex(cmap(i % 20))
    return {k: fixed[k] for k in names}

def draw_orf_track(ax, feats: pd.DataFrame, x0: int, x1: int) -> Dict[str, str]:
    """Hallmark ORFs and TIR markers on a single row."""
    hall = feats[feats["feature"] == "hallmark"].copy()
    hall = hall.sort_values(["start", "end"], kind="mergesort")
    names = sorted(hall["name"].fillna(".").astype(str).unique(), key=_natural_key)
    palette = build_hallmark_palette([n for n in names if n != "."])

    feature_lc = feats["feature"].fillna("").astype(str).str.lower()
    tir_left = feats[feature_lc.eq("tir_left")].copy()
    tir_right = feats[feature_lc.eq("tir_right")].copy()
    tir_left = tir_left[(tir_left["end"].astype(float) >= x0) & (tir_left["start"].astype(float) <= x1)]
    tir_right = tir_right[(tir_right["end"].astype(float) >= x0) & (tir_right["start"].astype(float) <= x1)]
    tir_left = tir_left.sort_values(["start", "end"], kind="mergesort")
    tir_right = tir_right.sort_values(["start", "end"], kind="mergesort")

    ax.set_xlim(x0, x1)
    ax.set_ylim(0, 1)
    y = 0.5

    tir_color = "black"
    tir_marker_size = 7
    tir_label_size = 8
    for rows, marker, label_offset, ha in [
        (tir_left, ">", (-6, 0), "right"),
        (tir_right, "<", (6, 0), "left"),
    ]:
        for _, r in rows.iterrows():
            start = int(r["start"])
            end = int(r["end"])
            cx = (start + end) / 2.0
            ax.plot(
                [start, end], [y, y],
                color=tir_color, lw=1.8, solid_capstyle="round",
                zorder=4,
            )
            ax.plot(
                cx, y,
                marker=marker, color=tir_color, markersize=tir_marker_size,
                markeredgewidth=0.55, markeredgecolor=PLOT_COLORS["hallmark_edge"],
                zorder=5, linestyle="none",
            )
            ax.annotate(
                "TIR", xy=(cx, y), xytext=label_offset, textcoords="offset points",
                ha=ha, va="center", fontsize=tir_label_size, color=tir_color,
                zorder=6, clip_on=True,
            )

    for _, r in hall.iterrows():
        color = palette.get(str(r["name"]), PLOT_COLORS["hallmark_default"])
        cx = (int(r["start"]) + int(r["end"])) / 2.0
        ax.plot(
            cx, y,
            marker="o", color=color, markersize=PLOT_STYLE["hallmark_star_size"],
            markeredgewidth=PLOT_STYLE["hallmark_star_edge_width"],
            markeredgecolor="none",
            zorder=7, linestyle="none",
        )

    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylabel(
        "Markers", rotation=0, ha="right", va="center", fontsize=10, labelpad=8
    )
    for s in ["left", "right", "bottom", "top"]:
        ax.spines[s].set_visible(False)
    return palette

def _format_axis_end_tick(v: float) -> str:
    if not np.isfinite(v):
        return ""
    av = abs(v)
    if av >= 100:
        return f"{v:.0f}"
    if av >= 10:
        return f"{v:.1f}"
    return f"{v:.2f}".rstrip("0").rstrip(".")

def draw_line_track(
    ax,
    x: np.ndarray,
    y: np.ndarray,
    color: str,
    label: str,
    baseline: Optional[float] = None,
    fill: bool = True,
    linewidth: float = 1.2,
    background_color: str = "#fbfbfb",
    ylim: Optional[Tuple[float, float]] = None,
    midline: Optional[float] = None,
    ytick_values: Optional[List[float]] = None,
) -> None:
    valid = np.isfinite(y)
    ax.set_xlim(x.min(), x.max())
    ax.set_facecolor(background_color)
    if ylim is not None:
        ax.set_ylim(*ylim)
    elif valid.any():
        yv = y[valid]
        ymin = np.nanmin(yv)
        ymax = np.nanmax(yv)
        if baseline is not None:
            ymin = min(ymin, baseline)
            ymax = max(ymax, baseline)
        if ymin == ymax:
            ymin -= 1.0
            ymax += 1.0
        pad = (ymax - ymin) * 0.12
        ax.set_ylim(ymin - pad, ymax + pad)
    else:
        ax.set_ylim(0, 1)

    if midline is not None:
        ax.axhline(midline, color="#d9d9d9", lw=0.7, ls="-", zorder=1)
    if baseline is not None:
        ax.axhline(baseline, color="#b0b0b0", lw=0.8, ls="--", zorder=1)

    if valid.any():
        ax.plot(x[valid], y[valid], color=color, lw=linewidth, zorder=2)
        if fill:
            fill_base = baseline if baseline is not None else 0.0
            ax.fill_between(x[valid], y[valid], fill_base, color=color, alpha=0.50, zorder=1)

    ymin, ymax = ax.get_ylim()
    if ytick_values is None:
        ticks = [ymin, ymax]
    else:
        ticks = [ymin]
        ticks.extend(float(t) for t in ytick_values if ymin <= float(t) <= ymax)
        ticks.append(ymax)
        deduped = []
        for t in ticks:
            if not any(np.isclose(t, old_t) for old_t in deduped):
                deduped.append(t)
        ticks = deduped

    ax.set_yticks(ticks)
    ax.set_yticklabels([_format_axis_end_tick(t) for t in ticks], fontsize=8, color="black")
    # ylabel is forced black regardless of line colour
    ax.set_ylabel(label, rotation=0, ha="right", va="center", fontsize=10, labelpad=8, color="black")
    ax.tick_params(axis="y", labelsize=8, colors="black", length=8, width=1.0, direction="out")
    ax.set_xticks([])
    for s in ["top", "right", "bottom"]:
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)

def draw_bar_track(
    ax,
    starts: np.ndarray,
    ends: np.ndarray,
    values: np.ndarray,
    color: str,
    label: str,
    background_color: str = "#fbfbfb",
) -> None:
    """Window-based bar track for annotation coverage distribution."""
    clean = np.array(values, dtype=float)
    valid = np.isfinite(clean)
    ax.set_xlim(starts.min(), ends.max())
    ax.set_facecolor(background_color)

    ymax = float(np.nanmax(clean[valid])) if valid.any() else 1.0
    if not np.isfinite(ymax) or ymax <= 0:
        ymax = 1.0
    ymax = min(1.0, max(0.05, ymax))
    ax.set_ylim(0, ymax * 1.12)

    widths = np.maximum(1.0, ends - starts + 1)
    if valid.any():
        ax.bar(
            starts[valid], clean[valid], width=widths[valid], align="edge",
            color=color, edgecolor="none", alpha=0.85, zorder=2,
        )

    ymin, ymax_axis = ax.get_ylim()
    ticks = [ymin, ymax_axis]
    ax.set_yticks(ticks)
    ax.set_yticklabels([_format_axis_end_tick(t) for t in ticks], fontsize=8, color="black")
    ax.set_ylabel(label, rotation=0, ha="right", va="center", fontsize=10, labelpad=8, color="black")
    ax.tick_params(axis="y", labelsize=8, colors="black", length=8, width=1.0, direction="out")
    ax.set_xticks([])
    for s in ["top", "right", "bottom"]:
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)

def draw_strip_track(
    ax,
    starts: np.ndarray,
    ends: np.ndarray,
    values: np.ndarray,
    cmap: str,
    label: str,
    vmax: Optional[float] = None,
    vmin: Optional[float] = 0.0,
) -> None:
    """1-D heatmap strip; each window is a coloured rectangle."""
    ax.set_xlim(starts.min(), ends.max())
    ax.set_facecolor("#ffffff")
    clean = np.array(values, dtype=float)
    valid = np.isfinite(clean)
    if vmax is None:
        vmax = np.nanpercentile(clean[valid], 98) if valid.any() else 1.0
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = 1.0
    if vmin is None:
        vmin = np.nanmin(clean[valid]) if valid.any() else 0.0

    base_colors = {"Greens": PLOT_COLORS["gvog"], "Purples": PLOT_COLORS["pfam"],
                   "Intron": PLOT_COLORS["intron"], "Exon": PLOT_COLORS["exon"],
                   "Repeat": PLOT_COLORS["repeat"]}
    color = base_colors.get(cmap, PLOT_COLORS["hallmark_default"])

    for s, e, v in zip(starts, ends, clean):
        if not np.isfinite(v) or v <= 0:
            continue
        ax.add_patch(patches.Rectangle(
            (float(s), 0), float(e) - float(s) + 1, 1,
            facecolor=color, edgecolor="none", alpha=0.85,
        ))

    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylabel(label, rotation=0, ha="right", va="center", fontsize=10, labelpad=8, color="black")
    for s in ["left", "right", "bottom", "top"]:
        ax.spines[s].set_visible(False)

def get_geve_interval(feats: pd.DataFrame, x0: int, x1: int) -> Tuple[int, int]:
    geve = feats[feats["feature"] == "GEVE"].sort_values("start").head(1)
    if geve.empty:
        return x0, x1
    r = geve.iloc[0]
    return int(r["start"]), int(r["end"])

def highlight_geve_region(ax, geve_start: int, geve_end: int) -> None:
    ax.axvspan(
        geve_start, geve_end,
        facecolor=PLOT_COLORS["geve_highlight"],
        alpha=PLOT_STYLE["geve_highlight_alpha"],
        edgecolor="none",
        zorder=0,
    )

def add_legends(fig, hallmark_palette: Dict[str, str]) -> None:
    if not hallmark_palette:
        return
    handles = []
    for name in sorted(hallmark_palette, key=_natural_key):
        handles.append(plt.Line2D(
            [0], [0], marker="o", linestyle="none",
            markerfacecolor=hallmark_palette[name],
            markeredgecolor="none",
            markeredgewidth=PLOT_STYLE["hallmark_star_edge_width"],
            markersize=10, label=name,
        ))
    fig.legend(
        handles=handles,
        title="Hallmark genes",
        loc="lower center",
        bbox_to_anchor=(0.5, PLOT_STYLE["legend_y"]),
        ncol=len(handles),
        frameon=False,
        fontsize=9,
        title_fontsize=10,
        handletextpad=0.5,
        columnspacing=1.2,
        borderaxespad=0.0,
    )

def _coord_hover() -> str:
    """Compact hover text for coordinate review."""
    return "coord: %{x:,.0f}<extra></extra>"


def _feature_rows_in_view(m: pd.DataFrame, feature: str, x0: int, x1: int, prefix: bool = False) -> pd.DataFrame:
    if m.empty or "feature" not in m:
        return m.iloc[0:0].copy()
    f = m["feature"].astype(str).str.lower()
    target = feature.lower()
    mask = f.str.startswith(target) if prefix else f.eq(target)
    q = m[mask].copy()
    if q.empty:
        return q
    q = q[(q["end"].astype(float) >= x0) & (q["start"].astype(float) <= x1)]
    return q.sort_values(["start", "end"], kind="mergesort")

def _add_html_evidence_background(
    fig: go.Figure,
    geve_start: int,
    geve_end: int,
    color: str = "rgba(230,230,230,0.45)",
) -> None:
    """Draw the GEVE region as a true background on the HTML evidence track.
    """
    fig.add_trace(
        go.Scatter(
            x=[geve_start, geve_end, geve_end, geve_start, geve_start],
            y=[0, 0, 1, 1, 0],
            mode="lines",
            fill="toself",
            fillcolor=color,
            line=dict(width=0, color=color),
            hoverinfo="skip",
            showlegend=False,
            name="GEVE region",
        ),
        row=3, col=1,
    )

def _html_bp_tick_formatter(x: float) -> str:
    """HTML coordinate labels that always use explicit Kb/Mb units."""
    if not np.isfinite(x):
        return ""
    ax = abs(float(x))
    if ax >= 1_000_000:
        value = float(x) / 1_000_000.0
        return f"{value:.2f} Mb"
    value = float(x) / 1_000.0
    if ax < 10_000:
        return f"{value:.1f} Kb"
    return f"{value:.0f} Kb"

def _count_round_ticks(start: int, end: int, step: int) -> int:
    """Number of rounded tick positions inside [start, end]."""
    if step <= 0 or end <= start:
        return 0
    first = math.ceil(start / step) * step
    last = math.floor(end / step) * step
    if first > last:
        return 0
    return int((last - first) // step) + 1

def _html_coordinate_ticks(
    start: int,
    end: int,
    geve_start: Optional[int] = None,
    geve_end: Optional[int] = None,
    max_ticks: int = 13,
) -> List[int]:
    """Adaptive x-axis ticks for the HTML coordinate review plot.
    """
    start = int(start)
    end = int(end)
    if end < start:
        start, end = end, start

    span = max(1, end - start)
    if geve_start is not None and geve_end is not None:
        geve_len = max(1, int(geve_end) - int(geve_start) + 1)
    else:
        geve_len = span
    candidate_steps = [
        10_000, 20_000, 25_000, 50_000,
        100_000, 200_000, 250_000, 500_000,
        1_000_000, 2_000_000, 5_000_000, 10_000_000,
    ]

    preferred_step = None
    short_reference = min(span, geve_len)
    if short_reference <= 350_000:
        preferred_step = 50_000     
    elif short_reference <= 700_000:
        preferred_step = 100_000 
    elif short_reference <= 1_200_000:
        preferred_step = 200_000
    elif span <= 5_000_000:
        preferred_step = 500_000

    if preferred_step is not None:
        step = preferred_step
        if _count_round_ticks(start, end, step) > max_ticks:
            start_idx = candidate_steps.index(step) if step in candidate_steps else 0
            for cand in candidate_steps[start_idx + 1:]:
                count = _count_round_ticks(start, end, cand)
                if 2 <= count <= max_ticks:
                    step = cand
                    break
    else:
        step = candidate_steps[-1]
        for cand in candidate_steps:
            count = _count_round_ticks(start, end, cand)
            if 2 <= count <= max_ticks:
                step = cand
                break

    first = math.ceil(start / step) * step
    ticks = []
    v = first
    while v <= end:
        ticks.append(int(v))
        v += step

    if len(ticks) < 2:
        ticks = nice_ticks(start, end, n=5)
    return sorted(set(ticks))

def _add_html_window_strip(
    fig: go.Figure,
    starts: np.ndarray,
    ends: np.ndarray,
    values: np.ndarray,
    y0: float,
    y1: float,
    color: str,
    label: str,
) -> None:
    """Draw BED-window hit strips in the same style used by the PDF plot.
    """
    clean = np.array(values, dtype=float)
    valid = np.isfinite(clean) & (clean > 0)
    if not valid.any():
        return

    s = np.array(starts, dtype=float)[valid]
    e = np.array(ends, dtype=float)[valid]
    v = clean[valid]
    centers = (s + e) / 2.0
    widths = np.maximum(1.0, e - s + 1.0)
    customdata = np.column_stack([s, e, v])

    fig.add_trace(
        go.Bar(
            x=centers,
            y=np.full(len(centers), y1 - y0, dtype=float),
            base=np.full(len(centers), y0, dtype=float),
            width=widths,
            marker=dict(color=color, opacity=1.0, line=dict(width=0)),
            customdata=customdata,
            hovertemplate=(
                f"{label}: %{{customdata[2]:.0f}}<br>"
                "window: %{customdata[0]:,.0f}-%{customdata[1]:,.0f}<br>"
                "coord: %{x:,.0f}<extra></extra>"
            ),
            showlegend=False,
            name=label,
        ),
        row=3, col=1,
    )

def _add_html_feature_segments(
    fig: go.Figure,
    rows: pd.DataFrame,
    y: float,
    color: str,
    width: int,
) -> None:
    if rows.empty:
        return
    xvals = []
    yvals = []
    for s, e in rows[["start", "end"]].itertuples(index=False):
        xvals.extend([float(s), float(e), None])
        yvals.extend([y, y, None])
    fig.add_trace(
        go.Scatter(
            x=xvals,
            y=yvals,
            mode="lines",
            line=dict(color=color, width=width),
            hovertemplate=_coord_hover(),
            showlegend=False,
            connectgaps=False,
        ),
        row=3, col=1,
    )

def _add_html_hallmark_stars(fig: go.Figure, rows: pd.DataFrame, y: float, color: str, size: int = 12) -> None:
    if rows.empty:
        return
    xs = ((rows["start"].astype(float) + rows["end"].astype(float)) / 2.0).tolist()
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=[y] * len(xs),
            mode="markers",
            marker=dict(symbol="circle", size=size, color=color, line=dict(width=0)),
            hovertemplate=_coord_hover(),
            showlegend=False,
        ),
        row=3, col=1,
    )

def plot_one_geve_html(marker: pd.DataFrame, bed: pd.DataFrame, geve_name: str) -> go.Figure:
    """Build a lightweight Plotly coordinate-review figure for one GEVE."""
    m = marker[marker["geve_name"] == geve_name].copy()
    b = bed[bed["geve_name"] == geve_name].copy().sort_values("window_start")
    if m.empty or b.empty:
        raise ValueError(f"Missing data for {geve_name}")

    contig = (
        str(m["contig"].dropna().iloc[0])
        if "contig" in m and not m["contig"].dropna().empty
        else str(b["contig_id"].dropna().iloc[0])
    )
    x0 = int(min(b["window_start"].min(), m["start"].dropna().min()))
    x1 = int(max(b["window_end"].max(), m["end"].dropna().max()))
    geve_start, geve_end = get_geve_interval(m, x0, x1)

    centers = (b["window_start"].to_numpy(dtype=float) + b["window_end"].to_numpy(dtype=float)) / 2.0
    viral = b["rolling_score_mean"].to_numpy(dtype=float)
    gc = b["gc"].to_numpy(dtype=float)

    tir_left = _feature_rows_in_view(m, "TIR_left", x0, x1)
    tir_right = _feature_rows_in_view(m, "TIR_right", x0, x1)
    hallmark = _feature_rows_in_view(m, "hallmark", x0, x1)
    starts = b["window_start"].to_numpy(dtype=float)
    ends = b["window_end"].to_numpy(dtype=float)

    fig = make_subplots(
        rows=3,
        cols=1,
        shared_xaxes=True,
        row_heights=[0.44, 0.20, 0.36],
        vertical_spacing=0.09,
        subplot_titles=("Viral score", "GC (%)", "Evidence marker"),
    )

    fig.add_trace(go.Scatter(
        x=centers, y=viral, mode="lines",
        line=dict(color=PLOT_COLORS["viral_score"], width=1.8),
        fill="tozeroy", hovertemplate=_coord_hover(), showlegend=False,
        connectgaps=True, name="Viral score"), row=1, col=1)
    fig.add_hline(y=0, line_dash="dash", line_width=1, line_color="#999999", row=1, col=1)

    fig.add_trace(go.Scatter(
        x=centers, y=gc, mode="lines",
        line=dict(color=PLOT_COLORS["gc_content"], width=1.6),
        hovertemplate=_coord_hover(), showlegend=False,
        connectgaps=True, name="GC (%)"), row=2, col=1)
    fig.add_hline(y=50, line_width=1, line_color=PLOT_COLORS["gc_midline"], row=2, col=1)
    fig.update_yaxes(range=[0, 100], row=2, col=1)

    for rr in [1, 2]:
        fig.add_vrect(
            x0=geve_start, x1=geve_end,
            fillcolor=PLOT_COLORS["geve_highlight"], opacity=0.55,
            line_width=0, layer="below", row=rr, col=1,
        )
    _add_html_evidence_background(fig, geve_start, geve_end)

    GVOG_Y0, GVOG_Y1 = 0.72, 0.93
    PFAM_Y0, PFAM_Y1 = 0.47, 0.68
    HALLMARK_Y = 0.32
    TIR_Y = 0.13
    EV_THIN = 4
    HALLMARK_SIZE = 9
    TIR_MARKER_SIZE = 12
    TIR_LABEL_SIZE = 11
    EV_LEGEND_WIDTH = 18

    _add_html_window_strip(
        fig, starts, ends, b["gvog_hits"].to_numpy(dtype=float),
        GVOG_Y0, GVOG_Y1, PLOT_COLORS["gvog"], "GVOG hits",
    )
    _add_html_window_strip(
        fig, starts, ends, b["pfam_hits"].to_numpy(dtype=float),
        PFAM_Y0, PFAM_Y1, PLOT_COLORS["pfam"], "Pfam hits",
    )
    _add_html_hallmark_stars(fig, hallmark, HALLMARK_Y, PLOT_COLORS["hallmark_default"], size=HALLMARK_SIZE)
    _add_html_feature_segments(fig, tir_left, TIR_Y, "#238b45", EV_THIN)
    _add_html_feature_segments(fig, tir_right, TIR_Y, "#238b45", EV_THIN)

    if not tir_left.empty:
        xs = ((tir_left["start"].astype(float) + tir_left["end"].astype(float)) / 2.0).tolist()
        fig.add_trace(go.Scatter(
            x=xs,
            y=[TIR_Y] * len(xs),
            mode="markers+text",
            marker=dict(symbol="triangle-right", size=TIR_MARKER_SIZE, color="#238b45", line=dict(color="black", width=0.7)),
            text=["TIR"] * len(xs),
            textposition="middle right",
            textfont=dict(size=TIR_LABEL_SIZE, color="#238b45"),
            hovertemplate=_coord_hover(),
            showlegend=False,
        ), row=3, col=1)
    if not tir_right.empty:
        xs = ((tir_right["start"].astype(float) + tir_right["end"].astype(float)) / 2.0).tolist()
        fig.add_trace(go.Scatter(
            x=xs,
            y=[TIR_Y] * len(xs),
            mode="markers+text",
            marker=dict(symbol="triangle-left", size=TIR_MARKER_SIZE, color="#238b45", line=dict(color="black", width=0.7)),
            text=["TIR"] * len(xs),
            textposition="middle left",
            textfont=dict(size=TIR_LABEL_SIZE, color="#238b45"),
            hovertemplate=_coord_hover(),
            showlegend=False,
        ), row=3, col=1)

    fig.add_trace(go.Scatter(
        x=[None], y=[None], mode="markers",
        marker=dict(
            symbol="square", size=12,
            color=PLOT_COLORS["geve_highlight"],
            line=dict(width=0.5, color="#bdbdbd"),
        ),
        name="GEVE", hoverinfo="skip", showlegend=True,
    ), row=3, col=1)
    fig.add_trace(go.Scatter(
        x=[None, None], y=[None, None], mode="lines",
        line=dict(color=PLOT_COLORS["gvog"], width=EV_LEGEND_WIDTH),
        name="GVOG", hoverinfo="skip", showlegend=True,
    ), row=3, col=1)
    fig.add_trace(go.Scatter(
        x=[None, None], y=[None, None], mode="lines",
        line=dict(color=PLOT_COLORS["pfam"], width=EV_LEGEND_WIDTH),
        name="Pfam", hoverinfo="skip", showlegend=True,
    ), row=3, col=1)
    fig.add_trace(go.Scatter(
        x=[None], y=[None], mode="markers",
        marker=dict(symbol="circle", size=HALLMARK_SIZE, color=PLOT_COLORS["hallmark_default"], line=dict(width=0)),
        name="Hallmark", hoverinfo="skip", showlegend=True,
    ), row=3, col=1)

    geve_len = int(geve_end - geve_start + 1)
    fig.update_layout(
        title=f"{geve_name} ({_html_bp_tick_formatter(geve_len)}) | {contig} ({geve_start:,}-{geve_end:,})",
        height=720, template="plotly_white", hovermode="closest", dragmode="zoom",
        barmode="overlay",
        margin=dict(l=70, r=30, t=85, b=100), uirevision=geve_name,
        legend=dict(
            orientation="h", yanchor="top", y=-0.17, xanchor="center", x=0.5,
            bgcolor="rgba(255,255,255,0)", borderwidth=0, font=dict(size=11),
            traceorder="normal",
            itemsizing="constant",
        ),
    )
    fig.update_xaxes(range=[x0, x1], title_text="Genomic coordinate (Kb/Mb)", row=3, col=1)
    for rr in range(1, 4):
        fig.update_xaxes(showgrid=True, gridcolor="#eeeeee", row=rr, col=1)
        fig.update_yaxes(showgrid=False, row=rr, col=1)
    x_ticks = _html_coordinate_ticks(x0, x1, geve_start=geve_start, geve_end=geve_end, max_ticks=13)
    x_tick_text = [_html_bp_tick_formatter(t) for t in x_ticks]
    fig.update_xaxes(
        showline=True, linewidth=1.0, linecolor="black", mirror=False,
        ticks="outside", ticklen=3, tickwidth=0.7, tickcolor="black",
        tickmode="array", tickvals=x_ticks, ticktext=x_tick_text,
        showexponent="none", exponentformat="none",
        row=3, col=1,
    )
    fig.update_yaxes(range=[0, 1], showticklabels=False, title_text="", row=3, col=1)
    return fig

def write_html_report(marker: pd.DataFrame, bed: pd.DataFrame, geves: List[str], outpath: Path) -> None:
    """Write one dropdown-based Plotly HTML coordinate-review report."""
    figures = {}
    options = []
    for geve_name in geves:
        fig = plot_one_geve_html(marker, bed, geve_name)
        figures[geve_name] = fig.to_plotly_json()
        options.append(f'<option value="{geve_name}">{geve_name}</option>')

    figures_json = json.dumps(figures, cls=PlotlyJSONEncoder, separators=(",", ":"))
    options_html = "\n".join(options)
    html = f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>findGEVE coordinate review</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
body {{ font-family: Arial, sans-serif; margin: 20px; background: #ffffff; color: #222; }}
h1 {{ font-size: 22px; margin: 0 0 8px 0; }}
.controls {{ display: flex; flex-wrap: wrap; gap: 12px; align-items: center; margin: 12px 0 10px 0; }}
select {{ font-size: 14px; padding: 4px 8px; min-width: 260px; }}
button {{ font-size: 14px; padding: 5px 10px; cursor: pointer; }}
#coordBox {{ font-size: 16px; font-weight: bold; padding: 6px 10px; background: #f3f3f3; border-radius: 4px; }}
.note {{ color: #555; max-width: 950px; line-height: 1.35; }}
#plot {{ width: 100%; min-height: 680px; }}
</style>
</head>
<body>
<h1>findGEVE coordinate review</h1>
<div class="note">This page is designed for GEVE boundary review. Hover over a data point to view its coordinate; use this coordinate if you want to adjust the GEVE boundary.</div>
<div class="controls">
  <label for="geveSelect"><b>GEVE:</b></label>
  <select id="geveSelect">{options_html}</select>
  <span id="coordBox">Coordinate: none</span>
  <button id="copyButton" type="button">Copy</button>
</div>
<div id="plot"></div>
<script>
const figures = {figures_json};
let currentCoord = "";
function renderSelected() {{
  const key = document.getElementById('geveSelect').value;
  const fig = figures[key];
  const config = {{responsive: true, displaylogo: false, scrollZoom: true}};
  Plotly.react('plot', fig.data, fig.layout, config).then(() => {{
    const plotDiv = document.getElementById('plot');
    if (plotDiv.removeAllListeners) {{ plotDiv.removeAllListeners('plotly_click'); }}
    plotDiv.on('plotly_click', function(data) {{
      if (!data.points || data.points.length === 0) return;
      const rounded = Math.round(Number(data.points[0].x));
      if (!Number.isFinite(rounded)) return;
      currentCoord = String(rounded);
      document.getElementById('coordBox').textContent = 'Clicked coordinate: ' + rounded.toLocaleString();
      if (navigator.clipboard && navigator.clipboard.writeText) {{ navigator.clipboard.writeText(currentCoord).catch(() => {{}}); }}
    }});
  }});
}}
document.getElementById('geveSelect').addEventListener('change', renderSelected);
document.getElementById('copyButton').addEventListener('click', function() {{
  if (!currentCoord) return;
  if (navigator.clipboard && navigator.clipboard.writeText) {{ navigator.clipboard.writeText(currentCoord).catch(() => {{}}); }}
}});
renderSelected();
</script>
</body>
</html>
"""
    outpath.write_text(html, encoding="utf-8")

def plot_one_geve(marker: pd.DataFrame, bed: pd.DataFrame, geve_name: str,
                  gene_gff: Optional[pd.DataFrame] = None,
                  repeat_gff: Optional[pd.DataFrame] = None):
    """Build and return a Figure for one GEVE."""
    m = marker[marker["geve_name"] == geve_name].copy()
    b = bed[bed["geve_name"] == geve_name].copy().sort_values("window_start")
    if m.empty or b.empty:
        raise ValueError(f"Missing data for {geve_name}")

    contig = (
        str(m["contig"].dropna().iloc[0])
        if "contig" in m and not m["contig"].dropna().empty
        else str(b["contig_id"].dropna().iloc[0])
    )

    x0 = int(min(b["window_start"].min(), m["start"].dropna().min()))
    x1 = int(max(b["window_end"].max(), m["end"].dropna().max()))

    centers = (b["window_start"].to_numpy(dtype=float) + b["window_end"].to_numpy(dtype=float)) / 2.0
    starts  = b["window_start"].to_numpy(dtype=float)
    ends    = b["window_end"].to_numpy(dtype=float)

    track_specs = [
        ("line", "viral_score"),
        ("line", "gc"),
    ]
    if _has_feature(gene_gff, "intron"):
        track_specs.append(("bar", "intron"))
    if _has_feature(gene_gff, "exon"):
        track_specs.append(("bar", "exon"))
    if repeat_gff is not None and not repeat_gff.empty:
        track_specs.append(("strip", "repeat"))
    track_specs.extend([
        ("strip", "gvog"),
        ("strip", "pfam"),
        ("orf", "hallmark"),
    ])

    height_by_key = {
        "viral_score": 1.12, "gc": 0.50,
        "intron": 0.32, "exon": 0.32,
        "repeat": 0.26, "gvog": 0.26, "pfam": 0.26, "hallmark": 0.26,
    }
    fig_height = FIG_HEIGHT + 0.32 * max(0, len(track_specs) - 5)
    fig = plt.figure(figsize=(FIG_WIDTH, fig_height), constrained_layout=False)
    gs = fig.add_gridspec(
        len(track_specs), 1,
        height_ratios=[height_by_key[key] for _, key in track_specs],
        hspace=0.14,
    )
    axs = [fig.add_subplot(gs[i, 0]) for i in range(len(track_specs))]

    geve_start, geve_end = get_geve_interval(m, x0, x1)

    palette = {}
    for ax, (_, key) in zip(axs, track_specs):
        if key == "viral_score":
            draw_line_track(
                ax, centers, b["rolling_score_mean"].to_numpy(dtype=float),
                PLOT_COLORS["viral_score"], "Viral score",
                baseline=0.0, fill=True, linewidth=1.2, ytick_values=[0],
            )
        elif key == "gc":
            draw_line_track(
                ax, centers, b["gc"].to_numpy(dtype=float),
                PLOT_COLORS["gc_content"], "GC (%)",
                baseline=None, fill=False, linewidth=1.1,
                ylim=(0, 100), midline=50, ytick_values=[0, 50, 100],
            )
        elif key == "intron":
            vals = _window_overlap_fraction(b, gene_gff, contig, {"intron"})
            draw_bar_track(ax, starts, ends, vals, PLOT_COLORS["intron"], "Introns")
        elif key == "exon":
            vals = _window_overlap_fraction(b, gene_gff, contig, {"exon"})
            draw_bar_track(ax, starts, ends, vals, PLOT_COLORS["exon"], "Exons")
        elif key == "repeat":
            vals = _window_overlap_fraction(b, repeat_gff, contig, None)
            draw_strip_track(ax, starts, ends, vals, "Repeat", "TEs", vmin=0.0)
        elif key == "gvog":
            draw_strip_track(ax, starts, ends, b["gvog_hits"].to_numpy(dtype=float),
                             "Greens", "GVOG hits", vmin=0.0)
        elif key == "pfam":
            draw_strip_track(ax, starts, ends, b["pfam_hits"].to_numpy(dtype=float),
                             "Purples", "Pfam hits", vmin=0.0)
        elif key == "hallmark":
            palette = draw_orf_track(ax, m, x0, x1)

    for ax in axs:
        highlight_geve_region(ax, geve_start, geve_end)
        ax.set_xlim(x0, x1)

    ticks = nice_ticks(x0, x1, n=6)
    axs[-1].set_xticks(ticks)
    axs[-1].set_xticklabels([bp_tick_formatter(t) for t in ticks], fontsize=9)
    axs[-1].tick_params(axis="x", which="both", length=4, pad=5, bottom=True)
    axs[-1].spines["bottom"].set_visible(True)
    axs[-1].spines["bottom"].set_linewidth(0.6)

    geve_rows = m[m["feature"] == "GEVE"].copy()
    if not geve_rows.empty:
        geve_len = int((geve_rows["end"] - geve_rows["start"] + 1).max())
    else:
        geve_len = x1 - x0 + 1

    title = f"{geve_name} ({bp_tick_formatter(geve_len)})  |  {contig}"
    fig.suptitle(title, fontsize=16, y=PLOT_STYLE["title_y"], va="top")

    add_legends(fig, palette)
    fig.subplots_adjust(top=0.895, left=0.12, right=0.98, bottom=0.14)
    return fig

def main(argv: List[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Plot findGEVE windowed viral score/GC/domain tracks.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Gene GFF adds exon coverage and intron coverage below GC. If intron rows are absent, introns are inferred from exon gaps per transcript. Repeat GFF adds a TE-only strip below GC.",
    )
    parser.add_argument("-g", "--gff", dest="gene_gff", type=Path,
                        help="Gene annotation GFF3. Uses exon rows and explicit introns when present; otherwise infers introns from exon gaps per transcript.")
    parser.add_argument("-r", "--repeat", dest="repeat_gff", type=Path,
                        help="Repeat/TE annotation GFF3. Uses TE-like records only, e.g. LTR, LINE, SINE, DNA transposons, RC/Helitron, etc.")
    parser.add_argument("markerout", type=Path, help="findGEVE markerout file")
    parser.add_argument("bed", type=Path, help="findGEVE geve.bed window file")
    args = parser.parse_args(argv)

    markerout = args.markerout
    bed_path = args.bed_path if hasattr(args, "bed_path") else args.bed

    for label, path in [
        ("markerout", markerout),
        ("bed", bed_path),
        ("gene GFF", args.gene_gff),
        ("repeat GFF", args.repeat_gff),
    ]:
        if path is not None and not path.is_file():
            print(f"Error: {label} file not found: {path}", file=sys.stderr)
            return 1

    marker = _read_markerout(markerout)
    bed = _read_bed(bed_path)
    gene_gff = _prepare_gene_gff(_read_gff(args.gene_gff)) if args.gene_gff is not None else None
    repeat_gff = _filter_te_repeats(_read_gff(args.repeat_gff)) if args.repeat_gff is not None else None

    if marker.empty or bed.empty:
        print(f"Error: input file is empty (markerout={markerout}, bed={bed_path})", file=sys.stderr)
        return 1

    if gene_gff is not None and gene_gff.empty:
        print(f"Warning: gene GFF has no readable records: {args.gene_gff}", file=sys.stderr)
    if repeat_gff is not None and repeat_gff.empty:
        print(f"Warning: repeat GFF has no readable TE records after filtering: {args.repeat_gff}", file=sys.stderr)

    geves = list_geves(marker, bed)
    if not geves:
        print("Error: no GEVE entries found to plot.", file=sys.stderr)
        return 1

    prefix = infer_output_prefix(geves)
    outpath = Path.cwd() / f"{prefix}.plot.pdf"
    html_outpath = Path.cwd() / f"{prefix}.plot.html"

    with PdfPages(outpath) as pdf:
        for geve_name in geves:
            fig = plot_one_geve(marker, bed, geve_name, gene_gff=gene_gff, repeat_gff=repeat_gff)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
            print(f"Plotted {geve_name}")

    print(f"Wrote {outpath} ({len(geves)} page(s))")

    write_html_report(marker, bed, geves, html_outpath)
    print(f"Wrote {html_outpath} ({len(geves)} interactive figure(s))")
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
