#!/usr/bin/env python3
from __future__ import annotations

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

USAGE = "Usage: python findGEVE_plot.py <prefix.markerout> <prefix.geve.bed>"

PLOT_COLORS = {
    "viral_score": "#ec3028",
    "gc_content": "#fec44a",
    "gc_midline": "#d9d9d9",
    "gvog": "#1497a5",
    "pfam": "#756bb1",
    "geve_highlight": "#e6e6e6",
    "track_background": "#ffffff",
    "line_track_background": "#fbfbfb",
    "hallmark_default": "#636363",
    "hallmark_edge": "black",
}

PLOT_STYLE = {
    "hallmark_star_size": 13,
    "hallmark_star_edge_width": 0.55,
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
        return f"{x/1_000:.0f} kb"
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
    """Hallmark ORFs as star markers on a single row."""
    hall = feats[feats["feature"] == "hallmark"].copy()
    hall = hall.sort_values(["start", "end"], kind="mergesort")
    names = sorted(hall["name"].fillna(".").astype(str).unique(), key=_natural_key)
    palette = build_hallmark_palette([n for n in names if n != "."])

    ax.set_xlim(x0, x1)

    if hall.empty:
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_ylabel(
            "Hallmark ORFs", rotation=0, ha="right", va="center", fontsize=10, labelpad=8
        )
        for s in ["left", "right", "bottom", "top"]:
            ax.spines[s].set_visible(False)
        return {}

    ax.set_ylim(0, 1)
    y = 0.5
    for _, r in hall.iterrows():
        color = palette.get(str(r["name"]), PLOT_COLORS["hallmark_default"])
        cx = (int(r["start"]) + int(r["end"])) / 2.0
        ax.plot(
            cx, y,
            marker="*", color=color, markersize=PLOT_STYLE["hallmark_star_size"],
            markeredgewidth=PLOT_STYLE["hallmark_star_edge_width"],
            markeredgecolor=PLOT_COLORS["hallmark_edge"],
            zorder=5, linestyle="none",
        )

    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylabel(
        "Hallmark ORFs", rotation=0, ha="right", va="center", fontsize=10, labelpad=8
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

    base_colors = {"Greens": PLOT_COLORS["gvog"], "Purples": PLOT_COLORS["pfam"]}
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
            [0], [0], marker="*", linestyle="none",
            markerfacecolor=hallmark_palette[name],
            markeredgecolor=PLOT_COLORS["hallmark_edge"],
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


def plot_one_geve(marker: pd.DataFrame, bed: pd.DataFrame, geve_name: str):
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

    # Rows: 0 viral score, 1 GC%, 2 GVOG, 3 Pfam, 4 hallmark ORFs + x-axis
    fig = plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT), constrained_layout=False)
    gs = fig.add_gridspec(
        5, 1,
        height_ratios=[1.12, 0.50, 0.35, 0.35, 0.26],
        hspace=0.14,
    )
    axs = [fig.add_subplot(gs[i, 0]) for i in range(5)]

    geve_start, geve_end = get_geve_interval(m, x0, x1)

    draw_line_track(
        axs[0], centers, b["rolling_score_mean"].to_numpy(dtype=float),
        PLOT_COLORS["viral_score"], "Viral score",
        baseline=0.0, fill=True, linewidth=1.2, ytick_values=[0],
    )
    draw_line_track(
        axs[1], centers, b["gc"].to_numpy(dtype=float),
        PLOT_COLORS["gc_content"], "GC (%)",
        baseline=None, fill=False, linewidth=1.1,
        ylim=(0, 100), midline=50, ytick_values=[0, 50, 100],
    )
    draw_strip_track(axs[2], starts, ends, b["gvog_hits"].to_numpy(dtype=float),
                     "Greens", "GVOG hits", vmin=0.0)
    draw_strip_track(axs[3], starts, ends, b["pfam_hits"].to_numpy(dtype=float),
                     "Purples", "Pfam hits", vmin=0.0)
    palette = draw_orf_track(axs[4], m, x0, x1)

    for ax in axs:
        highlight_geve_region(ax, geve_start, geve_end)
        ax.set_xlim(x0, x1)

    ticks = nice_ticks(x0, x1, n=6)
    axs[4].set_xticks(ticks)
    axs[4].set_xticklabels([bp_tick_formatter(t) for t in ticks], fontsize=9)
    axs[4].tick_params(axis="x", which="both", length=4, pad=5, bottom=True)
    axs[4].spines["bottom"].set_visible(True)
    axs[4].spines["bottom"].set_linewidth(0.6)

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
    if len(argv) == 1 and argv[0].lower() in ("help", "-h", "--help"):
        print(USAGE)
        return 0
    if len(argv) != 2:
        print(USAGE, file=sys.stderr)
        return 1

    markerout = Path(argv[0])
    bed_path  = Path(argv[1])

    if not markerout.is_file():
        print(f"Error: markerout file not found: {markerout}", file=sys.stderr)
        return 1
    if not bed_path.is_file():
        print(f"Error: bed file not found: {bed_path}", file=sys.stderr)
        return 1

    marker = _read_markerout(markerout)
    bed    = _read_bed(bed_path)

    if marker.empty or bed.empty:
        print(f"Error: input file is empty (markerout={markerout}, bed={bed_path})", file=sys.stderr)
        return 1

    geves = list_geves(marker, bed)
    if not geves:
        print("Error: no GEVE entries found to plot.", file=sys.stderr)
        return 1

    prefix  = infer_output_prefix(geves)
    outpath = Path.cwd() / f"{prefix}.plot.pdf"

    with PdfPages(outpath) as pdf:
        for geve_name in geves:
            fig = plot_one_geve(marker, bed, geve_name)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
            print(f"Plotted {geve_name}")

    print(f"Wrote {outpath} ({len(geves)} page(s))")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
