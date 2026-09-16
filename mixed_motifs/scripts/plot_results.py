#!/usr/bin/env python3
"""Create publication-friendly summaries of mixed_motifs C++ result tables.

All displayed groups use members.csv; no nearest-neighbor reconstruction or
inferred cell identities are used. Each specimen receives its own spatial map.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import tempfile
import textwrap

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "mixed_motifs_mpl"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from matplotlib.patches import Circle, Patch
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

BLUE = "#245B86"
TEAL = "#127C72"
ORANGE = "#D67A26"
GRAY = "#C7CDD3"
INK = "#253646"

SCHEMAS = {
    "cells": "cell_index cell_id x y cell_type sample stratum component",
    "neighborhoods": "neighborhood_id center_index x y sample size radius n_types dominance entropy effective_types cross_edge_fraction eligible",
    "motifs": "motif_id composition support raw_matches n_regions n_samples expected_support fold_enrichment excess_support p_value p_value_holm p_search significant inference_scope",
    "motif_composition": "motif_id cell_type count",
    "occurrences": "motif_id occurrence_id neighborhood_id center_index x y sample radius",
    "members": "motif_id occurrence_id cell_index",
    "null": "permutation max_support eligible_neighborhoods",
    "motif_null": "permutation motif_id support",
}
TEXT_COLUMNS = {"cell_id", "cell_type", "sample", "stratum", "component", "motif_id", "occurrence_id", "composition", "inference_scope"}


def setup_style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 12,
        "axes.titleweight": "bold", "axes.labelsize": 10, "axes.spines.top": False,
        "axes.spines.right": False, "axes.edgecolor": "#A5AFB8", "axes.labelcolor": INK,
        "text.color": INK, "xtick.color": INK, "ytick.color": INK,
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
        "pdf.fonttype": 42, "ps.fonttype": 42, "legend.frameon": False,
        "grid.color": "#E3E7EA", "grid.alpha": 0.8, "axes.axisbelow": True,
    })


def load_tables(results):
    try:
        run = json.loads((results / "run.json").read_text())
    except (OSError, ValueError) as exc:
        raise ValueError(f"Cannot read {results / 'run.json'}: {exc}") from exc
    if isinstance(run.get("config"), dict):
        run = {**run["config"], **run}
    tables = {}
    for name, schema in SCHEMAS.items():
        path = results / f"{name}.csv"
        if not path.exists():
            raise ValueError(f"Required result table is missing: {path}")
        try:
            df = pd.read_csv(path, dtype=str, keep_default_na=False)
        except pd.errors.EmptyDataError:
            df = pd.DataFrame(columns=schema.split())
        missing = set(schema.split()) - set(df.columns)
        if missing:
            raise ValueError(f"{path.name} is missing columns: {', '.join(sorted(missing))}")
        for col in df:
            if col not in TEXT_COLUMNS and col not in {"eligible", "significant"}:
                df[col] = pd.to_numeric(df[col], errors="coerce")
        for col in ("eligible", "significant"):
            if col in df:
                df[col] = df[col].astype(str).str.lower().isin(["true", "1", "yes"])
        tables[name] = df
    cells = tables["cells"]
    if cells.cell_index.duplicated().any():
        raise ValueError("cells.csv has duplicate cell_index values; atlas identities are ambiguous")
    for col in ("x", "y", "cell_index"):
        if cells[col].isna().any() or not np.isfinite(cells[col].to_numpy(float)).all():
            raise ValueError(f"cells.csv contains invalid {col} values")
    missing_members = set(tables["members"].cell_index) - set(cells.cell_index)
    if missing_members:
        raise ValueError("members.csv refers to cells missing from cells.csv")
    return run, tables


def palette_for(labels):
    """Use a stable shared palette, extending the categorical colors when needed."""
    labels = sorted(set(labels))
    base = list(plt.get_cmap("tab20").colors)
    # Alternate the dark tab20 colors before the corresponding lighter shades.
    base = base[::2] + base[1::2]
    colors = []
    for i in range(len(labels)):
        if i < len(base):
            colors.append(to_hex(base[i]))
        else:
            colors.append(to_hex(plt.get_cmap("hsv")(((i - 20) * 0.61803398875 + 0.13) % 1)))
    return dict(zip(labels, colors))


def p_text(value):
    try:
        value = float(value)
        return f"{value:.3g}" if np.isfinite(value) else "not available"
    except (ValueError, TypeError):
        return "not available"


def motif_labels(motifs):
    return [str(v) for v in motifs.motif_id]


def empty_axis(ax, message):
    ax.text(0.5, 0.5, textwrap.fill(message, 48), ha="center", va="center", transform=ax.transAxes, color="#677785")
    ax.set_xticks([])
    ax.set_yticks([])


def legend_height(labels, cols=3):
    # Full labels are wrapped, never truncated. Reserve enough figure height.
    entries = [textwrap.wrap(str(v), width=34) or [""] for v in labels]
    rows = math.ceil(len(entries) / cols)
    return max(0.55, rows * 0.24 * max((len(v) for v in entries), default=1) + 0.25)


def type_legend(fig, palette, cols=3):
    if not palette:
        return
    handles = [Patch(facecolor=color, label=textwrap.fill(label, width=34)) for label, color in palette.items()]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.005), ncol=cols,
               fontsize=8.5, columnspacing=2.0, handlelength=1.1, handletextpad=0.5, title="Cell type")


def save(fig, name, output, manifest, description, dpi):
    paths = []
    for extension in ("png", "pdf"):
        path = output / f"{name}.{extension}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight", pad_inches=0.22)
        paths.append(str(path.resolve()))
        print(path)
    manifest.append({"name": name, "description": description, "files": paths})
    plt.close(fig)


def overview(run, tables, motifs, palette, output, manifest, dpi):
    legend_h = legend_height(palette)
    plot_h = max(8.5, 0.48 * len(motifs) + 4.5)
    fig, axs = plt.subplots(2, 2, figsize=(15, plot_h + legend_h))
    mode = run.get("mode", "discovery")
    fig.suptitle("Mixed-cell motifs | recurrence, composition and statistical evidence", fontsize=17, fontweight="bold", y=0.99)
    subtitle = f"{mode.capitalize()} · {run.get('k', '?')} cells/group · at least {run.get('min_types', '?')} types · dominance ≤ {run.get('max_dominance', '?')} · {run.get('permutations', len(tables['null']))} null replicates"
    fig.text(0.5, 0.955, subtitle, ha="center", fontsize=10, color="#586A77")
    ids = motif_labels(motifs)
    y = np.arange(len(motifs))
    ax = axs[0, 0]
    ax.set_title("A  Representative cell counts", loc="left")
    if len(motifs):
        comp = tables["motif_composition"].pivot_table(index="motif_id", columns="cell_type", values="count", aggfunc="sum", fill_value=0).reindex(ids, fill_value=0)
        left = np.zeros(len(ids))
        for label, color in palette.items():
            vals = comp[label].to_numpy() if label in comp else np.zeros(len(ids))
            ax.barh(y, vals, left=left, color=color, height=0.72, edgecolor="white", linewidth=0.5)
            for j, val in enumerate(vals):
                if val >= 1:
                    ax.text(left[j] + val / 2, j, str(int(val)), ha="center", va="center", fontsize=8,
                            bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none", "pad": 0.6})
            left += vals
        ax.set_yticks(y, ids)
        ax.invert_yaxis()
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_xlabel("Cells in template")
    else:
        empty_axis(ax, "No reportable mixed motifs. Check neighborhood eligibility and the minimum support.")
    ax = axs[0, 1]
    ax.set_title("B  Recurrence after removing overlap", loc="left")
    if len(motifs):
        ax.barh(y - 0.18, motifs.raw_matches, height=0.34, color=GRAY, label="Matching neighborhoods")
        ax.barh(y + 0.18, motifs.support, height=0.34, color=BLUE, label="Selected cell-disjoint occurrences")
        ax.set_yticks(y, ids)
        ax.set_ylim(len(ids) + 0.4, -0.6)
        ax.set_xlabel("Count")
        ax.grid(axis="x")
        ax.legend(loc="lower right", fontsize=8)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    else:
        empty_axis(ax, "No reportable motifs")
    ax = axs[1, 0]
    ax.set_title("C  Same-template support against the null", loc="left")
    if len(motifs):
        means = motifs.expected_support.to_numpy(float)
        ax.barh(y - 0.18, means, height=0.34, color=GRAY, label="Mean null support, same template")
        ax.barh(y + 0.18, motifs.support, height=0.34, color=TEAL, label="Observed support")
        nulls = tables["motif_null"]
        for j, mid in enumerate(ids):
            vals = nulls.loc[nulls.motif_id == mid, "support"].to_numpy(float)
            vals = vals[np.isfinite(vals)]
            if len(vals):
                lo, hi = np.quantile(vals, [0.025, 0.975])
                ax.plot([lo, hi], [j - 0.18] * 2, color="#566979", lw=1.3)
        ax.set_yticks(y, ids)
        ax.set_ylim(len(ids) + 0.4, -0.6)
        ax.set_xlabel("Cell-disjoint support; lines: central 95% null range")
        ax.grid(axis="x")
        ax.legend(loc="lower right", fontsize=8)
    else:
        empty_axis(ax, "No reportable motifs")
    if mode == "discovery":
        ax.text(0.0, -0.24, "Post-selection descriptive comparison; this panel is not a motif-specific test.", transform=ax.transAxes, fontsize=8.5, color="#677785")
    ax = axs[1, 1]
    if mode == "validate":
        ax.set_title("D  Fixed-template validation", loc="left")
        if len(motifs):
            pvals = motifs.p_value_holm.to_numpy(float)
            valid = np.isfinite(pvals)
            if valid.any():
                colors = [TEAL if p <= float(run.get("alpha", 0.05)) else GRAY for p in pvals[valid]]
                ax.barh(y[valid], -np.log10(np.maximum(pvals[valid], np.finfo(float).tiny)), color=colors, height=0.65)
                alpha = float(run.get("alpha", 0.05))
                if 0 < alpha < 1:
                    ax.axvline(-np.log10(alpha), color=ORANGE, ls="--", label=f"α = {alpha:g}")
                    ax.legend(fontsize=8)
                ax.set_yticks(y, ids)
                ax.invert_yaxis()
                ax.set_xlabel("−log₁₀(Holm-adjusted p-value)")
                ax.grid(axis="x")
            else:
                empty_axis(ax, "No fixed-template p-values available")
        else:
            empty_axis(ax, "No templates to validate")
    else:
        ax.set_title("D  Complete-search null calibration", loc="left")
        vals = tables["null"].max_support.to_numpy(float)
        vals = vals[np.isfinite(vals)]
        if len(vals):
            ax.hist(vals, bins=integer_bins(vals), color=GRAY, edgecolor="white")
            observed = run.get("observed_max_support", float(motifs.support.max()) if len(motifs) else 0)
            if observed is not None:
                ax.axvline(float(observed), color=ORANGE, lw=2, label=f"Observed maximum: {observed:g}")
            ax.set_xlabel("Maximum support after repeating the search")
            ax.set_ylabel("Null replicates")
            ax.text(0.98, 0.96, f"Global-null p = {p_text(run.get('global_p_value'))}",
                    transform=ax.transAxes, ha="right", va="top", fontweight="bold",
                    bbox={"facecolor": "white", "edgecolor": "none", "pad": 3}, zorder=10)
            ax.legend(loc="upper left", bbox_to_anchor=(0, 0.85), fontsize=8)
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        else:
            empty_axis(ax, "No null replicates available; statistical evidence cannot be displayed")
    type_legend(fig, palette)
    fig.subplots_adjust(left=0.075, right=0.98, top=0.90, bottom=(legend_h + 0.45) / (plot_h + legend_h), hspace=0.52, wspace=0.30)
    save(fig, "overview", output, manifest,
         "Template compositions, raw and cell-disjoint recurrence, same-template null support, and selection-aware global-null evidence (or fixed-template validation).", dpi)


def integer_bins(values):
    lo, hi = int(np.floor(np.min(values))), int(np.ceil(np.max(values)))
    return np.arange(lo - 0.5, hi + 1.5) if hi - lo <= 60 else 30


def diagnostics(run, tables, output, manifest, dpi):
    nh = tables["neighborhoods"]
    fig, axs = plt.subplots(2, 2, figsize=(12, 8.7))
    fig.suptitle("Neighborhood diagnostics | geometric scale and mixed composition", fontsize=16, fontweight="bold", y=0.98)
    fig.text(0.5, 0.935, f"{len(nh):,} geometric neighborhoods · {int(nh.eligible.sum()):,} eligible mixed groups · coordinate distances use input units", ha="center", color="#586A77")
    for ax in axs.flat:
        ax.grid(axis="y")
    ax = axs[0, 0]
    ax.set_title("A  Physical neighborhood radius", loc="left")
    vals = nh.radius.to_numpy(float)
    finite = np.isfinite(vals)
    if finite.any():
        edges = np.histogram_bin_edges(vals[finite], bins=35)
        ax.hist(vals[finite], bins=edges, color=GRAY, label="All geometric neighborhoods")
        ax.hist(nh.loc[nh.eligible, "radius"], bins=edges, color=TEAL, alpha=0.9, label="Eligible mixed groups")
        ax.set_xlabel("Distance from center to farthest member")
        ax.set_ylabel("Neighborhoods")
        ax.legend(fontsize=8)
    else:
        empty_axis(ax, "No neighborhoods passed the geometric constraints")
    ax = axs[0, 1]
    ax.set_title("B  Cell-type diversity", loc="left")
    if len(nh):
        bins = integer_bins(nh.n_types.to_numpy(float))
        ax.hist(nh.n_types, bins=bins, color=GRAY, label="All")
        ax.hist(nh.loc[nh.eligible, "n_types"], bins=bins, color=TEAL, label="Eligible")
        ax.axvline(float(run.get("min_types", 2)) - 0.5, color=ORANGE, ls="--", lw=1)
        ax.set_xlabel("Number of cell types in group")
        ax.set_ylabel("Neighborhoods")
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    else:
        empty_axis(ax, "No geometric neighborhoods")
    ax = axs[1, 0]
    ax.set_title("C  Dominant-type fraction", loc="left")
    if len(nh):
        ax.hist(nh.dominance, bins=np.linspace(0, 1.00001, 21), color=GRAY)
        ax.hist(nh.loc[nh.eligible, "dominance"], bins=np.linspace(0, 1.00001, 21), color=TEAL)
        cutoff = float(run.get("max_dominance", 0.7))
        ax.axvline(cutoff, color=ORANGE, ls="--", label=f"Maximum allowed: {cutoff:g}")
        ax.legend(fontsize=8)
        ax.set_xlabel("Fraction of cells assigned to the most common type")
        ax.set_ylabel("Neighborhoods")
        ax.set_xlim(0, 1.03)
    else:
        empty_axis(ax, "No geometric neighborhoods")
    ax = axs[1, 1]
    ax.set_title("D  Diversity and local mixing", loc="left")
    good = np.isfinite(nh.effective_types) & np.isfinite(nh.cross_edge_fraction)
    contact_radius = float(run.get("contact_radius", 0) or 0)
    if contact_radius <= 0:
        empty_axis(ax, "Local edge-mixing diagnostic disabled. Set a positive contact radius during analysis to summarize nearby cross-type cell pairs.")
    elif good.any():
        plot = nh.loc[good]
        ax.hexbin(plot.effective_types, plot.cross_edge_fraction, gridsize=25, mincnt=1,
                  cmap="Blues", linewidths=0.1, bins="log")
        ax.set_xlabel("Effective number of cell types, exp(entropy)")
        ax.set_ylabel("Cross-type fraction of pairs within contact radius")
        ax.set_ylim(-0.03, 1.03)
        ax.text(0.03, 0.97, f"Contact radius: {contact_radius:g}; darker bins = more groups", transform=ax.transAxes, va="top", fontsize=8.5)
    else:
        empty_axis(ax, "No edge-mixing diagnostics available")
    fig.subplots_adjust(left=0.08, right=0.98, top=0.86, bottom=0.09, hspace=0.38, wspace=0.28)
    save(fig, "neighborhood_diagnostics", output, manifest,
         "Input-coordinate radii, type counts, dominance, and effective diversity versus cross-type edge fraction. Gray includes all geometric groups; teal shows eligible groups.", dpi)


def sample_maps(run, tables, motifs, palette, output, manifest, dpi):
    cells, occurrences, members = tables["cells"], tables["occurrences"], tables["members"]
    if cells.empty:
        return
    ids = motif_labels(motifs)
    for number, sample in enumerate(sorted(cells["sample"].unique()), 1):
        part = cells[cells["sample"] == sample]
        cols = min(3, len(ids) + 1)
        rows = math.ceil((len(ids) + 1) / cols)
        legend_h = legend_height(palette)
        xspan = float(part.x.max() - part.x.min())
        yspan = float(part.y.max() - part.y.min())
        spatial_ratio = yspan / max(xspan, 1e-9)
        row_height = min(6.2, max(2.4, 4.3 * spatial_ratio + 0.6))
        total_h = rows * row_height + 1.0 + legend_h
        fig, axs = plt.subplots(rows, cols, figsize=(max(9, 5.0 * cols), total_h), squeeze=False)
        fig.suptitle(textwrap.fill(f"Spatial recurrence | sample {sample}", 100), fontsize=16, fontweight="bold", y=0.985)
        fig.text(0.5, 0.942, "Each motif panel highlights actual member cells in its selected cell-disjoint occurrences.", ha="center", fontsize=9.5, color="#586A77")
        ax = axs.flat[0]
        # Dense scatter is rasterized inside the PDF; labels and axes remain vector.
        ax.scatter(part.x, part.y, c=part.cell_type.map(palette), s=4, linewidths=0, rasterized=True)
        ax.set_title(f"All cells ({len(part):,})", loc="left")
        for slot, mid in enumerate(ids, 1):
            ax = axs.flat[slot]
            ax.scatter(part.x, part.y, c="#D9DEE3", s=3, linewidths=0, rasterized=True)
            selected = occurrences[(occurrences.motif_id == mid) & (occurrences["sample"] == sample)]
            chosen = members[(members.motif_id == mid) & members.occurrence_id.isin(selected.occurrence_id)]
            selected_cells = part[part.cell_index.isin(chosen.cell_index)]
            if len(selected_cells):
                ax.scatter(selected_cells.x, selected_cells.y, c=selected_cells.cell_type.map(palette), s=10, linewidths=0, rasterized=True)
                # Outline centers makes occurrence counts visible without treating radii as boundaries.
                ax.scatter(selected.x, selected.y, facecolors="none", edgecolors=INK, s=19, linewidths=0.45, rasterized=True)
            else:
                ax.text(0.5, 0.5, "No selected occurrences in this sample", transform=ax.transAxes, ha="center", fontsize=9)
            ax.set_title(f"{mid} · {len(selected):,} occurrences", loc="left")
        margin = max(xspan, yspan, 1) * 0.025
        for ax in list(axs.flat)[:len(ids) + 1]:
            ax.set_xlim(part.x.min() - margin, part.x.max() + margin)
            ax.set_ylim(part.y.min() - margin, part.y.max() + margin)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlabel("x (input units)")
            ax.set_ylabel("y (input units)")
            ax.ticklabel_format(style="sci", axis="both", scilimits=(-3, 4))
        for ax in list(axs.flat)[len(ids) + 1:]:
            ax.set_visible(False)
        type_legend(fig, palette)
        fig.subplots_adjust(left=0.07, right=0.98, top=0.89, bottom=(legend_h + 0.45) / total_h, hspace=0.38, wspace=0.28)
        slug = re.sub(r"[^A-Za-z0-9_-]+", "_", str(sample)).strip("_")[:70] or "unnamed"
        save(fig, f"spatial_{number:02d}_{slug}", output, manifest,
             f"Sample {sample}: cell-type reference and selected occurrences of each displayed motif; different specimens never share a spatial coordinate panel.", dpi)


def representative_occurrences(occurrences, count):
    """Select deterministically: cover specimens, then favor separated centers."""
    if len(occurrences) <= count:
        return occurrences
    groups = [g.sort_values(["x", "y", "occurrence_id"]).copy() for _, g in occurrences.groupby("sample", sort=True)]
    selected_indices = []
    while len(selected_indices) < count and groups:
        remaining = []
        for group in groups:
            if len(selected_indices) == count:
                break
            previous = occurrences.loc[[idx for idx in selected_indices if occurrences.loc[idx, "sample"] == group.iloc[0]["sample"]]]
            if previous.empty:
                chosen = group.index[len(group) // 2]
            else:
                xy = group[["x", "y"]].to_numpy(float)
                old = previous[["x", "y"]].to_numpy(float)
                distance = ((xy[:, None, :] - old[None, :, :]) ** 2).sum(axis=2).min(axis=1)
                chosen = group.index[int(np.argmax(distance))]
            selected_indices.append(chosen)
            group = group.drop(chosen)
            if len(group):
                remaining.append(group)
        groups = remaining
    return occurrences.loc[selected_indices]


def occurrence_atlas(run, tables, motifs, palette, output, manifest, dpi, atlas_matches):
    if motifs.empty:
        return
    cells = tables["cells"].set_index("cell_index", drop=False)
    occurrences, members = tables["occurrences"], tables["members"]
    selected = {mid: representative_occurrences(occurrences[occurrences.motif_id == mid], atlas_matches) for mid in motif_labels(motifs)}
    # Paginate to keep PDF pages readable, including for many requested motifs.
    ids = motif_labels(motifs)
    for page in range(math.ceil(len(ids) / 4)):
        page_ids = ids[page * 4:(page + 1) * 4]
        legend_h = legend_height(palette)
        total_h = 3.8 * len(page_ids) + 1.4 + legend_h
        fig, axs = plt.subplots(len(page_ids), atlas_matches, figsize=(max(9, 4.6 * atlas_matches), total_h), squeeze=False)
        fig.suptitle("Occurrence atlas | actual cells in recurring groups", fontsize=16, fontweight="bold", y=0.985)
        fig.text(0.5, 0.95, "Groups are centered for comparison; distance and aspect ratio are preserved. Ring: neighborhood radius. +: center cell.", ha="center", fontsize=9, color="#586A77")
        radii = [float(row.radius) for mid in page_ids for row in selected[mid].itertuples() if np.isfinite(row.radius)]
        # Include the true member extents as well as the exported radius. This
        # keeps every cell visible even if metadata was rounded or edited.
        for mid in page_ids:
            for occurrence in selected[mid].itertuples():
                indices = members.loc[(members.motif_id == mid) & (members.occurrence_id == occurrence.occurrence_id), "cell_index"]
                group = cells.loc[indices]
                if len(group):
                    extents = np.hypot(group.x.to_numpy(float) - occurrence.x, group.y.to_numpy(float) - occurrence.y)
                    radii.append(float(extents.max()))
        half = max(radii, default=1) * 1.17
        half = max(half, 0.1)
        for rowno, mid in enumerate(page_ids):
            matches = selected[mid]
            for col in range(atlas_matches):
                ax = axs[rowno, col]
                if col >= len(matches):
                    empty_axis(ax, f"{mid}: fewer than {col + 1} selected occurrences")
                    continue
                occurrence = matches.iloc[col]
                indices = members.loc[(members.motif_id == mid) & (members.occurrence_id == occurrence.occurrence_id), "cell_index"]
                group = cells.loc[indices]
                if group.empty:
                    empty_axis(ax, f"{mid}: missing occurrence members")
                    continue
                dx = group.x.to_numpy(float) - float(occurrence.x)
                dy = group.y.to_numpy(float) - float(occurrence.y)
                ax.add_patch(Circle((0, 0), float(occurrence.radius), fill=False, edgecolor="#B7C0C8", ls="--", lw=0.8))
                ax.scatter(dx, dy, c=group.cell_type.map(palette), s=115, edgecolors="white", linewidths=0.7, zorder=3)
                ax.scatter([0], [0], marker="+", c=INK, s=32, linewidths=1.0, zorder=4)
                title = f"{mid} · sample {occurrence['sample']}\nOccurrence {occurrence.occurrence_id} · {len(group)} cells / {group.cell_type.nunique()} types"
                ax.set_title(textwrap.fill(title.split('\n')[0], 43) + '\n' + title.split('\n')[1], fontsize=10, loc="left")
                ax.set_xlim(-half, half)
                ax.set_ylim(-half, half)
                ax.set_aspect("equal", adjustable="box")
                ax.set_xlabel("Δx (input units)")
                ax.set_ylabel("Δy (input units)")
                ax.ticklabel_format(style="sci", axis="both", scilimits=(-3, 4))
        type_legend(fig, palette)
        fig.subplots_adjust(left=0.07, right=0.98, top=0.88, bottom=(legend_h + 0.4) / total_h, hspace=0.48, wspace=0.3)
        save(fig, f"occurrence_atlas_{page + 1:02d}", output, manifest,
             "Actual occurrence member cells, selected to span samples and spatial positions. All panels on a page use the same spatial scale; '+' marks the center cell.", dpi)


def sample_recurrence(tables, motifs, output, manifest, dpi):
    cells, occurrences = tables["cells"], tables["occurrences"]
    samples = sorted(cells["sample"].unique())
    if motifs.empty or not samples:
        return
    ids = motif_labels(motifs)
    pivot = occurrences[occurrences.motif_id.isin(ids)].groupby(["motif_id", "sample"]).size().unstack(fill_value=0).reindex(index=ids, columns=samples, fill_value=0).fillna(0)
    # Counts and exposure-normalized descriptive rates are both provided.
    exposures = cells.groupby("sample").size().reindex(samples).to_numpy(float)
    rates = pivot.to_numpy(float) / exposures[None, :] * 1000
    width = max(11, 1.0 * len(samples) + 7)
    height = max(4.5, 0.48 * len(ids) + 2.3)
    fig, axs = plt.subplots(1, 2, figsize=(width, height))
    fig.suptitle("Recurrence across samples", fontsize=16, fontweight="bold", y=0.98)
    for ax, values, title, fmt in zip(axs, [pivot.to_numpy(), rates], ["Selected cell-disjoint occurrences", "Occurrences per 1,000 input cells (descriptive)"], [".0f", ".1f"]):
        mesh = ax.imshow(values, aspect="auto", cmap="Blues", vmin=0)
        ax.set_title(title, loc="left", fontsize=11)
        ax.set_xticks(np.arange(len(samples)), [textwrap.fill(s, 18) for s in samples], rotation=45 if len(samples) > 4 else 0, ha="right" if len(samples) > 4 else "center")
        ax.set_yticks(np.arange(len(ids)), ids)
        maximum = max(float(np.max(values)), 1)
        for i in range(len(ids)):
            for j in range(len(samples)):
                ax.text(j, i, format(values[i, j], fmt), ha="center", va="center", fontsize=9,
                        color="white" if values[i, j] > maximum * 0.6 else INK)
        fig.colorbar(mesh, ax=ax, shrink=0.78, pad=0.025)
    fig.subplots_adjust(top=0.82, bottom=0.26 if len(samples) > 4 else 0.17, left=0.08, right=0.98, wspace=0.30)
    save(fig, "sample_recurrence", output, manifest,
         "Occurrence counts and descriptive rates per 1,000 input cells by sample. Rates do not adjust for cell-type abundance or tissue compartments.", dpi)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", type=Path, required=True, help="Directory containing run.json and result CSV tables")
    parser.add_argument("--output", type=Path, help="Figure directory (default: RESULTS/plots)")
    parser.add_argument("--top", type=int, default=8, help="Number of motifs to display in result-table order (default: 8)")
    parser.add_argument("--atlas-matches", type=int, default=3, help="Example occurrences per motif (default: 3; maximum: 6)")
    parser.add_argument("--dpi", type=int, default=180, help="PNG resolution (default: 180)")
    args = parser.parse_args()
    if args.top < 1 or not 1 <= args.atlas_matches <= 6 or args.dpi < 50:
        parser.error("--top must be positive, --atlas-matches must be 1–6, and --dpi must be at least 50")
    results = args.results.expanduser().resolve()
    output = (args.output or results / "plots").expanduser().resolve()
    try:
        run, tables = load_tables(results)
    except ValueError as exc:
        parser.error(str(exc))
    output.mkdir(parents=True, exist_ok=True)
    setup_style()
    palette = palette_for(tables["cells"].cell_type)
    # Preserve C++ reporting order rather than reranking selected results in Python.
    motifs = tables["motifs"].head(args.top)
    manifest = []
    overview(run, tables, motifs, palette, output, manifest, args.dpi)
    diagnostics(run, tables, output, manifest, args.dpi)
    sample_recurrence(tables, motifs, output, manifest, args.dpi)
    sample_maps(run, tables, motifs, palette, output, manifest, args.dpi)
    occurrence_atlas(run, tables, motifs, palette, output, manifest, args.dpi, args.atlas_matches)
    pd.DataFrame({"cell_type": list(palette), "color_hex": list(palette.values())}).to_csv(output / "cell_type_palette.csv", index=False)
    notes = [
        "Cell-disjoint selections can still be spatially or biologically dependent; they are not independent biological replicates.",
        "Coordinate distances are in the original input units. Cell types use one shared palette across all figures.",
        "Different motifs may share member cells. Spatial disjointness is enforced within each motif's reported support.",
        "Cross-type edges describe local mixing, not a test of a higher-order interaction.",
    ]
    if run.get("mode", "discovery") == "discovery":
        notes += ["Discovery p-values compare with maxima from complete null searches and describe evidence against the specified global null.",
                  "Same-template null means, intervals, and enrichment are descriptive after selection; they are not motif-specific confirmatory tests."]
    else:
        notes += ["Fixed-template p-values are meaningful when the templates and analysis settings were fixed independently of the validation tissue."]
    if run.get("null_assumption"):
        notes.append(str(run["null_assumption"]))
    record = {"results": str(results), "mode": run.get("mode"), "displayed_motif_ids": motif_labels(motifs), "figures": manifest, "notes": notes}
    (output / "plots_manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    captions = ["# Mixed-cell motif figures", "", *[f"- {note}" for note in notes], ""]
    for item in manifest:
        captions += [f"## {item['name']}", "", item["description"], ""]
    (output / "figure_captions.md").write_text("\n".join(captions))
    print(f"Wrote {len(manifest)} figure pairs and plots_manifest.json to {output}")


if __name__ == "__main__":
    main()
