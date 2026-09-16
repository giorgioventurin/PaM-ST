#!/usr/bin/env python3
"""Compare completed discovery runs with identical analysis settings.

Writes CSV/JSON and PNG/PDF null-distribution panels. Holm adjustment concerns
the supplied dataset-level global tests, never individual motif discoveries.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "mixed_motifs_mpl"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

SETTINGS = ["k", "max_radius", "min_types", "max_dominance", "max_relabels",
            "same_types", "min_support", "max_motifs", "candidate_pool",
            "null_model", "block_size", "region_size", "contact_radius",
            "min_cross_edge_fraction", "permutations", "alpha", "seed"]


def rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def compare(paths, output):
    loaded = []
    for path in paths:
        run = json.loads((path / "run.json").read_text())
        if run["mode"] != "discovery":
            raise ValueError(f"{path}: comparison requires discovery runs")
        if loaded:
            differences = [key for key in SETTINGS if run.get(key) != loaded[0][1].get(key)]
            if differences:
                raise ValueError(f"{path}: settings differ: {', '.join(differences)}")
        null = [int(row["max_support"]) for row in rows(path / "null.csv")]
        if len(null) != run["permutations"]:
            raise ValueError(f"{path}: incomplete permutation output")
        computed_p = (1 + sum(v >= run["observed_max_support"] for v in null)) / (len(null) + 1)
        if not math.isclose(computed_p, run["global_p_value"], rel_tol=1e-9):
            raise ValueError(f"{path}: stored p-value does not match null maxima")
        loaded.append((path, run, null, rows(path / "motifs.csv")))
    adjusted = [1.0] * len(loaded)
    order = sorted(range(len(loaded)), key=lambda i: loaded[i][1]["global_p_value"])
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (len(order) - rank) * loaded[i][1]["global_p_value"]))
        adjusted[i] = running
    summary = []
    for i, (path, run, null, motifs) in enumerate(loaded):
        top = motifs[0] if motifs else {}
        summary.append({
            "dataset": Path(run["input"]).name.split("_")[0],
            "classified_cells": run["cells"],
            "distinct_neighborhoods": run["n_neighborhoods"],
            "eligible_mixed_neighborhoods": run["n_eligible"],
            "eligible_fraction": run["n_eligible"] / run["n_neighborhoods"] if run["n_neighborhoods"] else 0,
            "recurring_candidate_templates": run["n_candidates"],
            "reported_templates": run["n_reported_motifs"],
            "maximum_disjoint_support": run["observed_max_support"],
            "mean_null_maximum": sum(null) / len(null),
            "global_p_value": run["global_p_value"],
            "global_p_holm_across_datasets": adjusted[i],
            "dataset_significant_after_holm": adjusted[i] <= run["alpha"],
            "reported_search_exceedances": sum(m["significant"].lower() == "true" for m in motifs),
            "top_descriptive_template": top.get("composition", ""),
            "top_descriptive_support": top.get("support", 0),
            "top_descriptive_expected_support": top.get("expected_support", 0),
            "elapsed_seconds": run["elapsed_seconds"],
            "results_directory": str(path),
        })
    output.mkdir(parents=True, exist_ok=True)
    with (output / "dataset_summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    record = {
        "settings": {key: loaded[0][1].get(key) for key in SETTINGS},
        "dataset_count": len(loaded),
        "holm_family": [row["dataset"] for row in summary],
        "notes": [
            "Each panel compares the observed maximum disjoint support with a repeated complete-search null.",
            "Holm adjustment covers the dataset-level tests listed here; it is not motif-level FDR control.",
            "Different files do not necessarily represent independent patients.",
            "Raw support is not directly comparable across tissues with different cell counts or composition.",
            "Null shuffles preserve tile composition, not within-tile homotypic clustering.",
            "The top descriptive template is ranked by excess support and need not maximize raw support.",
        ],
        "results": summary,
    }
    (output / "dataset_summary.json").write_text(json.dumps(record, indent=2) + "\n")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42})
    ncols = min(2, len(loaded))
    nrows = math.ceil(len(loaded) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 4.2 * nrows), squeeze=False)
    for ax, item, row in zip(axes.flat, loaded, summary):
        _, run, null, _ = item
        bins = min(30, max(1, len(set(null))))
        ax.hist(null, bins=bins, color="#BAC6D1", edgecolor="white", label="Maximum from each null search")
        ax.axvline(run["observed_max_support"], color="#D67A26", linewidth=2,
                   label=f"Observed maximum: {run['observed_max_support']}")
        ax.set_title(f"{row['dataset']}  |  {run['cells']:,} classified cells", loc="left", fontweight="bold")
        ax.text(.97, .95, f"Global p = {run['global_p_value']:.3g}\nHolm across {len(loaded)} datasets = {row['global_p_holm_across_datasets']:.3g}",
                transform=ax.transAxes, ha="right", va="top", fontsize=10,
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 4}, zorder=10)
        ax.set_xlabel("Maximum cell-disjoint support")
        ax.set_ylabel("Null replicates")
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.legend(loc="lower center", bbox_to_anchor=(.5, -.37), frameon=False, fontsize=9)
    for ax in list(axes.flat)[len(loaded):]:
        ax.set_visible(False)
    cfg = record["settings"]
    null_label = "blockwise label shuffles" if cfg["null_model"] == "blocks" else f"{cfg['null_model']} permutations"
    fig.suptitle("Mixed cell-type patterns | comparison across tissues", fontsize=17, fontweight="bold", y=.98)
    fig.text(.5, .93, f"Identical settings: {cfg['k']} cells/group, at least {cfg['min_types']} types, dominance ≤ {cfg['max_dominance']}; {cfg['permutations']} {null_label}", ha="center", fontsize=11)
    fig.subplots_adjust(top=.85, bottom=.12, hspace=.7, wspace=.22)
    for suffix in ("png", "pdf"):
        fig.savefig(output / f"null_comparison.{suffix}", dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    for row in summary:
        print(f"{row['dataset']}: {row['eligible_mixed_neighborhoods']:,} eligible neighborhoods; maximum support {row['maximum_disjoint_support']}; p={row['global_p_value']:.3g}; Holm={row['global_p_holm_across_datasets']:.3g}")
    print(f"Comparison outputs: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", nargs="+", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    resolved = [p.expanduser().resolve() for p in args.results]
    if len(set(resolved)) != len(resolved):
        parser.error("Each result directory must be supplied once")
    try:
        compare(resolved, args.output.expanduser().resolve())
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, f"Comparison failed: {error}\n")
