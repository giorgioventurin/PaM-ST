"""Where the time goes, for one run or many.

    python scripts/plot_runtime.py RUN_DIR
    python scripts/plot_runtime.py output/            # every run below it
    python scripts/plot_runtime.py RUN_A RUN_B --out runtime.png

Writes a figure, and a CSV with the same numbers beside it:

* each run's time split into its phases: loading the input, setting up
  (neighbourhoods, tiles, p-values), the three passes over the m+1 labellings
  (the universe of patterns, their recurrences, the scores) and the output;
* seconds per labelling against the number of cells and against the size of
  the universe, which the recurrence pass scales with.

With one run the figure goes to RUN_DIR/plots/runtime.png; with several, to
runtime.png in the current directory unless --out says otherwise.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter, LogLocator

from pam_results import DARK, INK, MUTED, PHASES, SERIES, find_runs, load_run, use_style

PHASE_NAMES = {
    "load_seconds": "load input",
    "setup_seconds": "set-up",
    "universe_seconds": "pass 1: universe",
    "recurrence_seconds": "pass 2: recurrences",
    "score_seconds": "pass 3: scores",
    "output_seconds": "write output",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("paths", nargs="+", help="pam_st output directories, or folders holding them")
    p.add_argument("--out", help="figure path")
    p.add_argument("--dpi", type=int, default=200)
    return p.parse_args()


def gather(paths: list[Path]) -> pd.DataFrame:
    rows = []
    for path in paths:
        run = load_run(path)
        n = run.number
        row = {
            "run": run.name,
            "cells": n("cells"),
            "samples": n("samples", 1),
            "labels": n("labels"),
            "radius": n("radius"),
            "rho": n("rho"),
            "min_recurrence": n("min_recurrence"),
            "permutations": n("permutations"),
            "threads": n("threads"),
            "observed_patterns": n("observed_patterns"),
            "universe_size": n("universe_size"),
            "significant_patterns": n("significant_patterns"),
            **{key: n(key) for key in ["load_seconds", "analysis_seconds", "universe_seconds",
                                       "recurrence_seconds", "score_seconds", "output_seconds",
                                       "total_seconds"]},
        }
        passes = row["universe_seconds"] + row["recurrence_seconds"] + row["score_seconds"]
        # Runs from before the pass timings count all their analysis as set-up.
        row["setup_seconds"] = (row["analysis_seconds"] - passes if np.isfinite(passes)
                                else row["analysis_seconds"])
        row["seconds_per_labelling"] = row["analysis_seconds"] / (row["permutations"] + 1)
        rows.append(row)
    return pd.DataFrame(rows).sort_values("total_seconds", ascending=False, ignore_index=True)


def readable(seconds: float) -> str:
    if seconds < 120:
        return f"{seconds:.1f} s"
    if seconds < 7200:
        return f"{seconds / 60:.1f} min"
    return f"{seconds / 3600:.1f} h"


def phases(ax, table: pd.DataFrame) -> None:
    y = np.arange(len(table))[::-1]
    left = np.zeros(len(table))
    for colour, (key, name) in zip(PHASES, PHASE_NAMES.items()):
        width = table[key].fillna(0).to_numpy()
        ax.barh(y, width, left=left, color=colour, edgecolor="white", linewidth=1, height=0.6,
                label=name)
        left += width
    for yi, total in zip(y, table["total_seconds"]):
        ax.text(total, yi, f"  {readable(total)}", va="center", fontsize=8, color=INK)
    ax.set_yticks(y, table["run"])
    ax.set_xlabel("seconds")
    ax.set_xlim(0, table["total_seconds"].max() * 1.18)
    ax.set_title("Time by phase")
    ax.legend(loc="lower right", fontsize=8, ncol=1)


def scaling(ax, table: pd.DataFrame, column: str, label: str) -> None:
    """Seconds per labelling against `column`, coloured by rho: a tolerance
    above zero costs far more than exact matching."""
    rhos = sorted(table["rho"].unique())
    for i, rho in enumerate(rhos):
        part = table[table["rho"] == rho]
        ax.scatter(part[column], part["seconds_per_labelling"], s=36, zorder=3,
                   c=SERIES[i] if i < len(SERIES) else DARK, edgecolors="white", linewidths=1,
                   label=f"ρ = {rho:g}")
    for _, row in table.iterrows():
        ax.annotate(row["run"], (row[column], row["seconds_per_labelling"]),
                    xytext=(5, 3), textcoords="offset points", fontsize=7, color=MUTED)
    ax.set_xscale("log")
    ax.set_yscale("log")
    plain = FuncFormatter(lambda v, _: f"{v:g}")
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_minor_locator(LogLocator(base=10, subs=(2, 5)))
        axis.set_major_formatter(plain)
        axis.set_minor_formatter(plain)
    ax.tick_params(which="minor", labelsize=7)
    ax.set_xlabel(label)
    ax.set_ylabel("seconds per labelling")
    ax.set_title(f"Time per labelling against {label.split(' |')[0]}")
    ax.grid(True, which="major", lw=0.5)
    ax.legend(loc="upper left", fontsize=8)


def main() -> None:
    args = parse_args()
    use_style()
    paths = find_runs(args.paths)
    if not paths:
        raise SystemExit("No pam_st output directory (with a summary.csv) found.")
    table = gather(paths)

    if args.out:
        out = Path(args.out).expanduser()
    elif len(paths) == 1:
        out = paths[0] / "plots" / "runtime.png"
    else:
        out = Path.cwd() / "runtime.png"
    out.parent.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(15, max(4.8, 1.6 + 0.42 * len(table))), layout="constrained")
    grid = fig.add_gridspec(1, 3, width_ratios=[1.5, 1, 1])
    phases(fig.add_subplot(grid[0]), table)
    scaling(fig.add_subplot(grid[1]), table, "cells", "cells")
    scaling(fig.add_subplot(grid[2]), table, "universe_size", "universe size |U|")

    total = table["total_seconds"].sum()
    share = table["recurrence_seconds"].sum() / table["analysis_seconds"].sum()
    note = f", {share:.0%} of the analysis in pass 2" if np.isfinite(share) else ""
    fig.suptitle(f"PaM-ST running time: {len(table)} run{'s' if len(table) > 1 else ''}, "
                 f"{readable(total)} in all{note}", x=0.01, ha="left", fontsize=12,
                 fontweight="bold")
    fig.savefig(out, dpi=args.dpi)
    table.to_csv(out.with_suffix(".csv"), index=False)
    print(f"{len(table)} runs -> {out} and {out.with_suffix('.csv').name}")


if __name__ == "__main__":
    main()
