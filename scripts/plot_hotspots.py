"""Where the significant patterns are, all at once.

    python scripts/plot_hotspots.py RUN_DIR

Writes RUN_DIR/plots/hotspots.png, one row per sample:

* left, how many distinct significant patterns have an occurrence centred in
  each hexagon: the places where the tissue's structure departs most from
  the shuffles;
* right, every occurrence centre coloured by the main cell type of the best
  pattern centred there, so the hotspots can be told apart by kind.

occurrences.csv lists the occurrences of significant patterns only. The y axis
points down, as in the images the coordinates come from.
"""

from __future__ import annotations

import argparse
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

from pam_results import (DARK, FAINT, MUTED, label_colours, load_cells, load_run, output_path,
                         parse_composition, short_names, use_style)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("run", help="a pam_st --output-dir")
    p.add_argument("--input", action="append",
                   help="input table, repeated in pam_st's order, if summary.csv's paths moved")
    p.add_argument("--hex", type=float,
                   help="hexagon width in data units (default: twice the radius r)")
    p.add_argument("--out", help="figure path (default: RUN_DIR/plots/hotspots.png)")
    p.add_argument("--dpi", type=int, default=200)
    return p.parse_args()


def tissue(ax, cells) -> None:
    ax.scatter(cells["x"], cells["y"], s=0.3, c=FAINT, linewidths=0, rasterized=True)
    ax.set_aspect("equal")
    ax.set_xlabel("x", color=MUTED)
    ax.set_ylabel("y", color=MUTED)


def main() -> None:
    args = parse_args()
    use_style()
    run = load_run(args.run)
    occ = run.occurrences
    if occ.empty:
        raise SystemExit("occurrences.csv is empty: no pattern is significant.")
    cells = load_cells(run, args.input)
    names = list(cells)
    short = short_names(names)

    # The main cell type of each significant pattern, and the best pattern
    # (lowest rank) centred at each spot.
    main_type = {row["rank"]: max(parse_composition(row["mean_composition"]).items(),
                                  key=lambda kv: kv[1])[0]
                 for _, row in run.patterns[run.patterns["significant"]].iterrows()}
    centres = occ.sort_values("rank").drop_duplicates(["sample", "spot"]).copy()
    centres["type"] = centres["rank"].map(main_type)
    types = centres["type"].value_counts().index.tolist()
    colours = label_colours(types)

    hex_width = args.hex or 2 * run.number("radius")
    fig, axes = plt.subplots(len(names), 2, figsize=(14, 6.6 * len(names)), squeeze=False,
                             layout="constrained")
    for (left, right), name in zip(axes, names):
        here = occ[occ["sample"] == name]
        c = cells[name]
        tissue(left, c)
        tissue(right, c)
        if not here.empty:
            span = max(np.ptp(c["x"]), hex_width)
            grid = int(np.clip(math.ceil(span / hex_width), 10, 200))
            hexes = left.hexbin(here["x"], here["y"], C=here["rank"], gridsize=grid,
                                reduce_C_function=lambda ranks: len(np.unique(ranks)),
                                cmap="Blues", mincnt=1, linewidths=0.2, edgecolors="white",
                                extent=(c["x"].min(), c["x"].max(), c["y"].min(), c["y"].max()))
            bar = fig.colorbar(hexes, ax=left, shrink=0.6, label="distinct significant patterns")
            bar.locator = MaxNLocator(integer=True)
            bar.update_ticks()
            spots = centres[centres["sample"] == name]
            for kind in reversed(types):
                part = spots[spots["type"] == kind]
                right.scatter(part["x"], part["y"], s=4, c=colours[kind], linewidths=0,
                              rasterized=True)
        for ax in (left, right):
            ax.invert_yaxis()
        n_patterns = here["rank"].nunique()
        left.set_title(f"{short[name]}: {n_patterns} significant patterns occur here",
                       fontweight="normal")
        right.set_title(f"{short[name]}: {here['spot'].nunique()} centres, by main cell type",
                        fontweight="normal")

    dot = dict(marker="o", linestyle="", markersize=6)
    handles = [Line2D([], [], color=colours[t], label=t, **dot) for t in types if colours[t] != DARK]
    if DARK in colours.values():
        handles.append(Line2D([], [], color=DARK, label="other main types", **dot))
    handles.append(Line2D([], [], color=FAINT, label="cells", **dot))
    fig.legend(handles=handles, loc="outside lower right", ncol=len(handles), fontsize=9)
    fig.suptitle(f"Where the significant patterns are: {run.name}\n{run.parameters()}   ·   "
                 f"hexagons {hex_width:g} wide", x=0.01, ha="left", fontsize=12,
                 fontweight="bold")
    out = output_path(run, "hotspots.png", args.out)
    fig.savefig(out, dpi=args.dpi)
    print(f"{occ['rank'].nunique()} significant patterns, {len(centres)} centres -> {out}")


if __name__ == "__main__":
    main()
