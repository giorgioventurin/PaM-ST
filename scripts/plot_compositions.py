"""What the significant patterns are made of, side by side.

    python scripts/plot_compositions.py RUN_DIR
    python scripts/plot_compositions.py RUN_DIR --top 80 --all

Writes RUN_DIR/plots/compositions.png: one row per pattern, the highest
scores first, showing the mean share of each cell type in its matching
neighbourhoods. Rows are grouped by their most common cell type, so families
of near-identical patterns (the same clump at different sizes, say) sit
together. Beside the heatmap: each pattern's recurrence and z-score.
"""

from __future__ import annotations

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pam_results import (INK, MUTED, NOT_SIGNIFICANT, SIGNIFICANT, load_run, output_path,
                         parse_composition, use_style)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("run", help="a pam_st --output-dir")
    p.add_argument("--top", type=int, default=50, help="patterns shown (default: 50)")
    p.add_argument("--all", action="store_true",
                   help="draw from every observed pattern, not only significant ones")
    p.add_argument("--out", help="figure path (default: RUN_DIR/plots/compositions.png)")
    p.add_argument("--dpi", type=int, default=200)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    use_style()
    run = load_run(args.run)
    pool = run.patterns if args.all else run.patterns[run.patterns["significant"]]
    if pool.empty:
        raise SystemExit("No significant pattern to show; --all draws from every pattern.")
    rows = pool.sort_values("rank").head(args.top).reset_index(drop=True)

    shares = pd.DataFrame([parse_composition(t) for t in rows["mean_composition"]]).fillna(0.0)
    columns = shares.mean().sort_values(ascending=False).index.tolist()
    shares = shares[columns]
    dominant = shares.idxmax(axis=1)
    # Groups in the column order, each in rank order.
    order = sorted(range(len(rows)), key=lambda i: (columns.index(dominant[i]), rows["rank"][i]))
    rows, shares, dominant = rows.iloc[order], shares.iloc[order], dominant.iloc[order]
    values = shares.to_numpy()

    height = 1.8 + 0.24 * len(rows)
    width = 4.5 + 0.42 * len(columns)
    fig, (heat, rec_ax, z_ax) = plt.subplots(
        1, 3, figsize=(width + 3.6, height), sharey=True, layout="constrained",
        gridspec_kw={"width_ratios": [0.42 * len(columns), 1.4, 1.4]})

    image = heat.imshow(values, aspect="auto", cmap="Blues", vmin=0, vmax=100)
    for (i, j), v in np.ndenumerate(values):
        if v >= 5:
            heat.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=6.5,
                      color="white" if v > 55 else INK)
    # A white rule between groups of the same dominant cell type.
    groups = dominant.to_numpy()
    for i in np.flatnonzero(groups[1:] != groups[:-1]):
        heat.axhline(i + 0.5, color="white", lw=2.5)
    heat.set_xticks(range(len(columns)), columns, rotation=50, ha="right", fontsize=8)
    heat.set_yticks(range(len(rows)), [f"#{r}" for r in rows["rank"]], fontsize=7)
    heat.tick_params(length=0)
    for spine in heat.spines.values():
        spine.set_visible(False)
    heat.set_title("Mean share of each cell type (%)")
    fig.colorbar(image, ax=heat, location="bottom", shrink=0.5, pad=0.01, label="%")

    y = np.arange(len(rows))
    colours = np.where(rows["significant"], SIGNIFICANT, NOT_SIGNIFICANT)
    rec_ax.barh(y, rows["recurrence"], color=colours, height=0.7)
    rec_ax.set_xscale("log")
    rec_ax.set_title("recurrence R", fontweight="normal")
    rec_ax.tick_params(axis="y", length=0)
    z_ax.scatter(rows["z_score"], y, s=14, c=colours)
    z_ax.axvline(run.ceiling, color=INK, lw=1, ls="--")
    critical = run.critical_value()
    if np.isfinite(critical):
        z_ax.axvline(critical, color=INK, lw=1, ls=":")
    z_ax.set_title("z-score", fontweight="normal")
    z_ax.text(1, -0.02, "-- ceiling √m   ··· critical value", transform=z_ax.transAxes,
              ha="right", va="top", fontsize=7, color=MUTED)
    for ax in (rec_ax, z_ax):
        ax.grid(True, axis="x", lw=0.5)
    heat.set_ylim(len(rows) - 0.5, -0.5)

    which = "patterns" if args.all else "significant patterns"
    fig.suptitle(f"Composition of the top {len(rows)} {which}, averaged over the neighbourhoods "
                 f"matching each: {run.name}\n{run.parameters()}",
                 x=0.01, ha="left", fontsize=12, fontweight="bold")
    out = output_path(run, "compositions.png", args.out)
    fig.savefig(out, dpi=args.dpi)
    print(f"{len(rows)} patterns in {len(set(groups))} groups -> {out}")


if __name__ == "__main__":
    main()
