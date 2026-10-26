"""The test's statistics for one run: scores, the null, p-values.

    python scripts/plot_statistics.py RUN_DIR

Writes RUN_DIR/plots/statistics.png with six panels:

* the z-scores of the observed patterns, significant ones in blue;
* the best score T of every shuffled labelling, the null the observed
  scores are judged against, with the observed labelling's own T;
* the p-values, with how many sit at the smallest possible, 1/(m+1);
* z-score against recurrence;
* recurrence against frequency: how many of a pattern's occurrences survive
  when they must not share a cell;
* how many patterns each labelling has, observed against shuffled.

The critical value is the score a pattern must exceed to be significant, and
the ceiling sqrt(m) the largest any score can be.
"""

from __future__ import annotations

import argparse
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pam_results import (INK, MUTED, NOT_SIGNIFICANT, NULL, SIGNIFICANT, load_run, output_path,
                         use_style)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("run", help="a pam_st --output-dir")
    p.add_argument("--out", help="figure path (default: RUN_DIR/plots/statistics.png)")
    p.add_argument("--dpi", type=int, default=200)
    return p.parse_args()


def mark(ax, x, text, style, colour=INK, left=False) -> None:
    """A vertical reference line with its name written along it, on its right
    or, with left=True, on its left."""
    if not np.isfinite(x):
        return
    ax.axvline(x, color=colour, lw=1, ls=style)
    ax.text(x, 0.97, f"{text} " if left else f" {text}", transform=ax.get_xaxis_transform(),
            color=colour, fontsize=8, va="top", ha="right" if left else "left", rotation=90)


def scatter_by_significance(ax, x, y, significant) -> None:
    ax.scatter(x[~significant], y[~significant], s=9, c=NOT_SIGNIFICANT, linewidths=0,
               label="not significant")
    ax.scatter(x[significant], y[significant], s=9, c=SIGNIFICANT, linewidths=0,
               label="significant")


def no_data(ax, text) -> None:
    ax.text(0.5, 0.5, text, transform=ax.transAxes, ha="center", va="center", color=MUTED)
    ax.set_xticks([])
    ax.set_yticks([])


def main() -> None:
    args = parse_args()
    use_style()
    run = load_run(args.run)
    pats = run.patterns
    if pats.empty:
        raise SystemExit("patterns.csv is empty: the observed tissue has no pattern.")
    z = pats["z_score"].to_numpy(float)
    p = pats["p_value"].to_numpy(float)
    sig = pats["significant"].to_numpy(bool)
    rec = pats["recurrence"].to_numpy(float)
    freq = pats["frequency"].to_numpy(float)
    ceiling, critical = run.ceiling, run.critical_value()
    pool = run.m + 1
    alpha = run.number("alpha")
    lab = run.labellings
    best = lab["best_z"].to_numpy(float) if lab is not None else None

    fig, axes = plt.subplots(2, 3, figsize=(15, 9), layout="constrained")
    (a, b, c), (d, e, f) = axes

    # The observed scores.
    finite_best = best[np.isfinite(best)] if best is not None else np.array([])
    low = min(z.min(), finite_best.min() if finite_best.size else z.min())
    bins = np.linspace(low, ceiling, 41)
    a.hist([z[~sig], z[sig]], bins=bins, stacked=True, color=[NOT_SIGNIFICANT, SIGNIFICANT],
           label=[f"not significant ({(~sig).sum()})", f"significant ({sig.sum()})"])
    mark(a, ceiling, "ceiling √m", "--")
    mark(a, critical, "critical value", ":")
    a.set_xlabel("z-score")
    a.set_ylabel("observed patterns")
    a.set_title("Scores of the observed patterns")
    a.legend(loc="upper left", fontsize=8)

    # The null: each labelling's best score.
    if best is None:
        no_data(b, "labellings.csv is missing:\nthe run predates it")
    else:
        shuffled = best[1:][np.isfinite(best[1:])]
        b.hist(shuffled, bins=bins, color=NULL, label=f"shuffled labellings ({len(best) - 1})")
        mark(b, best[0], "observed T⁽⁰⁾", "-", SIGNIFICANT, left=True)
        mark(b, ceiling, "ceiling √m", "--")
        mark(b, critical, "critical value", ":")
        empty = int((~np.isfinite(best[1:])).sum())
        if empty:
            b.text(0.02, 0.8, f"{empty} shuffles have no pattern (T = −∞)",
                   transform=b.transAxes, fontsize=8, color=MUTED)
        b.set_xlim(a.get_xlim())
        b.set_xlabel("best score T⁽ʲ⁾ of a labelling")
        b.set_ylabel("labellings")
        b.legend(loc="upper left", fontsize=8)
    b.set_title("The null: each labelling's best score")

    # p-values.
    pbins = np.linspace(0, 1, 51)
    c.hist([p[~sig], p[sig]], bins=pbins, stacked=True, color=[NOT_SIGNIFICANT, SIGNIFICANT])
    mark(c, alpha, f"α = {alpha:g}", ":")
    at_floor = int(np.isclose(p, 1 / pool).sum())
    c.text(0.98, 0.95, f"{at_floor} patterns at the smallest p,\n1/(m+1) = {1 / pool:.2g}",
           transform=c.transAxes, ha="right", va="top", fontsize=8)
    c.set_xlabel("p-value")
    c.set_ylabel("observed patterns")
    c.set_title("p-values")

    # Score against recurrence.
    scatter_by_significance(d, rec, z, sig)
    d.axhline(ceiling, color=INK, lw=1, ls="--")
    if np.isfinite(critical):
        d.axhline(critical, color=INK, lw=1, ls=":")
    d.set_xscale("log")
    d.set_xlabel("recurrence R")
    d.set_ylabel("z-score")
    d.set_title("Score against recurrence")
    d.legend(loc="lower right", fontsize=8)

    # Recurrence against frequency.
    scatter_by_significance(e, freq, rec, sig)
    top = max(freq.max(), rec.max()) * 1.2
    e.plot([1, top], [1, top], color=INK, lw=1, ls="--")
    e.text(top, top, "R = f ", ha="right", va="bottom", fontsize=8)
    e.set_xscale("log")
    e.set_yscale("log")
    e.set_xlabel("frequency f (all occurrences)")
    e.set_ylabel("recurrence R (sharing no cell)")
    e.set_title("How many occurrences stay apart")

    # Patterns per labelling.
    if lab is None:
        no_data(f, "labellings.csv is missing:\nthe run predates it")
    else:
        counts = lab["patterns"].to_numpy()
        f.hist(counts[1:], bins=30, color=NULL, label="shuffled labellings")
        mark(f, counts[0], f"observed: {counts[0]}", "-", SIGNIFICANT, left=True)
        share = (counts[1:] >= counts[0]).mean()
        f.text(0.02, 0.86, f"{share:.0%} of shuffles have\nat least as many",
               transform=f.transAxes, ha="left", va="top", fontsize=8)
        f.set_xlabel("patterns |P⁽ʲ⁾| in a labelling")
        f.set_ylabel("labellings")
        f.legend(loc="upper left", fontsize=8)
    f.set_title("How many patterns recur")

    fig.suptitle(f"PaM-ST statistics: {run.name}\n{run.parameters()}   ·   "
                 f"{len(pats)} observed patterns, {sig.sum()} significant, "
                 f"universe {int(run.number('universe_size')):,}",
                 x=0.01, ha="left", fontsize=12, fontweight="bold")
    out = output_path(run, "statistics.png", args.out)
    fig.savefig(out, dpi=args.dpi)
    print(f"{len(pats)} patterns, {sig.sum()} significant, critical value "
          f"{critical:.3f}, ceiling {ceiling:.3f} -> {out}")


if __name__ == "__main__":
    main()
