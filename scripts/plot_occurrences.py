"""Where a pattern occurs: every matching neighbourhood circled in red.

    python scripts/plot_occurrences.py RUN_DIR                # the top pattern
    python scripts/plot_occurrences.py RUN_DIR --rank 1 2 5
    python scripts/plot_occurrences.py RUN_DIR --pattern "Other_cancer:27;Proliferating_Other_cancer:1"

For each pattern it writes two figures to RUN_DIR/plots/:

* occurrences_rank<N>.png: every sample, the cells of the pattern's labels in
  colour and all other cells in grey, with a red circle of radius r around the
  centre of every occurrence;
* occurrences_rank<N>_zoom.png: some occurrences up close, picked far apart
  (--zoom sets how many; 0 skips it).

occurrences.csv lists the occurrences of significant patterns only. The y axis
points down, as in the images the coordinates come from.
"""

from __future__ import annotations

import argparse
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from pam_results import (DARK, FAINT, INK, MUTED, RED, label_colours, load_cells, load_run,
                         output_path, parse_counts, pattern_title, short_names, use_style)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("run", help="a pam_st --output-dir")
    which = p.add_mutually_exclusive_group()
    which.add_argument("--rank", type=int, nargs="+", help="ranks in patterns.csv (default: 1)")
    which.add_argument("--pattern", nargs="+", help="patterns as written in patterns.csv")
    p.add_argument("--input", action="append",
                   help="input table, repeated in pam_st's order, if summary.csv's paths moved")
    p.add_argument("--zoom", type=int, default=12, help="occurrences shown up close (default: 12)")
    p.add_argument("--out-dir", help="where to write the figures (default: RUN_DIR/plots)")
    p.add_argument("--dpi", type=int, default=200)
    return p.parse_args()


def chosen_ranks(run, args) -> list[int]:
    if args.pattern:
        wanted = [parse_counts(text) for text in args.pattern]
        found = {i: row["rank"] for i, row in run.patterns.iterrows()
                 if parse_counts(row["pattern"]) in wanted}
        if len(found) < len(wanted):
            raise SystemExit("Some --pattern is not in patterns.csv; copy it from the pattern column.")
        return list(found.values())
    return args.rank or [1]


HALO = [pe.withStroke(linewidth=3.2, foreground="white")]
MIN_RING = 9.0  # points: rings smaller than this on the page are drawn at this size


def legend_handles(colours: dict[str, str], radius: float, enlarged: bool) -> list[Line2D]:
    dot = dict(marker="o", linestyle="", markersize=6)
    handles = [Line2D([], [], color=c, label=label, **dot)
               for label, c in colours.items() if c != DARK]
    if DARK in colours.values():
        handles.append(Line2D([], [], color=DARK, label="its other labels", **dot))
    handles.append(Line2D([], [], color=FAINT, label="other cells", **dot))
    note = f"occurrence, radius r = {radius:g}" + (" (drawn larger to be seen)" if enlarged else "")
    handles.append(Line2D([], [], color=RED, marker="o", markerfacecolor="none", linestyle="",
                          markersize=9, markeredgewidth=1.4, label=note))
    return handles


def draw_cells(ax, cells, colours, size_other, size_pattern, alpha=1.0) -> None:
    pattern_cells = cells["label"].isin(colours.keys())
    other = cells[~pattern_cells]
    ax.scatter(other["x"], other["y"], s=size_other, c=FAINT, linewidths=0, rasterized=True)
    # The leading labels last, so they stay on top.
    for label in reversed(list(colours)):
        part = cells[cells["label"] == label]
        ax.scatter(part["x"], part["y"], s=size_pattern, c=colours[label], linewidths=0,
                   alpha=alpha, rasterized=True)


def circles(ax, xs, ys, radius, width) -> bool:
    """Red rings of radius `radius` around the centres, or of MIN_RING points
    across when that would be smaller. Call once the limits are final; says
    whether the rings were enlarged."""
    ax.figure.canvas.draw()
    x0, x1 = ax.get_xlim()
    points_per_unit = ax.get_window_extent().width * 72 / ax.figure.dpi / abs(x1 - x0)
    diameter = 2 * radius * points_per_unit
    ax.scatter(xs, ys, s=max(diameter, MIN_RING) ** 2, marker="o", facecolors="none",
               edgecolors=RED, linewidths=width, path_effects=HALO, zorder=3)
    return diameter < MIN_RING


def spread_out(points: np.ndarray, k: int) -> list[int]:
    """k indices of points far apart: farthest-point sampling from the first."""
    chosen = [0]
    distance = np.linalg.norm(points - points[0], axis=1)
    while len(chosen) < min(k, len(points)):
        nxt = int(np.argmax(distance))
        chosen.append(nxt)
        distance = np.minimum(distance, np.linalg.norm(points - points[nxt], axis=1))
    return chosen


def dot_size(cells, panel_points: float) -> float:
    """Marker area for a tissue panel: about half the mean cell spacing across."""
    xr = max(np.ptp(cells["x"]), 1.0)
    yr = max(np.ptp(cells["y"]), 1.0)
    spacing = math.sqrt(xr * yr / max(len(cells), 1))
    diameter = 0.5 * spacing * panel_points / max(xr, yr)
    return float(np.clip(diameter, 0.6, 4.0)) ** 2


def overview(run, row, occ, cells, colours, path, dpi) -> None:
    radius = run.number("radius")
    names = list(cells)
    short = short_names(names)
    ncols = 1 if len(names) == 1 else 2
    nrows = math.ceil(len(names) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(7.5 * ncols, 7.6 * nrows), squeeze=False,
                             layout="constrained")
    enlarged = False
    for ax, name in zip(axes.flat, names):
        here = occ[occ["sample"] == name]
        size = dot_size(cells[name], 7.0 * 72)
        draw_cells(ax, cells[name], colours, size, 2.2 * size, alpha=0.6)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        enlarged |= circles(ax, here["x"], here["y"], radius, 1.3)
        ax.set_xlabel("x", color=MUTED)
        ax.set_ylabel("y", color=MUTED)
        n = len(here)
        ax.set_title(f"{short[name]}: {n} occurrence{'s' if n != 1 else ''}", fontweight="normal")
    for ax in axes.flat[len(names):]:
        ax.axis("off")
    fig.suptitle(f"Pattern #{row['rank']}: {pattern_title(row)}\n"
                 f"frequency f = {row['frequency']}   recurrence R = {row['recurrence']}   "
                 f"z = {row['z_score']:.2f}   p = {row['p_value']:.4g}   "
                 f"{len(occ)} occurrences, each circled in red",
                 x=0.02, ha="left", fontsize=12, fontweight="bold")
    fig.legend(handles=legend_handles(colours, radius, enlarged), loc="outside lower left",
               ncol=3, fontsize=9, handletextpad=0.3)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def zoom(run, row, occ, cells, colours, k, path, dpi) -> None:
    radius = run.number("radius")
    names = list(cells)
    short = short_names(names)
    # Samples have their own coordinates: set them far apart before spreading.
    offset = occ["sample"].map({n: i for i, n in enumerate(names)}).to_numpy() * 1e12
    points = np.column_stack([occ["x"].to_numpy() + offset, occ["y"].to_numpy()])
    picks = occ.iloc[spread_out(points, k)]
    ncols = min(4, len(picks))
    nrows = math.ceil(len(picks) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.1 * ncols, 3.5 * nrows + 1.0),
                             squeeze=False, layout="constrained")
    half = 1.6 * radius
    for ax, (_, o) in zip(axes.flat, picks.iterrows()):
        c = cells[o["sample"]]
        window = c[(c["x"].between(o["x"] - half, o["x"] + half)) &
                   (c["y"].between(o["y"] - half, o["y"] + half))]
        draw_cells(ax, window, colours, 14, 18)
        ax.plot(o["x"], o["y"], marker="+", color=INK, markersize=7)
        inside = ((window["x"] - o["x"]) ** 2 + (window["y"] - o["y"]) ** 2) <= radius ** 2
        ax.set_xlim(o["x"] - half, o["x"] + half)
        ax.set_ylim(o["y"] + half, o["y"] - half)
        ax.set_aspect("equal")
        circles(ax, [o["x"]], [o["y"]], radius, 1.6)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color(FAINT)
        ax.set_title(f"{short[o['sample']]} ({o['x']:.0f}, {o['y']:.0f})\n"
                     f"{int(inside.sum())} cells in the circle", fontsize=9, fontweight="normal")
    for ax in axes.flat[len(picks):]:
        ax.axis("off")
    fig.suptitle(f"Pattern #{row['rank']} up close: {len(picks)} of its {len(occ)} occurrences, "
                 f"spread out\n{pattern_title(row)}", x=0.02, ha="left", fontsize=12,
                 fontweight="bold")
    fig.legend(handles=legend_handles(colours, radius, False), loc="outside lower left", ncol=3,
               fontsize=9, handletextpad=0.3)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    use_style()
    run = load_run(args.run)
    ranks = chosen_ranks(run, args)
    significant = run.patterns[run.patterns["significant"]]
    cells = None
    for rank in ranks:
        rows = run.patterns[run.patterns["rank"] == rank]
        if rows.empty:
            raise SystemExit(f"No pattern of rank {rank}: patterns.csv has ranks 1-{len(run.patterns)}.")
        row = rows.iloc[0]
        occ = run.occurrences[run.occurrences["rank"] == rank]
        if occ.empty:
            raise SystemExit(
                f"Pattern #{rank} is not significant (p = {row['p_value']:.4g}), so occurrences.csv "
                f"does not list its occurrences; the significant ones are ranks 1-{len(significant)}.")
        if cells is None:
            cells = load_cells(run, args.input)
        counts = parse_counts(row["pattern"])
        colours = label_colours(sorted(counts, key=lambda label: -counts[label]))
        out = output_path(run, f"occurrences_rank{rank}.png",
                          f"{args.out_dir}/occurrences_rank{rank}.png" if args.out_dir else None)
        overview(run, row, occ, cells, colours, out, args.dpi)
        print(f"Pattern #{rank}: {len(occ)} occurrences -> {out}")
        if args.zoom > 0:
            zoom_out = out.with_name(f"occurrences_rank{rank}_zoom.png")
            zoom(run, row, occ, cells, colours, args.zoom, zoom_out, args.dpi)
            print(f"  up close -> {zoom_out}")


if __name__ == "__main__":
    main()
