"""Reading PaM-ST results, and the plot style the scripts share.

A run directory is what ``pam_st --output-dir DIR`` writes: summary.csv,
patterns.csv, labellings.csv and occurrences.csv. The input tables are found
from the paths in summary.csv, or given again with --input.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd

# Black text on white; colour only marks data. The categorical hues are the
# first slots of the reference palette, which stay apart for every pair and
# under colour-vision deficiencies. Red is kept for the occurrence circles.
INK = "#0b0b0b"
MUTED = "#6d6c66"
GRID = "#e4e3de"
FAINT = "#d6d5cf"   # cells that play no part
DARK = "#8a8983"    # labels folded into "other"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]
RED = "#e34948"
SIGNIFICANT = SERIES[0]
NOT_SIGNIFICANT = "#b9b8b2"
NULL = SERIES[1]
PHASES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]


def use_style() -> None:
    mpl.rcParams.update({
        "font.size": 10,
        "text.color": INK,
        "axes.labelcolor": INK,
        "axes.edgecolor": MUTED,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "grid.color": GRID,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "legend.frameon": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",
    })


@dataclass
class Run:
    path: Path
    summary: dict[str, str]
    patterns: pd.DataFrame
    occurrences: pd.DataFrame
    labellings: pd.DataFrame | None

    def number(self, key: str, default: float = math.nan) -> float:
        try:
            return float(self.summary[key])
        except (KeyError, ValueError):
            return default

    @property
    def name(self) -> str:
        return self.path.name

    @property
    def m(self) -> int:
        return int(self.number("permutations"))

    @property
    def ceiling(self) -> float:
        """sqrt(m), the largest score a pattern can reach."""
        return math.sqrt(self.m)

    @property
    def inputs(self) -> list[str]:
        return [p for p in self.summary.get("input", "").split(";") if p]

    def critical_value(self) -> float:
        """The score an observed pattern must exceed to be significant.

        p = #{j : T(j) >= z} / (m+1) <= alpha exactly when z is above the
        (K+1)-th largest T, K = floor(alpha (m+1)).
        """
        if self.labellings is None:
            return math.nan
        best = np.sort(self.labellings["best_z"].to_numpy())[::-1]
        k = int(math.floor(self.number("alpha") * len(best) + 1e-9))
        return float(best[k]) if k < len(best) else -math.inf

    def parameters(self) -> str:
        """The run's settings in one line, for figure subtitles."""
        fields = [
            f"{int(self.number('cells')):,} cells",
            f"{int(self.number('labels'))} labels",
            f"r = {self.number('radius'):g}",
            f"ρ = {self.number('rho'):g}",
            f"κ = {int(self.number('min_recurrence', 1))}",
            f"L = {self.number('tile_size'):g}",
            f"m = {self.m}",
            f"α = {self.number('alpha'):g}",
        ]
        return "   ".join(fields)


def load_run(path: str | Path) -> Run:
    path = Path(path).expanduser().resolve()
    if not (path / "summary.csv").is_file():
        raise FileNotFoundError(f"No summary.csv in {path}: is it a pam_st --output-dir?")
    summary = pd.read_csv(path / "summary.csv", dtype=str, keep_default_na=False)
    patterns = pd.read_csv(path / "patterns.csv", keep_default_na=False)
    patterns["significant"] = patterns["significant"].astype(str).str.lower() == "true"
    occurrences_file = path / "occurrences.csv"
    occurrences = (pd.read_csv(occurrences_file) if occurrences_file.is_file()
                   else pd.DataFrame(columns=["rank", "sample", "spot", "x", "y"]))
    labellings_file = path / "labellings.csv"
    labellings = pd.read_csv(labellings_file) if labellings_file.is_file() else None
    return Run(path, dict(zip(summary["metric"], summary["value"])), patterns,
               occurrences, labellings)


def find_runs(paths: list[str]) -> list[Path]:
    """Run directories among the paths, searching below any that is not one."""
    runs: list[Path] = []
    for p in map(lambda s: Path(s).expanduser().resolve(), paths):
        if (p / "summary.csv").is_file():
            runs.append(p)
        else:
            runs.extend(sorted(f.parent for f in p.rglob("summary.csv")))
    return list(dict.fromkeys(runs))


def parse_counts(text: str) -> dict[str, int]:
    """'A:2;B:1' -> {'A': 2, 'B': 1}, the format of the pattern column."""
    out: dict[str, int] = {}
    for part in filter(None, str(text).split(";")):
        label, value = part.rsplit(":", 1)
        out[label] = int(value)
    return out


def parse_composition(text: str) -> dict[str, float]:
    """'A:66.67;B:33.33' -> percentages, the format of mean_composition."""
    out: dict[str, float] = {}
    for part in filter(None, str(text).split(";")):
        label, value = part.rsplit(":", 1)
        out[label] = float(value)
    return out


def short_names(names: list[str]) -> dict[str, str]:
    """Sample names cut at the first underscore when that keeps them apart."""
    short = {n: n.split("_", 1)[0] for n in names}
    return short if len(set(short.values())) == len(names) else {n: n for n in names}


def _resolve_input(raw: str, run: Run) -> Path:
    p = Path(raw).expanduser()
    if p.is_absolute():
        candidates = [p]
    else:
        # The binary may have been started anywhere; try the usual places.
        candidates = [Path.cwd() / p, run.path / p, *[a / p for a in run.path.parents[:4]]]
    for c in candidates:
        if c.is_file():
            return c.resolve()
    raise FileNotFoundError(f"Cannot find the input table {raw!r}; pass it again with --input")


def load_cells(run: Run, inputs: list[str] | None = None) -> dict[str, pd.DataFrame]:
    """Every sample's cells, filtered as pam_st filters them, by sample name.

    Columns x, y, label. Rows are in the order pam_st numbers spots, so the
    spot column of occurrences.csv indexes them.
    """
    paths = [Path(p).expanduser().resolve() for p in inputs] if inputs else \
            [_resolve_input(p, run) for p in run.inputs]
    cells: dict[str, pd.DataFrame] = {}
    for path in paths:
        table = pd.read_csv(path, usecols=["X_centroid", "Y_centroid", "Cell_Type"])
        label = table["Cell_Type"].astype(str).str.strip()
        keep = label.str.lower() != "unclassified"
        frame = pd.DataFrame({"x": table["X_centroid"][keep].to_numpy(),
                              "y": table["Y_centroid"][keep].to_numpy(),
                              "label": label[keep].to_numpy()})
        name, copy = path.stem, 2
        while name in cells:  # the same table twice, named as pam_st names it
            name, copy = f"{path.stem}#{copy}", copy + 1
        cells[name] = frame
    return cells


def label_colours(labels: list[str]) -> dict[str, str]:
    """The first labels get the categorical hues in order; the rest share DARK."""
    return {label: SERIES[i] if i < len(SERIES) else DARK for i, label in enumerate(labels)}


def output_path(run: Run, name: str, out: str | None) -> Path:
    """`out` if given, else plots/<name> inside the run directory."""
    path = Path(out).expanduser() if out else run.path / "plots" / name
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def pattern_title(row: pd.Series) -> str:
    counts = parse_counts(row["pattern"])
    return ", ".join(f"{label} {n}" for label, n in sorted(counts.items(), key=lambda kv: -kv[1]))
