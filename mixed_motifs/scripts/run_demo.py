#!/usr/bin/env python3
"""Generate spatial controls, mine motifs, validate, and draw result figures.

The orchestration and synthetic generation use only the Python standard
library. Figure generation uses the dependencies in requirements.txt; pass
--skip-plots to run the C++ benchmark without those packages.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys


PROJECT = Path(__file__).resolve().parents[1]
SCRIPTS = PROJECT / "scripts"
CASES = ("random_labels", "homotypic_blobs", "adjacent_domains", "planted_mixed")


def call(command: list[str], description: str) -> None:
    print(f"\n{description}", flush=True)
    subprocess.run(command, check=True)


def number(value: str | None) -> float | None:
    return float(value) if value not in (None, "") else None


def summarize(name: str, result: Path, root: Path) -> dict:
    run = json.loads((result / "run.json").read_text(encoding="utf-8"))
    with (result / "motifs.csv").open(newline="", encoding="utf-8") as stream:
        motifs = list(csv.DictReader(stream))
    search_p = [number(row.get("p_search")) for row in motifs]
    search_p = [value for value in search_p if value is not None]
    adjusted_p = [number(row.get("p_value_holm")) for row in motifs]
    adjusted_p = [value for value in adjusted_p if value is not None]
    alpha = float(run.get("alpha", 0.05))
    return {
        "case": name,
        "mode": run["mode"],
        "cells": run["cells"],
        "unique_neighborhoods": run["n_neighborhoods"],
        "eligible_neighborhoods": run["n_eligible"],
        "recurring_candidates": run["n_candidates"],
        "reported_motifs": len(motifs),
        "significant_reported_motifs": sum(row.get("significant", "").lower() == "true" for row in motifs),
        "observed_max_support": run["observed_max_support"],
        "global_p_value": run.get("global_p_value"),
        "minimum_reported_search_p": min(search_p) if search_p else None,
        "strongest_search_significant": min(search_p) <= alpha if search_p else None,
        "minimum_validation_holm_p": min(adjusted_p) if adjusted_p else None,
        "top_reported_composition": motifs[0]["composition"] if motifs else "",
        "top_reported_support": int(motifs[0]["support"]) if motifs else 0,
        "results_directory": str(result.relative_to(root)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=PROJECT / "build" / "pam_mixed")
    parser.add_argument("--output", type=Path, default=PROJECT / "examples")
    parser.add_argument("--permutations", type=int, default=199)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--seed", type=int, default=37)
    parser.add_argument("--skip-plots", action="store_true", help="Do not invoke plot_results.py.")
    parser.add_argument("--overwrite", action="store_true", help="Replace generated files in an existing demo directory; unrelated files are retained.")
    args = parser.parse_args()
    binary, root = args.binary.expanduser().resolve(), args.output.expanduser().resolve()
    if not binary.is_file():
        parser.error(f"Executable not found: {binary}. Build pam_mixed first or pass --binary.")
    if args.permutations < 1 or args.threads < 1 or args.seed < 0:
        parser.error("--permutations and --threads must be positive; --seed must be nonnegative.")
    if root.exists() and (not root.is_dir() or (any(root.iterdir()) and not args.overwrite)):
        parser.error(f"Output is not an empty directory: {root}. Choose a new directory or pass --overwrite.")
    root.mkdir(parents=True, exist_ok=True)
    common = [str(binary), "--k", "10", "--max-radius", "20", "--min-types", "2",
              "--max-dominance", "0.7", "--max-relabels", "1", "--min-support", "5",
              "--null", "global", "--region-size", "80", "--contact-radius", "8",
              "--min-cross-edge-fraction", "0", "--candidate-pool", "100",
              "--permutations", str(args.permutations), "--threads", str(args.threads),
              "--seed", str(args.seed)]
    if args.overwrite:
        common.append("--overwrite")
    try:
        call([sys.executable, str(SCRIPTS / "generate_synthetic.py"), "--output", str(root / "data"),
              "--seed", str(args.seed)], "Generating discovery controls and planted groups")
        call([sys.executable, str(SCRIPTS / "generate_synthetic.py"), "--output", str(root / "validation_data"),
              "--seed", str(args.seed + 1)], "Generating independent validation data")
        results = []
        for name in CASES:
            result = root / "results" / name
            call([*common, "--input", str(root / "data" / f"{name}.csv"), "--output", str(result)], f"Mining {name}")
            results.append((name, result))
        validation = root / "results" / "validation_planted"
        call([*common, "--input", str(root / "validation_data" / "planted_mixed.csv"),
              "--output", str(validation), "--templates",
              str(root / "results" / "planted_mixed" / "motif_composition.csv")],
             "Testing frozen discovery templates on independent planted data")
        results.append(("validation_planted", validation))
        summary = [summarize(name, result, root) for name, result in results]
        with (root / "benchmark_summary.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(summary[0]))
            writer.writeheader()
            writer.writerows(summary)
        record = {
            "discovery_seed": args.seed,
            "validation_seed": args.seed + 1,
            "permutations": args.permutations,
            "minimum_p_value": 1 / (args.permutations + 1),
            "null_model": "global random labeling within sample/component",
            "notes": [
                "These are controlled synthetic examples, not a repeated-simulation false-positive calibration study.",
                "The boundary control can reject random labeling without demonstrating intermingled cells or irreducible higher-order interaction.",
                "Discovery search p-values test the complete random-label null; enrichment ranking remains descriptive.",
                "Validation templates are frozen before testing on a separately generated dataset. Holm correction covers all supplied templates.",
                "Reported motif count can be zero even when the maximum search support is positive but below the reporting threshold.",
            ],
            "results": summary,
        }
        (root / "benchmark_summary.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
        if not args.skip_plots:
            for name, result in results:
                call([sys.executable, str(SCRIPTS / "plot_results.py"), "--results", str(result), "--top", "6"], f"Drawing {name} figures")
    except (subprocess.CalledProcessError, OSError, ValueError, KeyError) as error:
        parser.exit(1, f"Demo stopped: {error}\nCompleted outputs are retained in {root}.\n")
    print(f"\nBenchmark summary: {root / 'benchmark_summary.csv'}")
    print(f"Interpretation and settings: {root / 'benchmark_summary.json'}")


if __name__ == "__main__":
    main()
