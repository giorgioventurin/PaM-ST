#!/usr/bin/env python3
"""Write reproducible spatial controls and a planted mixed-motif benchmark.

Only the Python standard library is required. Coordinates are arbitrary units;
the recommended physical cutoff for these fixtures is 20, not a tissue-derived
measurement. The truth manifest describes the construction, not test outcomes.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
import math
from pathlib import Path
import random


FIELDS = ["Cell_ID", "X_centroid", "Y_centroid", "Cell_Type", "Sample", "Stratum", "Component"]
TYPES = tuple("ABCDEF")


def compact_positions(rng: random.Random, x: float, y: float) -> list[tuple[float, float]]:
    """Ten distinct, mildly perturbed points, all at most 16 units apart."""
    rotation = rng.uniform(0.0, math.tau)
    positions = []
    for index in range(10):
        angle = rotation + math.tau * index / 10.0
        radius = rng.uniform(5.5, 7.5)
        positions.append((x + radius * math.cos(angle), y + radius * math.sin(angle)))
    return positions


def cell_record(cell_id: str, x: float, y: float, label: str, sample: str) -> dict:
    return {
        "Cell_ID": cell_id,
        "X_centroid": round(x, 6),
        "Y_centroid": round(y, 6),
        "Cell_Type": label,
        "Sample": sample,
        "Stratum": "tissue",
        "Component": "section",
    }


def compact_dataset(kind: str, seed: int, samples: int, groups: int, planted_per_sample: int) -> tuple[list[dict], list[dict]]:
    rng = random.Random(seed)
    records: list[dict] = []
    truth: list[dict] = []
    columns = math.ceil(math.sqrt(groups))
    for sample_index in range(samples):
        sample = f"sample_{sample_index + 1:02d}"
        planted_indices = set(rng.sample(range(groups), planted_per_sample)) if kind == "planted_mixed" else set()
        homotypic_labels = [TYPES[index % len(TYPES)] for index in range(groups)]
        rng.shuffle(homotypic_labels)
        for group_index in range(groups):
            x = 30.0 + 80.0 * (group_index % columns)
            y = 30.0 + 80.0 * (group_index // columns)
            is_planted = group_index in planted_indices
            if is_planted:
                labels = list("AAAABBBCCC")
                rng.shuffle(labels)
            elif kind == "homotypic_blobs":
                labels = [homotypic_labels[group_index]] * 10
            else:
                labels = [rng.choice(TYPES) for _ in range(10)]
            positions = compact_positions(rng, x, y)
            member_ids = []
            for cell_index, ((cell_x, cell_y), label) in enumerate(zip(positions, labels)):
                cell_id = f"{sample}_g{group_index + 1:03d}_c{cell_index + 1:02d}"
                member_ids.append(cell_id)
                records.append(cell_record(cell_id, cell_x, cell_y, label, sample))
            if is_planted:
                truth.append({
                    "sample": sample,
                    "group": group_index + 1,
                    "center": [x, y],
                    "composition": {"A": 4, "B": 3, "C": 3},
                    "member_ids": member_ids,
                })
    return records, truth


def adjacent_dataset(seed: int, samples: int, groups: int) -> list[dict]:
    """Six touching single-type domains; mixed neighborhoods occur at borders."""
    rng = random.Random(seed)
    count_per_sample = groups * 10
    columns = 30
    records = []
    for sample_index in range(samples):
        sample = f"sample_{sample_index + 1:02d}"
        for index in range(count_per_sample):
            column, row = index % columns, index // columns
            # Small jitter removes most distance ties without blurring domains.
            x = 30.0 + 4.0 * column + rng.uniform(-0.1, 0.1)
            y = 30.0 + 4.0 * row + rng.uniform(-0.1, 0.1)
            label = TYPES[min(len(TYPES) - 1, column // 5)]
            records.append(cell_record(f"{sample}_c{index + 1:05d}", x, y, label, sample))
    return records


def generate(output: Path, seed: int = 37, samples: int = 3, groups: int = 60, planted_per_sample: int = 24) -> dict:
    if samples < 1 or groups < 6:
        raise ValueError("Use at least one sample and six groups per sample.")
    if not 0 < planted_per_sample < groups:
        raise ValueError("Planted groups must be positive and fewer than all groups per sample.")
    output.mkdir(parents=True, exist_ok=True)
    descriptions = {
        "random_labels": "Independent uniform labels on isolated ten-cell clusters: a random-label negative control.",
        "homotypic_blobs": "Each isolated ten-cell cluster has one type: strong homotypic organization that fails mixed-type eligibility.",
        "adjacent_domains": "Six touching single-type domains on a jittered lattice: mixed groups describe boundaries without proving intermingling.",
        "planted_mixed": "Repeated isolated 4A+3B+3C groups, with independent six-type background clusters; three samples by default.",
    }
    manifest = {
        "seed": seed,
        "coordinate_units": "arbitrary synthetic units",
        "recommended_options": {
            "k": 10, "max_radius": 20, "min_types": 2,
            "max_dominance": 0.7, "max_relabels": 1,
            "min_support": 5, "null": "global", "permutations": 199,
            "region_size": 80, "contact_radius": 8,
        },
        "interpretation": [
            "The random-label control is one realization; statistical calibration requires repeated independently generated datasets.",
            "Homotypic blobs should yield no eligible observed motifs with the recommended physical cutoff.",
            "Boundary motifs in adjacent domains may be significant under random labeling; that does not establish intermingling or a higher-order interaction.",
            "Planted groups are separated by 80 units; each includes exactly ten cells with diameter below 16 units.",
            "A template near the planted composition may also match background groups; compare membership with the planted truth before estimating recovery.",
        ],
        "datasets": {},
    }
    for offset, kind in enumerate(descriptions):
        if kind == "adjacent_domains":
            records, truth = adjacent_dataset(seed + offset * 1009, samples, groups), []
        else:
            records, truth = compact_dataset(kind, seed + offset * 1009, samples, groups, planted_per_sample)
        destination = output / f"{kind}.csv"
        with destination.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(records)
        manifest["datasets"][kind] = {
            "file": destination.name,
            "description": descriptions[kind],
            "cells": len(records),
            "samples": sorted({row["Sample"] for row in records}),
            "type_counts": dict(sorted(Counter(row["Cell_Type"] for row in records).items())),
            "planted_occurrences": truth,
        }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="Directory for four CSV files and manifest.json.")
    parser.add_argument("--seed", type=int, default=37)
    parser.add_argument("--samples", type=int, default=3)
    parser.add_argument("--groups-per-sample", type=int, default=60)
    parser.add_argument("--planted-per-sample", type=int, default=24)
    args = parser.parse_args()
    try:
        manifest = generate(args.output, args.seed, args.samples, args.groups_per_sample, args.planted_per_sample)
    except ValueError as error:
        parser.error(str(error))
    for name, item in manifest["datasets"].items():
        print(f"{name}: {item['cells']} cells, {len(item['planted_occurrences'])} planted groups")
    print(f"Wrote {args.output.resolve() / 'manifest.json'}")


if __name__ == "__main__":
    main()
