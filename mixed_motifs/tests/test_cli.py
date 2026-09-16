#!/usr/bin/env python3
"""Integration tests for the compiled mixed-motif command.

Run with PAM_MIXED_BIN=/absolute/path/to/pam_mixed python3 -m unittest
discover -s mixed_motifs/tests -v. Tests use temporary directories only and do
not need the plotting dependencies. Alternatively pass --binary PATH here.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


FIELDS = ["Cell_ID", "X_centroid", "Y_centroid", "Cell_Type", "Sample", "Stratum", "Component"]


def group_rows(groups: int = 6, labels: str | list[str] = "AAAABBBCCC", sample: str = "sample_1", start: float = 0.0) -> list[dict]:
    rows = []
    for group in range(groups):
        for index, label in enumerate(labels):
            angle = math.tau * index / len(labels)
            rows.append({
                "Cell_ID": f"{sample}_g{group}_c{index}",
                "X_centroid": start + 100.0 * group + math.cos(angle),
                "Y_centroid": math.sin(angle),
                "Cell_Type": label,
                "Sample": sample,
                "Stratum": "tissue",
                "Component": "section",
            })
    return rows


class MixedMotifCLI(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        configured = os.environ.get("PAM_MIXED_BIN")
        default = Path(__file__).resolve().parents[1] / "build" / "pam_mixed"
        cls.binary = Path(configured).expanduser().resolve() if configured else default
        if not cls.binary.is_file():
            raise unittest.SkipTest(f"Build pam_mixed or set PAM_MIXED_BIN (looked for {cls.binary}).")

    def setUp(self) -> None:
        self.directory = tempfile.TemporaryDirectory(prefix="pam_mixed_test_")
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.counter = 0

    def write_input(self, rows: list[dict], fields: list[str] | None = None) -> Path:
        self.counter += 1
        path = self.root / f"input_{self.counter}.csv"
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields or FIELDS, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        return path

    def run_cli(self, rows: list[dict], extra: tuple[str, ...] = (), *, fields: list[str] | None = None, success: bool = True) -> Path:
        source = self.write_input(rows, fields)
        destination = self.root / f"result_{self.counter}"
        command = [str(self.binary), "--input", str(source), "--output", str(destination),
                   "--k", "10", "--max-radius", "5", "--min-types", "2",
                   "--max-dominance", "0.7", "--max-relabels", "0",
                   "--min-support", "1", "--permutations", "19", "--seed", "37",
                   "--threads", "1", "--null", "global", "--region-size", "25", *extra]
        result = subprocess.run(command, text=True, capture_output=True, timeout=90)
        details = f"Command: {' '.join(command)}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        if success:
            self.assertEqual(result.returncode, 0, details)
            self.assertTrue((destination / "run.json").is_file(), details)
        else:
            self.assertNotEqual(result.returncode, 0, details)
            self.assertTrue(result.stderr.strip() or result.stdout.strip(), "A rejected input should explain the error.")
        return destination

    @staticmethod
    def table(directory: Path, name: str) -> list[dict[str, str]]:
        with (directory / name).open(newline="", encoding="utf-8") as stream:
            return list(csv.DictReader(stream))

    def test_ten_cells_include_center_and_identical_sets_are_deduplicated(self) -> None:
        output = self.run_cli(group_rows())
        motifs = self.table(output, "motifs.csv")
        self.assertEqual(len(motifs), 1)
        self.assertEqual(int(motifs[0]["support"]), 6)
        self.assertEqual(int(motifs[0]["raw_matches"]), 6, "Ten centers viewing one identical group should count as one unique window.")
        self.assertEqual(len(self.table(output, "neighborhoods.csv")), 6)
        occurrences = self.table(output, "occurrences.csv")
        members = self.table(output, "members.csv")
        self.assertEqual(len(occurrences), 6)
        self.assertEqual(len(members), 60)
        for occurrence in occurrences:
            selected = [row["cell_index"] for row in members if row["motif_id"] == occurrence["motif_id"] and row["occurrence_id"] == occurrence["occurrence_id"]]
            self.assertEqual(len(selected), 10)
            self.assertEqual(len(set(selected)), 10)
            self.assertIn(occurrence["center_index"], selected)
        self.assertEqual(len({row["cell_index"] for row in members}), 60)

    def test_overlapping_windows_do_not_become_independent_occurrences(self) -> None:
        rows = group_rows(groups=1, labels="AAAABBBBCCC")
        output = self.run_cli(rows, ("--max-relabels", "1"))
        motifs = self.table(output, "motifs.csv")
        self.assertTrue(motifs)
        self.assertGreater(len(self.table(output, "neighborhoods.csv")), 1)
        self.assertTrue(all(int(motif["support"]) == 1 for motif in motifs), "Any two ten-cell sets from eleven cells overlap.")
        self.assertTrue(any(int(motif["raw_matches"]) > 1 for motif in motifs))

    def test_pure_groups_are_ineligible_even_when_frequent(self) -> None:
        rows = group_rows(groups=4, labels="AAAAAAAAAA")
        rows += group_rows(groups=4, labels="BBBBBBBBBB", sample="sample_2")
        output = self.run_cli(rows)
        self.assertEqual(self.table(output, "motifs.csv"), [])
        self.assertEqual(self.table(output, "occurrences.csv"), [])

    def test_every_match_must_satisfy_dominance_limit(self) -> None:
        rows = group_rows(groups=3, labels="AAAAAAABBB")
        rows += group_rows(groups=4, labels="AAAAAAAABB", start=1000.0)
        # The 8A+2B windows are one relabel away but fail the 70% rule.
        for index, row in enumerate(rows):
            row["Cell_ID"] = f"cell_{index}"
        output = self.run_cli(rows, ("--max-relabels", "1"))
        motifs = self.table(output, "motifs.csv")
        self.assertEqual(len(motifs), 1)
        self.assertEqual(int(motifs[0]["raw_matches"]), 3)
        self.assertEqual(int(motifs[0]["support"]), 3)

    def test_matching_tolerance_counts_relabelled_cells(self) -> None:
        rows = group_rows(groups=3, labels="AAAABBBCCC")
        rows += group_rows(groups=3, labels="AAAAAABCCC", start=1000.0)
        for index, row in enumerate(rows):
            row["Cell_ID"] = f"cell_{index}"
        one_change = self.run_cli(rows, ("--max-relabels", "1"))
        two_changes = self.run_cli(rows, ("--max-relabels", "2"))
        self.assertEqual(max(int(row["support"]) for row in self.table(one_change, "motifs.csv")), 3)
        self.assertEqual(max(int(row["support"]) for row in self.table(two_changes, "motifs.csv")), 6)

    def test_required_types_do_not_disappear_under_matching_tolerance(self) -> None:
        rows = group_rows(groups=3, labels="AAAAABBBBC")
        rows += group_rows(groups=3, labels="AAAAABBBBB", start=1000.0)
        for index, row in enumerate(rows):
            row["Cell_ID"] = f"cell_{index}"
        constrained = self.run_cli(rows, ("--max-relabels", "1"))
        permissive = self.run_cli(rows, ("--max-relabels", "1", "--allow-type-turnover"))
        self.assertEqual(max(int(row["support"]) for row in self.table(constrained, "motifs.csv")), 3)
        self.assertEqual(max(int(row["support"]) for row in self.table(permissive, "motifs.csv")), 6)

    def test_contact_filter_distinguishes_equal_composition_arrangements(self) -> None:
        rows = group_rows(groups=3, labels="AAAAABBBBB")
        rows += group_rows(groups=3, labels="ABABABABAB", start=1000.0)
        for index, row in enumerate(rows):
            row["Cell_ID"] = f"cell_{index}"
        # On these regular ten-cell rings only adjacent vertices are within 0.7.
        # Separated halves have 2/10 cross-type edges; alternating labels 10/10.
        unfiltered = self.run_cli(rows, ("--contact-radius", "0.7"))
        filtered = self.run_cli(rows, ("--contact-radius", "0.7", "--min-cross-edge-fraction", "0.5"))
        self.assertEqual(int(self.table(unfiltered, "motifs.csv")[0]["support"]), 6)
        self.assertEqual(int(self.table(filtered, "motifs.csv")[0]["support"]), 3)

    def test_physical_cutoff_excludes_undersized_groups(self) -> None:
        rows = group_rows(groups=1, labels="AAABBBCCC")
        distant = group_rows(groups=1, labels="A", start=1000.0)[0]
        distant["Cell_ID"] = "distant_cell"
        rows.append(distant)
        output = self.run_cli(rows)
        self.assertEqual(self.table(output, "motifs.csv"), [])
        self.assertEqual(self.table(output, "neighborhoods.csv"), [])

    def test_samples_and_components_are_geometric_barriers(self) -> None:
        for barrier in ("Sample", "Component"):
            with self.subTest(barrier=barrier):
                rows = group_rows(groups=1)
                for index, row in enumerate(rows):
                    row[barrier] = "first" if index < 5 else "second"
                output = self.run_cli(rows)
                self.assertEqual(self.table(output, "neighborhoods.csv"), [])
                self.assertEqual(self.table(output, "motifs.csv"), [])

    def test_stratified_shuffle_conditions_on_stratum_counts(self) -> None:
        rows = group_rows()
        for row in rows:
            row["Stratum"] = row["Cell_Type"]
        output = self.run_cli(rows, ("--null", "stratified"))
        motifs = self.table(output, "motifs.csv")
        self.assertEqual(len(motifs), 1)
        self.assertEqual(float(motifs[0]["p_search"]), 1.0)
        null = self.table(output, "null.csv")
        self.assertEqual(len(null), 19)
        self.assertTrue(all(int(row["max_support"]) == 6 for row in null))

    def test_block_shuffle_preserves_counts_in_spatial_blocks(self) -> None:
        # Ring centers lie on x-block boundaries. All cells left of a center
        # are A, and all cells right of it are B; each block is homogeneous.
        rows = group_rows(labels="BBBAAAAABB")
        output = self.run_cli(rows, ("--null", "blocks", "--block-size", "50"))
        motifs = self.table(output, "motifs.csv")
        self.assertEqual(len(motifs), 1)
        self.assertEqual(float(motifs[0]["p_search"]), 1.0)
        null = self.table(output, "null.csv")
        self.assertEqual(len(null), 19)
        self.assertTrue(all(int(row["max_support"]) == 6 for row in null))

    def test_permutation_results_do_not_depend_on_thread_count(self) -> None:
        rows = group_rows(groups=7)
        single = self.run_cli(rows)
        parallel = self.run_cli(rows, ("--threads", "3"))
        self.assertEqual(self.table(single, "null.csv"), self.table(parallel, "null.csv"))
        self.assertEqual(self.table(single, "motif_null.csv"), self.table(parallel, "motif_null.csv"))
        self.assertEqual(self.table(single, "motifs.csv"), self.table(parallel, "motifs.csv"))

    def test_search_p_values_use_complete_null_and_plus_one_correction(self) -> None:
        output = self.run_cli(group_rows())
        null = self.table(output, "null.csv")
        self.assertEqual(len(null), 19)
        for motif in self.table(output, "motifs.csv"):
            support = int(motif["support"])
            exceedances = sum(int(row["max_support"]) >= support for row in null)
            expected = (1 + exceedances) / 20
            self.assertAlmostEqual(float(motif["p_search"]), expected, places=8)
            self.assertGreaterEqual(float(motif["p_search"]), 0.05)

    def test_frozen_templates_are_evaluated_on_other_cells(self) -> None:
        discovery = self.run_cli(group_rows())
        validation_rows = group_rows(groups=4, sample="held_out", start=250.0)
        validation = self.run_cli(validation_rows, ("--templates", str(discovery / "motif_composition.csv")))
        discovered = self.table(discovery, "motifs.csv")
        validated = self.table(validation, "motifs.csv")
        self.assertEqual(len(discovered), 1)
        self.assertEqual(len(validated), 1)
        self.assertEqual(validated[0]["motif_id"], discovered[0]["motif_id"])
        self.assertEqual(int(validated[0]["support"]), 4)
        self.assertEqual(len(self.table(validation, "members.csv")), 40)

    def test_validation_retains_unsupported_frozen_templates(self) -> None:
        discovery = self.run_cli(group_rows())
        rows = []
        for index, label in enumerate("ABC"):
            rows.extend(group_rows(groups=2, labels=label * 10, sample=f"held_out_{index}"))
        validation = self.run_cli(rows, ("--templates", str(discovery / "motif_composition.csv"), "--min-support", "5"))
        motifs = self.table(validation, "motifs.csv")
        self.assertEqual(len(motifs), 1, "Fixed templates must not be screened out by their held-out support.")
        self.assertEqual(int(motifs[0]["support"]), 0)
        self.assertEqual(float(motifs[0]["p_value"]), 1.0)
        self.assertEqual(float(motifs[0]["p_value_holm"]), 1.0)

    def test_validation_handles_a_required_type_absent_from_held_out_data(self) -> None:
        discovery = self.run_cli(group_rows())
        validation = self.run_cli(group_rows(groups=4, labels="AAAAABBBBB", sample="held_out"),
                                  ("--templates", str(discovery / "motif_composition.csv")))
        motifs = self.table(validation, "motifs.csv")
        self.assertEqual(len(motifs), 1)
        self.assertEqual(int(motifs[0]["support"]), 0)
        self.assertEqual(float(motifs[0]["p_value"]), 1.0)

    def test_invalid_coordinates_and_missing_headers_fail_clearly(self) -> None:
        for coordinate in ("nan", "inf", "not_a_number"):
            with self.subTest(coordinate=coordinate):
                rows = group_rows(groups=1)
                rows[0]["X_centroid"] = coordinate
                self.run_cli(rows, success=False)
        self.run_cli(group_rows(groups=1), fields=[field for field in FIELDS if field != "Cell_Type"], success=False)
        self.run_cli(group_rows(groups=1), ("--k", "0"), success=False)

    def test_quoted_csv_labels_are_preserved(self) -> None:
        rows = group_rows(labels=["A, activated"] * 4 + ["B \"memory\""] * 3 + ["C"] * 3)
        output = self.run_cli(rows)
        motifs = self.table(output, "motifs.csv")
        self.assertEqual(len(motifs), 1)
        self.assertEqual(int(motifs[0]["support"]), 6)
        composition = self.table(output, "motif_composition.csv")
        all_values = {value for row in composition for value in row.values()}
        self.assertIn("A, activated", all_values)
        self.assertIn("B \"memory\"", all_values)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--binary", type=Path)
    arguments, remaining = parser.parse_known_args()
    if arguments.binary:
        os.environ["PAM_MIXED_BIN"] = str(arguments.binary)
    unittest.main(argv=[__file__, *remaining])
