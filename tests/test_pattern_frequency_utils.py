import os
import sys
import unittest
from pathlib import Path

os.environ.setdefault("GLOG_minloglevel", "2")
os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")

import numpy as np


scripts_dir = Path(
    os.environ.get(
        "PAM_ST_SCRIPTS_DIR",
        Path(__file__).resolve().parents[2] / "scripts",
    )
)
sys.path.insert(0, str(scripts_dir))

import utils


class NormalizedAnchorTests(unittest.TestCase):
    def setUp(self):
        self.count_vectors = np.array(
            [
                [1, 1],
                [1, 1],
                [2, 2],
                [4, 4],
                [1, 0],
                [2, 0],
                [0, 1],
            ],
            dtype=np.int64,
        )

    def test_proportional_vectors_form_one_anchor(self):
        expected_patterns = np.array([[1, 1], [1, 0], [0, 1]])
        for metric in ("l2", "js"):
            for rho in (0.0, 0.01):
                with self.subTest(metric=metric, rho=rho):
                    result = utils.count_vector_pattern_frequencies(
                        self.count_vectors, rho=rho, metric=metric)
                    np.testing.assert_array_equal(
                        result["patterns"], expected_patterns)
                    np.testing.assert_array_equal(
                        result["composition_multiplicities"], [4, 2, 1])
                    np.testing.assert_array_equal(
                        result["exact_multiplicities"], [2, 1, 1])
                    np.testing.assert_array_equal(
                        result["raw_anchor_counts"], [3, 2, 1])
                    np.testing.assert_array_equal(
                        result["within_rho_frequencies"], [4, 2, 1])
                    self.assertEqual(result["n_raw_anchors"], 6)
                    self.assertEqual(result["n_composition_anchors"], 3)

    def test_representative_and_size_ranges_follow_first_occurrence(self):
        result = utils.count_vector_pattern_frequencies(
            self.count_vectors, rho=0.0, metric="l2")
        np.testing.assert_array_equal(result["first_cell_indices"], [0, 4, 6])
        np.testing.assert_array_equal(
            result["min_neighborhood_sizes"], [2, 1, 1])
        np.testing.assert_array_equal(
            result["max_neighborhood_sizes"], [8, 2, 1])

    def test_permutation_helpers_use_the_same_composition_identity(self):
        expected_matches = np.array(
            [True, True, True, True, False, False, False])
        for metric in ("l2", "js"):
            with self.subTest(metric=metric):
                observed, pattern = utils._densest_motif(
                    self.count_vectors, rho=0.0, metric=metric)
                self.assertEqual(observed, 4)
                np.testing.assert_array_equal(pattern, [1, 1])
                np.testing.assert_array_equal(
                    utils._motif_matches(
                        self.count_vectors, pattern, rho=0.0, metric=metric),
                    expected_matches,
                )

    def test_non_integer_counts_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "integers"):
            utils.count_vector_pattern_frequencies(
                np.array([[0.5, 0.5]]), rho=0.0, metric="l2")

    def test_matched_mean_compositions_follow_metric_rho_balls(self):
        raw_proportions = (
            self.count_vectors
            / self.count_vectors.sum(axis=1, keepdims=True)
        )
        for metric in ("l2", "js"):
            for rho in (0.0, 0.2, 2.0):
                with self.subTest(metric=metric, rho=rho):
                    result = utils.count_vector_pattern_frequencies(
                        self.count_vectors,
                        rho=rho,
                        metric=metric,
                        include_match_mean=True,
                    )
                    for index, pattern in enumerate(result["patterns"]):
                        matches = utils._motif_matches(
                            self.count_vectors,
                            pattern,
                            rho=rho,
                            metric=metric,
                        )
                        self.assertEqual(
                            result["within_rho_frequencies"][index],
                            int(matches.sum()),
                        )
                        np.testing.assert_allclose(
                            result["mean_matched_proportions"][index],
                            raw_proportions[matches].mean(axis=0),
                        )


if __name__ == "__main__":
    unittest.main()
