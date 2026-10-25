"""Focused safeguards for holdout isolation and missing-response evaluation.

Run in the teaching environment: python -m unittest discover -s GRN_comparison
"""

from dataclasses import replace
import os
from pathlib import Path
import tempfile
import unittest

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "cosmos-grn-matplotlib"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(Path(tempfile.gettempdir()) / "cosmos-grn-numba"))

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

from grn_benchmark import Config, RNA_PATH, evaluate_prediction, infer_networks, prepare_data


class BenchmarkSafeguards(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = Config()
        cls.data = prepare_data(cls.config)
        cls.networks = infer_networks(cls.data, cls.config)

    def test_heldout_values_cannot_change_training_or_networks(self):
        expression = pd.read_csv(RNA_PATH).set_index("gene")
        test_names = self.data["test_raw"].index
        # Change both held-out values and their missingness pattern drastically.
        expression[test_names] = expression[test_names].fillna(-1000) + 10000
        with tempfile.TemporaryDirectory() as directory:
            altered_path = Path(directory) / "rna.csv"
            expression.reset_index().to_csv(altered_path, index=False)
            altered = prepare_data(self.config, rna_path=altered_path)
        self.assertEqual(self.data["targets"], altered["targets"])
        self.assertEqual(self.data["tfs"], altered["tfs"])
        assert_frame_equal(self.data["train"], altered["train"])
        assert_frame_equal(self.data["gene_info"], altered["gene_info"])
        for name, network in infer_networks(altered, self.config).items():
            assert_frame_equal(self.networks[name], network)

    def test_sampling_seed_changes_only_sampled_supports(self):
        changed = infer_networks(self.data, replace(self.config, network_seed=43))
        for name in ["Pearson", "Spearman"]:
            assert_frame_equal(self.networks[name], changed[name])
        for name in ["CollecTRI", "Random"]:
            self.assertFalse(self.networks[name].equals(changed[name]))
        for networks in [self.networks, changed]:
            for network in networks.values():
                self.assertEqual(len(network), 200)
                self.assertFalse(network.duplicated(["source", "target"]).any())
                self.assertFalse((network.source == network.target).any())
                self.assertTrue(network.groupby("source").size().eq(10).all())

    def test_missing_responses_are_excluded_and_parentless_targets_retained(self):
        train_raw = pd.DataFrame({"A": [1., 2., 3., 4., 5., 6.],
                                  "B": [1., 2., 3., 4., 5., 6.],
                                  "C": [1., 2., 3., 4., 5., 6.]})
        test_raw = pd.DataFrame({"A": [2., 3., np.nan, 4., 5., 6.],
                                 "B": [2., 3., np.nan, 4., 5., np.nan],
                                 "C": [2., 2., 2., 2., 2., 2.]})
        data = dict(train_raw=train_raw, test_raw=test_raw,
                    train=train_raw.fillna(9999), test=test_raw.fillna(9999),
                    targets=["A", "B", "C"])
        empty_network = pd.DataFrame(columns=["source", "target", "score"])
        scores, summary = evaluate_prediction(data, {"NoParents": empty_network})
        scores = scores.set_index("target")
        self.assertEqual(scores.loc["A", "n_test_observed"], 5)
        self.assertAlmostEqual(scores.loc["A", "mse"], np.mean((np.array([2, 3, 4, 5, 6]) - 3.5) ** 2))
        self.assertEqual(scores.loc["A", "mse_ratio_to_mean"], 1.)
        self.assertTrue(np.isnan(scores.loc["B", "mse"]))
        self.assertTrue(np.isnan(scores.loc["C", "r2"]))
        self.assertEqual(summary.iloc[0].n_targets, 3)
        self.assertEqual(summary.iloc[0].n_targets_scored, 2)
        self.assertEqual(summary.iloc[0].n_targets_with_parents, 0)
        self.assertEqual(summary.iloc[0].mean_skill_vs_train_mean, 0.)


if __name__ == "__main__":
    unittest.main()
