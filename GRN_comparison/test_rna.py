"""RNA-only safeguards: full prior retention, holdout isolation and verified reuse."""

from hashlib import sha256
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "cosmos-grn-matplotlib"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(Path(tempfile.gettempdir()) / "cosmos-grn-numba"))

import anndata as ad
import gretapy as gt
import gretapy.tl._eval as gt_eval
import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

from greta_benchmark import _knocktf_support, configure_cache, evaluation_compatibility
from rna_benchmark import RNA_PATH, RnaConfig, load_rna, build_rna_networks, run_rna


@unittest.skipUnless((Path(__file__).parent / "gretapy_data/hg38_gst_collectri.csv.gz").exists(),
                     "Run the RNA practical once to cache the native CollecTRI reference.")
class RnaSafeguards(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        configure_cache()
        cls.config = RnaConfig()
        cls.adata, cls.split, cls.features = load_rna(cls.config)
        cls.networks = build_rna_networks(cls.adata, cls.split, cls.config)

    def test_full_eligible_prior_and_degree_budgets_are_preserved(self):
        prior = gt.ds.read_db("hg38", "CollecTRI").rename(columns={"weight": "score"})
        actual = self.networks["CollecTRI"]
        expected = prior[prior.source.isin(actual.source) & prior.target.isin(self.adata.var_names)
                         & prior.source.ne(prior.target)]
        self.assertEqual(set(zip(actual.source, actual.target)), set(zip(expected.source, expected.target)))
        self.assertGreater(len(actual), 5000)
        self.assertGreater(actual.source.nunique(), 100)
        self.assertEqual(self.adata.n_obs, 60)
        self.assertNotIn("cre", actual.columns)
        degrees = actual.groupby("source").size()
        for graph in self.networks.values():
            assert_frame_equal(degrees.to_frame(), graph.groupby("source").size().to_frame())

    def test_heldout_rna_cannot_change_genes_or_networks(self):
        expression = pd.read_csv(RNA_PATH, sep="\t")
        heldout = self.split.loc[self.split.partition == "test", "sample"].tolist()
        expression[heldout] += 1000
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "changed.tsv.gz"
            expression.to_csv(path, sep="\t", index=False)
            changed, split, features = load_rna(self.config, path)
        self.assertEqual(list(self.adata.var_names), list(changed.var_names))
        assert_frame_equal(self.features, features)
        for name, network in build_rna_networks(changed, split, self.config).items():
            assert_frame_equal(self.networks[name], network)

    @unittest.skipUnless((Path(__file__).parent / "gretapy_data/hg38_gst_hall.csv.gz").exists(),
                         "Hallmarks reference must already be cached.")
    def test_matching_cache_is_reused_but_changed_output_is_recomputed(self):
        with tempfile.TemporaryDirectory() as folder:
            config = RnaConfig(max_genes=1000, max_tfs=15, profile="quick")
            first = run_rna(config, output_dir=folder, metrics=["Hallmarks"], verbose=False)
            reused = run_rna(config, output_dir=folder, metrics=["Hallmarks"], verbose=False)
            self.assertFalse(first["cache_hit"])
            self.assertTrue(reused["cache_hit"])
            path = Path(folder) / "evaluation_ranked.tsv"
            original_hash = sha256(path.read_bytes()).hexdigest()
            path.write_text(path.read_text() + "\n")
            repaired = run_rna(config, output_dir=folder, metrics=["Hallmarks"], verbose=False)
            self.assertFalse(repaired["cache_hit"])
            self.assertEqual(sha256(path.read_bytes()).hexdigest(), original_hash)


class KnocktfSupportSafeguard(unittest.TestCase):
    def test_nonzero_target_eligibility_preserves_original_recall_universe(self):
        reference = ad.AnnData(
            np.array([[-1., -2., -3.], [-1., 0., 0.], [-1., -2., -3.]]),
            obs=pd.DataFrame({"TF": ["TF1", "TF1", "TF2"], "logFC": [-1., -1., -1.],
                              "Tissue.Type": ["Example"] * 3}, index=["usable", "too_few", "absent_tf"]),
            var=pd.DataFrame(index=["A", "B", "C"]),
        )
        data = ad.AnnData(np.ones((2, 5)), var=pd.DataFrame(index=["TF1", "TF2", "A", "B", "C"]))
        graph = pd.DataFrame({"source": ["TF1"] * 3, "target": ["A", "B", "C"], "score": [1.] * 3})
        eligible, supported, _ = _knocktf_support(graph, data, reference, None)
        self.assertEqual(len(eligible), 3)
        self.assertEqual(supported, ["usable"])

        def native_stub(metric_type, db_name, grn, db, genes, peaks, cats, adata, verbose=True):
            self.assertEqual(list(db.obs_names), ["usable"])
            return 1.0, 1.0, 1.0

        with patch.object(gt_eval, "_run_metric", native_stub):
            with evaluation_compatibility({}):
                precision, recall, f01 = gt_eval._run_metric(
                    "TF Scoring", "KnockTF (scoring)", graph, reference,
                    list(data.var_names), [], None, data, verbose=False,
                )
        self.assertEqual(precision, 1.0)
        self.assertAlmostEqual(recall, 1 / 3)
        self.assertAlmostEqual(f01, 1.01 / 3 / (0.01 + 1 / 3))


if __name__ == "__main__":
    unittest.main()
