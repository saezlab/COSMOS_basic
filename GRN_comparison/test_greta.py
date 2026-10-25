"""Safeguards for the native GRETA practical; cached integration tests never download."""

import os
from pathlib import Path
import tempfile
import unittest

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "cosmos-grn-matplotlib"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(Path(tempfile.gettempdir()) / "cosmos-grn-numba"))

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

from greta_benchmark import (
    HERE, DIAGNOSTICS, TEACHING_METRICS, GretaConfig, build_networks,
    comparable_results, configure_cache, evaluate_networks, evaluation_compatibility, load_dataset, validate_grn,
)


class RankingSafeguards(unittest.TestCase):
    def test_constant_positive_response_is_unavailable_instead_of_a_crash(self):
        import anndata as ad
        import gretapy.tl._predictive as predictive
        rng = np.random.default_rng(1)
        x = rng.uniform(1, 4, 50)
        y = 2 * x + rng.normal(0, 0.2, 50)
        y[30:] = 2.0
        data = ad.AnnData(np.column_stack([x, y]),
                         obs=pd.DataFrame(index=[f"cell{i}" for i in range(50)]),
                         var=pd.DataFrame(index=["TF1", "G"]))
        graph = pd.DataFrame(dict(source=["TF1"], target=["G"], score=[1.0]))
        coverage = {}
        with evaluation_compatibility(coverage):
            frame = predictive._test_predictability(
                data=data, train_obs_names=list(data.obs_names[:30]), test_obs_names=list(data.obs_names[30:]),
                grn=graph, col_source="source", col_target="target", mod_source=None, mod_target=None, verbose=False,
            )
        self.assertTrue(frame.empty)
        self.assertEqual(coverage["Gene ~ TFs"]["n_tested"], 0)
        self.assertEqual(coverage["Gene ~ TFs"]["n_skipped_constant_response"], 1)

    def test_incomplete_columns_and_circular_checks_cannot_change_order(self):
        rows = []
        for name in ["A", "B"]:
            for db, score in [("Gene ~ TFs", 0.5), ("Hallmarks", 0.9 if name == "A" else np.nan),
                              ("CollecTRI", 1.0 if name == "A" else 0.0), ("Promoters", 1.0)]:
                rows.append(dict(name=name, organism="hg38", dataset="Example", task="Example", db=db,
                                 f01=score, **{"class": "Predictive"}))
        kept, excluded = comparable_results(pd.DataFrame(rows))
        self.assertEqual(set(kept.db), {"Gene ~ TFs"})
        self.assertEqual(set(kept.name), {"A", "B"})
        self.assertEqual(set(excluded.db), {"Hallmarks", "CollecTRI", "Promoters"})
        self.assertTrue(excluded.reason.str.len().gt(0).all())

    def test_student_network_must_have_real_matched_regions(self):
        mdata = load_dataset()
        bad = pd.DataFrame(dict(source=["PAX5"] * 5, target=["CD19"] * 5,
                                cre=["chr99-10-20"] * 5, score=[1.] * 5))
        with self.assertRaisesRegex(ValueError, "ATAC feature names"):
            validate_grn(bad, mdata)


def cache_ready():
    # Check all resources needed by these tests, including implicit context terms.
    names = ["hg38_cre_promoters.bed.gz", "hg38_gst_collectri.csv.gz", "hg38_cre_encode.bed.gz",
             "hg38_cre_blacklist.bed.gz", "hg38_tfm_hpa.tsv.gz", "hg38_tfm_tfmdb.tsv.gz",
             "hg38_tfp_intact.tsv.gz", "hg38_tfp_europmc.tsv.gz", "hg38_gst_hall.csv.gz",
             "hg38_prt_knocktf.h5ad.gz", "hg38_terms.csv.gz", "mm10_terms.csv.gz"]
    return all((HERE / "gretapy_data" / name).exists() for name in names)


@unittest.skipUnless(cache_ready(), "Run the teaching profile first to populate its reference cache.")
class CachedNativeIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import gretapy as gt
        configure_cache()
        cls.config = GretaConfig()
        cls.data = load_dataset(cls.config)
        cls.networks, cls.split, _ = build_networks(cls.data, cls.config)
        cls.results, cls.coverage = evaluate_networks(
            cls.data, cls.networks, gt.get_terms("PBMC"), list(TEACHING_METRICS) + list(DIAGNOSTICS),
        )

    def test_native_tasks_are_not_silently_omitted(self):
        self.assertEqual(len(self.results), 4 * (len(TEACHING_METRICS) + len(DIAGNOSTICS)))
        self.assertEqual(set(self.results.name), set(self.networks))
        self.assertEqual(set(self.results.db), set(TEACHING_METRICS + DIAGNOSTICS))
        # Gene Sets was silently omitted by the upstream 0.0.3 dispatcher.
        self.assertEqual(len(self.results[self.results.db == "Hallmarks"]), 4)
        kept, _ = comparable_results(self.results)
        self.assertEqual(kept["class"].nunique(), 4)
        self.assertFalse(kept.db.isin(DIAGNOSTICS).any())
        tested = self.coverage[self.coverage.db.str.contains(" ~ ", regex=False)]
        self.assertTrue(tested.n_tested.gt(0).all())

    def test_test_cells_cannot_select_demo_edges(self):
        altered = self.data.copy()
        heldout = self.split.loc[self.split.partition == "test", "cell"]
        index = altered.obs_names.get_indexer(heldout)
        altered.mod["rna"].X[index, :] += 1000
        changed, _, _ = build_networks(altered, self.config)
        for name in self.networks:
            assert_frame_equal(self.networks[name], changed[name])


if __name__ == "__main__":
    unittest.main()
