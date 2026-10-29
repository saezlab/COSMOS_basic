"""Scientific design checks for the small heart GRN comparison."""
import importlib.util
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

MODULE = Path(__file__).resolve().parents[1] / "02_grn_comparison.py"
spec = importlib.util.spec_from_file_location("heart_grn_comparison", MODULE)
grn = importlib.util.module_from_spec(spec)
spec.loader.exec_module(grn)


def networks(tmp_path):
    genes = ["TF1", "TF2", *[f"G{i}" for i in range(1, 7)]]
    data = ad.AnnData(np.ones((6, len(genes))), var=pd.DataFrame(index=genes))
    rows = [("TF1", "G1", "activator", .9, 3), ("TF1", "G1", "repressor", .8, 4)]
    rows += [("TF1", f"G{i}", "activator", 1 - .1*i, i) for i in range(2, 6)]
    rows += [("TF2", f"G{i}", "repressor" if i % 2 == 0 else "activator", 1 - .1*i, i) for i in range(1, 6)]
    scenic = pd.DataFrame(rows, columns=["TF", "gene", "type", "score", "n_topics"])
    prior = pd.DataFrame([(tf, f"G{i}", (-1)**i) for tf, count in [("TF1", 3), ("TF2", 4)] for i in range(1, count+1)],
                         columns=["source", "target", "weight"])
    scenic_path, prior_path = tmp_path / "scenic.tsv", tmp_path / "prior.tsv"
    scenic.to_csv(scenic_path, sep="\t", index=False)
    prior.to_csv(prior_path, sep="\t", index=False)
    return data, scenic_path, prior_path


def test_opposite_sign_pairs_are_removed_entirely_and_matched_budget_is_fair(tmp_path):
    args = networks(tmp_path)
    full, matched, audit = grn.prepare_networks(*args, seed=42, min_targets=3)
    assert audit["conflict_pairs"] == 1
    assert audit["conflict_rows"] == 2
    assert not ((full["SCENIC+"].source == "TF1") & (full["SCENIC+"].target == "G1")).any()
    assert full["SCENIC+"].groupby("source").size().to_dict() == {"TF1": 4, "TF2": 5}
    for network in matched.values():
        assert network.groupby("source").size().to_dict() == {"TF1": 3, "TF2": 4}
        assert not network.duplicated(["source", "target"]).any()
        assert (network.source != network.target).all()
    assert set(matched["SCENIC+"].loc[matched["SCENIC+"].source == "TF1", "target"]) == {"G2", "G3", "G4"}


def test_random_preserves_each_tfs_signed_score_multiset_and_is_reproducible(tmp_path):
    args = networks(tmp_path)
    full, matched, _ = grn.prepare_networks(*args, seed=42, min_targets=3)
    rerun_full, rerun_matched, _ = grn.prepare_networks(*args, seed=42, min_targets=3)
    for collection, rerun in ((full, rerun_full), (matched, rerun_matched)):
        for name in collection:
            pd.testing.assert_frame_equal(collection[name], rerun[name])
        for tf in collection["SCENIC+"].source.unique():
            scores = collection["SCENIC+"].loc[lambda x: x.source == tf, "score"]
            random_scores = collection["Random"].loc[lambda x: x.source == tf, "score"]
            np.testing.assert_array_equal(np.sort(scores), np.sort(random_scores))


def test_knocktf_small_common_coverage_is_diagnostic_and_gate_is_configurable():
    results = pd.DataFrame([dict(name=name, db=db, precision=.5, recall=.3, f01=.49)
                            for name in ("SCENIC+", "CollecTRI", "Random")
                            for db in ("Hallmarks", "KnockTF (scoring)", "TF-Marker")])
    coverage = pd.DataFrame([dict(name=name, db="KnockTF (scoring)", n_common_tested=2)
                             for name in ("SCENIC+", "CollecTRI", "Random")])
    ranked, excluded = grn.comparable(results, coverage, min_knocktf_tests=10)
    assert set(ranked.db) == {"Hallmarks"}
    assert excluded.loc[excluded.db == "KnockTF (scoring)", "reason"].str.contains("2 common tested experiments").all()
    native = results.copy()
    ranked, _ = grn.comparable(results, coverage, min_knocktf_tests=2)
    assert set(ranked.db) == {"Hallmarks", "KnockTF (scoring)"}
    pd.testing.assert_frame_equal(results, native)
