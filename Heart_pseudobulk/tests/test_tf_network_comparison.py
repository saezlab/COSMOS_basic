"""Guard source coverage, signed-network cleaning, and descriptive score subtraction."""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "01_tf_activity.py"
spec = importlib.util.spec_from_file_location("heart_tf_comparison", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_scenic_sign_conflicts_and_same_sign_evidence(tmp_path):
    raw = pd.DataFrame([
        ["UNMEASURED_TF", "g1", "activator", .2, 2],
        ["UNMEASURED_TF", "g1", "activator", .8, 1],
        ["UNMEASURED_TF", "g2", "repressor", .7, 3],
        ["CONFLICT", "g3", "activator", .9, 20],
        ["CONFLICT", "g3", "repressor", .01, 1],
    ], columns=["TF", "gene", "type", "score", "n_topics"])
    path = tmp_path / "grn.tsv"
    raw.to_csv(path, sep="\t", index=False)
    network, conflicts, audit = module.load_scenicplus(path)
    assert network.source.unique().tolist() == ["UNMEASURED_TF"]
    assert network.set_index("target").weight.to_dict() == {"g1": .8, "g2": -.7}
    assert len(conflicts) == 2  # Exclude both signs regardless of relative strength.
    assert audit["opposing_sign_pairs_excluded"] == 1
    assert audit["same_sign_duplicate_rows_excluded"] == 1
    raw.loc[0, "score"] = -1
    raw.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="positive"):
        module.load_scenicplus(path)


def test_full_scenic_prior_keeps_sources_without_measured_own_rna():
    network, conflicts, audit = module.load_scenicplus(SCRIPT.parent / "data/scenicplus_grn.tsv")
    assert audit["raw_rows"] == 35268
    assert audit["opposing_sign_pairs_excluded"] == 594
    assert len(conflicts) == 1188
    assert len(network) == 34080
    assert network.source.nunique() == 102
    assert not network.duplicated(["source", "target"]).any()
    # Activity comes from targets: do not apply the GRETA common-source filter.
    genes = pd.Index([f"g{i}" for i in range(40)])
    values = pd.DataFrame({"A": np.linspace(-2, 2, 40), "B": np.linspace(2, -2, 40)}, index=genes)
    prior = pd.DataFrame({"source": ["ABSENT_SOURCE"] * 20, "target": genes[:20], "weight": np.ones(20)})
    result = module.estimate_tf_activities(values, prior, min_targets=20)
    assert "ABSENT_SOURCE" in result["scores"].index
    assert (result["target_counts"].loc["ABSENT_SOURCE"] == 20).all()


def test_aligned_deltas_keep_missing_values_and_sign_concordance():
    a = pd.DataFrame({"A": [2., -1., np.nan, 9.], "B": [1., 2., 3., 4.]},
                     index=["TF1", "TF2", "TF3", "ONLY_COLLECTRI"])
    # Both axes deliberately reordered to guard accidental positional subtraction.
    b = pd.DataFrame({"B": [4., 2., 1.], "A": [4., 1., 3.]}, index=["TF3", "TF2", "TF1"])
    tables, agreement = module.comparison_tables({"scores": a}, {"scores": b})
    delta = tables["delta_scenicplus_minus_collectri"]
    assert delta.index.tolist() == ["TF1", "TF2", "TF3"]
    assert delta.columns.tolist() == ["A", "B"]
    assert delta.at["TF1", "A"] == 1
    assert delta.at["TF2", "A"] == 2
    assert np.isnan(delta.at["TF3", "A"])
    assert agreement.at["A", "n_shared_tested_tfs"] == 2
    assert agreement.at["A", "sign_concordance_nonzero"] == .5
    assert agreement.at["B", "sign_concordance_nonzero"] == 1
