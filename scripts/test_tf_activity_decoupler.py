"""Checks that guard the R-to-Python translation's scientific choices."""

import numpy as np
import pandas as pd
import pytest
from scipy.stats import false_discovery_control, linregress

from tf_activity_decoupler import estimate_tf_activities, load_regulons, scale_rna


def test_scaling_uses_sample_sd_and_retains_missingness():
    rna = pd.DataFrame(
        [[1.0, 2.0, 3.0], [2.0, np.nan, 4.0], [np.nan, np.nan, 3.0], [5.0, 5.0, 5.0]],
        index=["complete", "partial", "excluded", "constant"],
        columns=["a", "b", "c"],
    )
    scaled = scale_rna(rna, max_missing=1)
    assert list(scaled.index) == ["complete", "partial", "constant"]
    np.testing.assert_allclose(scaled.loc["complete"], [-1, 0, 1])
    np.testing.assert_allclose(scaled.loc["partial", ["a", "c"]], [-1 / np.sqrt(2), 1 / np.sqrt(2)])
    assert np.isnan(scaled.loc["partial", "b"])
    assert scaled.loc["constant"].isna().all()
    assert scaled.attrs["preprocessing"]["n_removed_missing"] == 1


def test_cellwise_ulm_keeps_zero_targets_and_uses_all_available_genes():
    genes = pd.Index([f"g{i}" for i in range(30)])
    values = np.sin(np.arange(30) * 1.1) + np.arange(30) * 0.03
    values[0] = 0.0  # A measured mean-level z-score, not missing expression.
    rna = pd.DataFrame({"a": values, "b": values * 0.7 + 0.2}, index=genes)
    rna.loc["g0", "b"] = np.nan
    net = pd.DataFrame({
        "source": ["A"] * 20 + ["B"] * 20,
        "target": list(genes[:20]) + list(genes[5:25]),
        "weight": [1.0, -1.0] * 10 + [-1.0, 1.0] * 10,
    })
    result = estimate_tf_activities(rna, net, min_targets=20)
    assert result["target_counts"].loc["A", "a"] == 20
    assert result["target_counts"].loc["A", "b"] == 19
    assert np.isfinite(result["scores"].loc["A", "a"])
    assert np.isnan(result["scores"].loc["A", "b"])
    assert result["cell_summary"].loc["a", "zero_zscore_genes"] == 1
    for cell in rna:
        available = rna[cell].dropna()
        tested = result["scores"][cell].dropna()
        for tf, score in tested.items():
            predictor = net[net.source == tf].set_index("target")["weight"].reindex(available.index, fill_value=0)
            regression = linregress(predictor.to_numpy(), available.to_numpy())
            np.testing.assert_allclose(score, regression.slope / regression.stderr, rtol=1e-12, atol=1e-12)
            np.testing.assert_allclose(result["p_values"].loc[tf, cell], regression.pvalue, rtol=1e-12)
        expected_bh = false_discovery_control(result["p_values"].loc[tested.index, cell])
        # Decoupler's native BH helper stores adjusted values at float32 precision.
        np.testing.assert_allclose(result["p_adjusted"].loc[tested.index, cell], expected_bh, rtol=2e-7)


def test_ambiguous_network_duplicates_are_not_silently_averaged(tmp_path):
    path = tmp_path / "conflicting.tsv"
    path.write_text("source\ttarget\tscore\nA\tg0\t1\nA\tg0\t-1\n")
    with pytest.raises(ValueError, match="Duplicate source"):
        load_regulons(path)
