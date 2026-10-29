"""Guard count pooling and the six-cell-type expression reference."""

import importlib.util
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

SCRIPT = Path(__file__).resolve().parents[1] / "01_tf_activity.py"
spec = importlib.util.spec_from_file_location("heart_tf", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_count_sums_normalization_and_group_sd():
    counts = np.array([[1, 1, 0], [3, 1, 0], [0, 2, 2], [0, 4, 4], [6, 0, 0], [2, 0, 0]])
    cells = ad.AnnData(sparse.csr_matrix(counts),
                      obs=pd.DataFrame({"cell_type": ["A", "A", "B", "B", "C", "C"],
                                        "donor": ["d1", "d2"] * 3}, index=[f"c{i}" for i in range(6)]),
                      var=pd.DataFrame(index=["g1", "g2", "g3"]))
    pb = module.pseudobulk_counts(cells, target_sum=100)
    expected = np.array([[4, 2, 0], [0, 6, 6], [8, 0, 0]])
    np.testing.assert_array_equal(pb.layers["counts"], expected)
    np.testing.assert_array_equal(pb.layers["counts"].sum(axis=0), counts.sum(axis=0))
    np.testing.assert_allclose(np.expm1(pb.X).sum(axis=1), 100, rtol=2e-7)
    assert pb.obs.n_cells.tolist() == [2, 2, 2]
    assert pb.obs.n_donors.tolist() == [2, 2, 2]
    assert pb.X[0, 2] == 0  # An observed zero remains zero before scaling.
    scaled = module.scale_rna(pd.DataFrame(pb.X.T), max_missing=0)
    np.testing.assert_allclose(scaled.mean(axis=1), 0, atol=1e-14)
    np.testing.assert_allclose(scaled.std(axis=1, ddof=1), 1, atol=1e-14)


def test_refuses_to_sum_logged_or_negative_expression():
    cells = ad.AnnData(np.array([[0.0, 1.2], [0.5, 3.0], [2., 1.]]),
                      obs=pd.DataFrame({"cell_type": ["A", "B", "C"]}, index=["a", "b", "c"]))
    with pytest.raises(ValueError, match="integer counts"):
        module.pseudobulk_counts(cells)


def test_bundled_counts_and_provenance_round_trip():
    pb = ad.read_h5ad(SCRIPT.parent / "data/heart_pseudobulk.h5ad")
    module.validate_prepared(pb)
    assert pb.shape == (6, 4070)
    assert pb.obs.n_cells.sum() == 600
    assert int(pb.layers["counts"].sum()) == 1110049
    assert set(pb.obs_names) == set(pb.obs.cell_type)
