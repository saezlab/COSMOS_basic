"""Guard descriptive support differences and compression-aware illustrations."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location("moon_comparison", Path(__file__).resolve().parents[1] / "04_compare_moon.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_negative_activity_support_and_incomparability():
    frame = pd.DataFrame({"cell_type": ["A"] * 4, "source": ["X", "Y", "TF", "missing"],
                          "score_collectri": [-1., 2., 3., np.nan], "score_scenicplus": [-3., -2., 5., 4.],
                          "level_collectri": [1, 1, 0, np.nan], "level_scenicplus": [1, 1, 0, 1]})
    result = module.annotate_comparison(frame).set_index("source")
    assert result.loc["X", "delta_score"] == -2
    assert result.loc["X", "delta_abs_score"] == 2
    assert result.loc["X", "stronger_alignment"] == "SCENIC+"
    assert not result.loc["X", "sign_flip"]
    assert result.loc["Y", "sign_flip"]
    assert result.loc["Y", "delta_abs_score"] == 0
    assert result.loc[["TF", "missing"], "delta_abs_score"].isna().all()


def test_aliases_are_not_independent_regulator_evidence():
    frame = pd.DataFrame({"cell_type": ["A"] * 3, "source": ["X", "XX", "Y"],
                          "score_collectri": [-1., -1., 2.], "score_scenicplus": [-3., -3., 1.],
                          "level_collectri": [1] * 3, "level_scenicplus": [1] * 3,
                          "source_scoring_collectri": ["group1", "group1", "Y"],
                          "source_scoring_scenicplus": ["group1", "group1", "Y"]})
    annotated = module.annotate_comparison(frame)
    reduced, aliases = module.collapse_aliases(annotated)
    assert set(reduced.source) == {"X", "Y"}
    assert aliases.set_index("scoring_group").loc["group1", "n_original_nodes"] == 2
    assert len(annotated) == 3
    assert module.selected_groups(reduced, 1, 20) == ["group1", "Y"]


def test_distinct_pkn_parents_can_share_the_same_observed_tf_footprint():
    frame = pd.DataFrame({"cell_type": ["A"] * 4, "source": ["HDAC4", "HDAC5", "TF1", "TF2"],
                          "score_collectri": [-3., -3., 1., 2.], "score_scenicplus": [-.5, -.5, 3., 4.],
                          "level_collectri": [1, 1, 0, 0], "level_scenicplus": [1, 1, 0, 0],
                          "source_scoring_collectri": ["HDAC4", "HDAC5", "TF1", "TF2"],
                          "source_scoring_scenicplus": ["HDAC4", "HDAC5", "TF1", "TF2"]})
    pkn = pd.DataFrame({"source": ["HDAC4", "HDAC4", "HDAC4", "HDAC5", "HDAC5", "HDAC5"],
                        "target": ["TF1", "TF2", "unmeasuredA", "TF1", "TF2", "unmeasuredB"],
                        "interaction": [-1] * 6})
    view, groups = module.direct_tf_patterns(module.annotate_comparison(frame), pkn)
    assert len(view) == len(groups) == 1
    assert view.iloc[0].plot_label == "HDAC4 (+1)"
    assert view.iloc[0].n_direct_observed_TFs == 2
    assert groups.iloc[0].original_nodes == "HDAC4;HDAC5"
