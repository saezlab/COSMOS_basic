"""Visualize how the TF regulon changes upstream MOON hypotheses.

Run 01_tf_activity.py and 03_moon.R --all first. The main figures use the
controlled comparison: identical unthresholded TF identities, PKN and layers.
Scores are compared only for nodes inferred at the same level in both runs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "cosmos-heart-matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DEFAULT_DIR = HERE / "results/moon_comparison"


def annotate_comparison(frame: pd.DataFrame, cutoff: float = 1.5) -> pd.DataFrame:
    """Distinguish signed activity shifts, support magnitude, and unavailable tests."""
    required = {"cell_type", "source", "score_collectri", "score_scenicplus", "level_collectri", "level_scenicplus"}
    if not required.issubset(frame) or frame.duplicated(["cell_type", "source"]).any():
        raise ValueError("Expected one MOON comparison record per cell type and original node.")
    result = frame.copy()
    c, s = result.score_collectri, result.score_scenicplus
    comparable = np.isfinite(c) & np.isfinite(s) & result.level_collectri.eq(result.level_scenicplus) & result.level_collectri.gt(0)
    if "comparable" in result:
        supplied = result.comparable.astype(str).str.lower().isin(["true", "1"])
        if not supplied.equals(comparable):
            raise ValueError("Saved comparability flags disagree with finite, same-level inferred scores.")
    result["comparable"] = comparable
    result["delta_score"] = (s - c).where(comparable)
    result["delta_abs_score"] = (s.abs() - c.abs()).where(comparable)
    result["sign_flip"] = pd.Series(pd.NA, index=result.index, dtype="boolean")
    result.loc[comparable, "sign_flip"] = np.sign(c[comparable]).ne(np.sign(s[comparable])) & c[comparable].ne(0) & s[comparable].ne(0)
    result["stronger_alignment"] = "not comparable"
    result.loc[comparable, "stronger_alignment"] = np.select(
        [result.loc[comparable, "delta_abs_score"] > 0, result.loc[comparable, "delta_abs_score"] < 0],
        ["SCENIC+", "CollecTRI"], default="equal magnitude")
    result["exceeds_display_cutoff_either"] = comparable & ((c.abs() > cutoff) | (s.abs() > cutoff))
    return result


def collapse_aliases(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Show one representative for nodes that have identical signed children.

    The original-ID table remains untouched. A shared compressed score does not
    distinguish its members, so a plot must not imply independent support for
    every name expanded from that group.
    """
    comparable = frame.loc[frame.comparable].copy()
    ckey, skey = "source_scoring_collectri", "source_scoring_scenicplus"
    if ckey not in comparable or skey not in comparable:
        raise ValueError("Comparison must include both scoring identifiers to audit compressed aliases.")
    if comparable[[ckey, skey]].isna().any().any() or not comparable[ckey].eq(comparable[skey]).all():
        raise ValueError("Controlled runs must share scoring identifiers and compression groups.")
    comparable["scoring_group"] = comparable[ckey]
    records = []
    for group, rows in comparable.groupby("scoring_group", sort=True):
        sources = sorted(rows.source.unique(), key=lambda x: (len(x), x))
        representative = sources[0]
        records.append({"scoring_group": group, "representative": representative,
                        "n_original_nodes": len(sources), "original_nodes": ";".join(sorted(sources))})
        for _, block in rows.groupby("cell_type"):
            for metric in ("score_collectri", "score_scenicplus", "level_collectri", "level_scenicplus"):
                if block[metric].nunique() != 1:
                    raise ValueError("Compressed aliases have inconsistent scores or levels.")
    aliases = pd.DataFrame(records)
    if aliases.empty:
        raise ValueError("No comparable inferred MOON nodes are available.")
    view = comparable.merge(aliases, on="scoring_group", validate="many_to_one")
    view = view.loc[view.source.eq(view.representative)].copy()
    view["plot_label"] = view.representative.str.slice(0, 30)
    mask = view.n_original_nodes > 1
    view.loc[mask, "plot_label"] += " (+" + (view.loc[mask, "n_original_nodes"] - 1).astype(str) + ")"
    return view, aliases


def selected_groups(frame: pd.DataFrame, level: int, n: int) -> list[str]:
    subset = frame.loc[frame.level_collectri.eq(level)]
    eligible = subset.groupby("scoring_group").exceeds_display_cutoff_either.any()
    score = subset.groupby("scoring_group").delta_abs_score.apply(lambda x: x.abs().mean())
    ranking = pd.DataFrame({"mean_absolute_delta": score, "eligible": eligible}).reset_index()
    ranking = ranking.loc[ranking.eligible].sort_values(["mean_absolute_delta", "scoring_group"], ascending=[False, True])
    return ranking.scoring_group.head(n).tolist()


def direct_tf_patterns(frame: pd.DataFrame, pkn: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Collapse direct regulators with the same signed effects on observed TFs.

    Full-PKN compression alone is insufficient: two parents can have different
    unobserved children but the same effective predictor on the measured TFs.
    """
    anchors = set(frame.loc[frame.level_collectri.eq(0) & frame.level_scenicplus.eq(0), "source_scoring_collectri"].dropna())
    direct = frame.loc[frame.comparable & frame.level_collectri.eq(1)].copy()
    edges = pkn.loc[pkn.target.isin(anchors)]
    signatures = {}
    for source, rows in edges.groupby("source", sort=True):
        signature = tuple(sorted((str(row.target), int(row.interaction)) for row in rows.itertuples()))
        signatures[source] = signature
    direct["direct_signature"] = direct.source_scoring_collectri.map(signatures)
    if direct.direct_signature.isna().any():
        raise ValueError("A level-1 regulator has no signed observed-TF footprint in the saved scoring PKN.")
    records, selected = [], []
    for signature, rows in direct.groupby("direct_signature", sort=True):
        original = sorted(rows.source.unique(), key=lambda x: (len(x), x))
        representative = original[0]
        identifier = "L1_" + hashlib.sha256(repr(signature).encode()).hexdigest()[:16]
        record = {"scoring_group": identifier, "representative": representative, "n_original_nodes": len(original),
                  "original_nodes": ";".join(sorted(original)), "n_direct_observed_TFs": len(signature),
                  "direct_signed_TF_targets": ";".join(f"{target}:{sign:+d}" for target, sign in signature),
                  "group_kind": "same effective signed observed-TF footprint"}
        records.append(record)
        for _, block in rows.groupby("cell_type"):
            for metric in ("score_collectri", "score_scenicplus"):
                np.testing.assert_allclose(block[metric], block[metric].iloc[0], rtol=1e-12, atol=1e-12)
        chosen = rows.loc[rows.source.eq(representative)].drop(columns="direct_signature").copy()
        for key, value in record.items():
            chosen[key] = value
        chosen["plot_label"] = representative[:30] + (f" (+{len(original)-1})" if len(original) > 1 else "")
        selected.append(chosen)
    if not selected:
        raise ValueError("No direct upstream comparison is available.")
    return pd.concat(selected, ignore_index=True), pd.DataFrame(records)


def matrix(frame, metric, groups, cells):
    return frame.pivot(index="scoring_group", columns="cell_type", values=metric).reindex(index=groups, columns=cells)


def color_limit(values) -> float:
    finite = np.asarray(values)[np.isfinite(values)]
    return max(float(np.quantile(np.abs(finite), .98)), .1) if finite.size else 1.0


def save_figure(fig, output, name):
    for extension in ("pdf", "png"):
        fig.savefig(output / f"{name}.{extension}", dpi=160, bbox_inches="tight")
    plt.close(fig)


def heatmaps(frame, groups, cells, output):
    labels = frame.drop_duplicates("scoring_group").set_index("scoring_group").loc[groups, "plot_label"]
    c = matrix(frame, "score_collectri", groups, cells)
    s = matrix(frame, "score_scenicplus", groups, cells)
    delta = matrix(frame, "delta_abs_score", groups, cells)
    limit = color_limit(np.concatenate([c.to_numpy().ravel(), s.to_numpy().ravel()]))
    cm = plt.get_cmap("RdBu_r").copy()
    cm.set_bad("#dddddd")
    fig, axes = plt.subplots(1, 3, figsize=(17, max(5, .29 * len(groups) + 2.3)), layout="constrained")
    for ax, values, title, bound, legend in zip(
        axes, (c, s, delta), ("CollecTRI-derived TF activities", "SCENIC+-derived TF activities", "Relative support: |SCENIC+| − |CollecTRI|"),
        (limit, limit, color_limit(delta)), ("MOON ULM t-statistic", "MOON ULM t-statistic", "Difference in absolute MOON scores"), strict=True):
        im = ax.imshow(values, aspect="auto", cmap=cm, vmin=-bound, vmax=bound)
        ax.set_xticks(range(len(cells)), [x.replace(" ", "\n") for x in cells], fontsize=8)
        ax.set_yticks(range(len(groups)), labels, fontsize=8)
        ax.set_title(title, fontsize=10, pad=10)
        fig.colorbar(im, ax=ax, fraction=.035, pad=.02, label=legend)
    fig.suptitle("Direct upstream explanations (level 1): the same TF inputs, network and scoring layers", fontsize=13)
    fig.supxlabel("Delta: red = stronger |score| with SCENIC+; blue = CollecTRI. (+n) = other regulators with the same observed-TF footprint. Gray = unavailable.", fontsize=9)
    save_figure(fig, output, "01_controlled_upstream_heatmaps_level1")


def by_level(frame, cells, output, n=12):
    levels = sorted(frame.level_collectri.astype(int).unique())
    fig, axes = plt.subplots(int(np.ceil(len(levels) / 2)), 2, figsize=(15, 5 * int(np.ceil(len(levels) / 2))), squeeze=False, layout="constrained")
    selections = []
    for ax, level in zip(axes.flat, levels):
        groups = selected_groups(frame, level, n)
        if not groups:
            ax.set_axis_off()
            continue
        values = matrix(frame, "delta_abs_score", groups, cells)
        labels = frame.drop_duplicates("scoring_group").set_index("scoring_group").loc[groups, "plot_label"]
        cm = plt.get_cmap("RdBu_r").copy()
        cm.set_bad("#dddddd")
        bound = color_limit(values)
        im = ax.imshow(values, aspect="auto", cmap=cm, vmin=-bound, vmax=bound)
        ax.set(title=f"Level {level}: difference in absolute MOON scores")
        ax.set_xticks(range(len(cells)), [x.replace(" ", "\n") for x in cells], fontsize=8)
        ax.set_yticks(range(len(groups)), labels, fontsize=8)
        fig.colorbar(im, ax=ax, fraction=.035, pad=.02)
        selections.extend({"level": level, "scoring_group": group, "order": i + 1} for i, group in enumerate(groups))
    for ax in list(axes.flat)[len(levels):]:
        ax.set_axis_off()
    fig.suptitle("Controlled comparison, inspected separately at each upstream layer", fontsize=14)
    fig.supxlabel("Red: larger |score| with SCENIC+; blue: CollecTRI. (+n): indistinguishable members (see tables). Panel color ranges differ; differences are descriptive.", fontsize=9)
    save_figure(fig, output, "02_controlled_support_by_level")
    return pd.DataFrame(selections)


def example_plot(frame, cell_type, output, n_each=5):
    subset = frame.loc[frame.cell_type.eq(cell_type) & frame.level_collectri.eq(1) & frame.exceeds_display_cutoff_either]
    ordered = subset.sort_values(["delta_abs_score", "scoring_group"])
    selected = pd.concat([ordered.loc[ordered.delta_abs_score < 0].head(n_each),
                          ordered.loc[ordered.delta_abs_score > 0].tail(n_each)]).drop_duplicates("scoring_group").sort_values("delta_abs_score")
    if selected.empty:
        return selected
    labels = [f"{row.plot_label} [nTF={int(row.n_direct_observed_TFs)}]{' [sign flip]' if row.sign_flip else ''}" for row in selected.itertuples()]
    fig, ax = plt.subplots(figsize=(11, max(4, .4 * len(selected) + 1.6)), layout="constrained")
    values = selected.delta_abs_score.to_numpy()
    ax.barh(labels, values, color=np.where(values > 0, "#B44445", "#3975A8"))
    ax.axvline(0, color="#666666", linewidth=.8)
    span = max(np.abs(values).max(), .2)
    for i, row in enumerate(selected.itertuples()):
        label = f"C={row.score_collectri:+.2f}; S={row.score_scenicplus:+.2f}"
        ax.text(row.delta_abs_score + (.035 * span if row.delta_abs_score >= 0 else -.035 * span), i, label,
                va="center", ha="left" if row.delta_abs_score >= 0 else "right", fontsize=8)
    ax.set_xlim(min(values.min() - .6 * span, -.2), max(values.max() + .6 * span, .2))
    ax.set(title=f"{cell_type}: illustrative direct upstream differences", xlabel="|SCENIC+ MOON score| − |CollecTRI MOON score|")
    ax.spines[["top", "right"]].set_visible(False)
    fig.supxlabel("C/S: signed scores. nTF: directly observed TF children. (+n): other indistinguishable regulators. Bars show magnitude differences, not significance.", fontsize=9)
    save_figure(fig, output, "03_example_upstream_regulators")
    return selected


def run(input_dir: Path, output: Path, cutoff=1.5, example="Ventricular Cardiomyocyte", n=20):
    input_file = input_dir / "controlled_shared_tf_comparison.csv"
    annotated = annotate_comparison(pd.read_csv(input_file), cutoff)
    view, aliases = collapse_aliases(annotated)
    scoring_files = sorted((input_dir / "controlled/collectri").glob("*/pkn_scoring.csv"))
    if not scoring_files:
        raise FileNotFoundError("Controlled scoring PKN files are needed to audit indistinguishable direct TF footprints.")
    graph_hash = hashlib.sha256(scoring_files[0].read_bytes()).hexdigest()
    if any(hashlib.sha256(path.read_bytes()).hexdigest() != graph_hash for path in scoring_files):
        raise ValueError("Controlled scoring PKNs differ between cell types; use a single consistent comparison.")
    direct, patterns = direct_tf_patterns(annotated, pd.read_csv(scoring_files[0]))
    view = view.loc[view.level_collectri.ne(1)].copy()
    view["group_kind"] = "same full-PKN compression group"
    view = pd.concat([direct, view], ignore_index=True)
    cells = sorted(view.cell_type.unique())
    if example not in cells:
        example = cells[0]
    output.mkdir(parents=True, exist_ok=True)
    annotated.to_csv(output / "upstream_support_all_original_ids.csv", index=False, na_rep="NA")
    view.to_csv(output / "upstream_support_display_groups.csv", index=False, na_rep="NA")
    (output / "upstream_support_unique_scoring_groups.csv").unlink(missing_ok=True)
    aliases.to_csv(output / "indistinguishable_regulators.csv", index=False)
    patterns.to_csv(output / "direct_tf_footprint_groups.csv", index=False)
    groups = selected_groups(view, 1, n)
    if not groups:
        raise ValueError("No level-1 nodes meet the illustrative display cutoff.")
    heatmaps(view, groups, cells, output)
    selection = by_level(view, cells, output)
    selection.to_csv(output / "display_selection_by_level.csv", index=False)
    selected = example_plot(view, example, output)
    selected.to_csv(output / "example_upstream_regulators.csv", index=False)
    summaries = []
    for (cell, level), rows in view.groupby(["cell_type", "level_collectri"]):
        summaries.append({"cell_type": cell, "level": int(level), "n_scoring_groups": len(rows),
                          "pearson": rows.score_collectri.corr(rows.score_scenicplus),
                          "spearman": rows.score_collectri.corr(rows.score_scenicplus, method="spearman"),
                          "n_sign_flips": int(rows.sign_flip.sum()),
                          "n_larger_absolute_scenicplus": int((rows.delta_abs_score > 0).sum()),
                          "n_larger_absolute_collectri": int((rows.delta_abs_score < 0).sum())})
    pd.DataFrame(summaries).to_csv(output / "agreement_by_cell_type_and_level.csv", index=False)
    manifest = {
        "input": input_file.name, "input_sha256": hashlib.sha256(input_file.read_bytes()).hexdigest(),
        "controlled_pkn_sha256": graph_hash,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "comparison": "common unthresholded downstream TFs and identical signaling PKN/layers; only inferred same-level nodes compared",
        "signed_delta": "SCENIC+ minus CollecTRI: change in inferred activity direction/magnitude",
        "support_delta": "absolute SCENIC+ score minus absolute CollecTRI score; larger footprint alignment, not model probability or accuracy",
        "grouping": "level 1: same signed effects on observed TFs, including distinct full-PKN parents; levels 2+: full-PKN compression aliases. Shortest original ID then alphabetical represents each group; all aliases exported",
        "row_selection": "within level, greatest mean absolute support delta across available cell types, requiring |score|>cutoff in at least one run/profile",
        "cutoff": cutoff, "top_n_level1": n, "level1_heatmap_groups": groups, "example_cell_type": example,
        "example_selection": "up to five positive and five negative support deltas at level 1 in the named cell type",
        "limits": "t-statistics are model-dependent; regulon membership, signs, weights and sizes all vary; no donor-level significance or superiority claim",
    }
    (output / "plot_run_info.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Saved controlled upstream comparison plots and annotated tables to {output}")
    print(f"Compared {view.scoring_group.nunique()} display groups ({len(patterns)} effective direct-TF footprints); example: {example}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_DIR / "figures")
    parser.add_argument("--example-cell-type", default="Ventricular Cardiomyocyte")
    parser.add_argument("--display-cutoff", type=float, default=1.5)
    parser.add_argument("--top-n", type=int, default=20)
    args = parser.parse_args()
    if not np.isfinite(args.display_cutoff) or args.display_cutoff < 0 or args.top_n < 1:
        parser.error("Use a nonnegative finite display cutoff and positive top-n.")
    run(args.input_dir, args.output_dir, args.display_cutoff, args.example_cell_type, args.top_n)


if __name__ == "__main__":
    main()
