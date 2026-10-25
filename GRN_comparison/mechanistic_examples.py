"""Show how the native KnockTF activity task uses a real perturbation.

The example is selected without looking at activity scores. All scores and
adjusted p-values are recomputed with the same experiment eligibility, decoupler
ULM calls, and across-experiment BH family as the benchmark's native TF task.
These helpers do not download a reference when the task was not requested.
"""

import textwrap

import decoupler as dc
import gretapy as gt
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import false_discovery_control

from greta_benchmark import _knocktf_support, configure_cache


COLORS = {"Pearson": "#2A6F97", "CollecTRI": "#D97732", "Random": "#777777"}


def _empty_example(message):
    fig, ax = plt.subplots(figsize=(11, 2.2))
    ax.axis("off")
    ax.text(0.5, 0.58, "KnockTF perturbation example", ha="center", va="center",
            fontsize=14, weight="bold", transform=ax.transAxes)
    ax.text(0.5, 0.30, textwrap.fill(message, 110), ha="center", va="center",
            fontsize=11, transform=ax.transAxes)
    tables = {key: pd.DataFrame() for key in
              ("experiment", "activities", "selected_activity", "targets", "gene_responses", "coverage")}
    return fig, tables


def _activity_table(result, db):
    """Reproduce native per-experiment ULM and per-network BH correction.

    The >=3 nonzero measured-target gate is the compatibility correction used
    by ``evaluate_networks``. The native ULM itself fits all nonzero-response
    reference genes, assigning zero weight to genes outside the regulon.
    """
    terms = result.get("run_info", {}).get("manifest", {}).get("terms", {})
    cats = (terms or {}).get("KnockTF (scoring)")
    activities, coverage = [], []
    for name, graph in result["networks"].items():
        graph = graph.loc[graph.target.isin(db.var_names)]
        graph = graph.groupby(["source", "target"], as_index=False).score.mean()
        eligible, supported, insufficient = _knocktf_support(graph, result["adata"], db, cats)
        rows = []
        for experiment in supported:
            tf = db.obs.at[experiment, "TF"]
            source_grn = graph.loc[graph.source == tf].rename(columns={"score": "weight"})
            # This is the call made by gretapy.tl._mechanistic._tfa.
            scores, pvals = dc.mt.ulm(data=db[[experiment], :].to_df(), net=source_grn, tmin=3)
            rows.append(dict(name=name, experiment=experiment, TF=tf,
                             activity=float(scores.iloc[0, 0]), pval=float(pvals.iloc[0, 0])))
        frame = pd.DataFrame(rows, columns=["name", "experiment", "TF", "activity", "pval"])
        frame["padj"] = false_discovery_control(frame.pval.to_numpy(), method="bh") if len(frame) else []
        frame["hit"] = (frame.activity < 0) & (frame.padj < 0.05)
        activities.append(frame)
        coverage.append(dict(name=name, n_eligible=len(eligible), n_tested=len(supported),
                             n_hits=int(frame.hit.sum()),
                             n_skipped_insufficient_targets=insufficient))
    return pd.concat(activities, ignore_index=True), pd.DataFrame(coverage)


def plot_knocktf_example(result, experiment=None):
    """Return a real KnockTF example figure and the underlying data tables.

    Parameters
    ----------
    result : dict
        The output of ``run_rna`` with the teaching profile. The quick profile
        returns an explanatory figure without accessing the KnockTF reference.
    experiment : str, optional
        A reference experiment supported by every compared network. By default,
        choose the upper median own-TF log fold change among common supported
        experiments, ordered by logFC then experiment ID. Activity is never used
        for selection. The displayed example therefore need not be a hit.

    Returns
    -------
    fig : matplotlib.figure.Figure
    tables : dict[str, pandas.DataFrame]
        ``experiment`` metadata, all ``activities`` with native adjusted p-values,
        ``selected_activity``, all selected-TF ``targets`` (including excluded zero
        responses), the full experiment's ``gene_responses``, and whole-task
        ``coverage``. No file is saved by this helper.
    """
    evaluated = result.get("results", pd.DataFrame())
    if "db" not in evaluated or "KnockTF (scoring)" not in set(evaluated.db):
        return _empty_example(
            "KnockTF was not included in this run. Use the teaching profile for the external "
            "perturbation example; the quick profile does not download or evaluate KnockTF."
        )
    if not result.get("networks"):
        return _empty_example("No networks are available for a perturbation example.")

    configure_cache()
    db = gt.ds.read_db("hg38", "KnockTF (scoring)", verbose=False)
    activities, coverage = _activity_table(result, db)
    names = list(result["networks"])
    supported_sets = [set(activities.loc[activities.name == name, "experiment"]) for name in names]
    common = set.intersection(*supported_sets)
    if not common:
        fig, tables = _empty_example(
            "No KnockTF experiment has at least three nonzero measured targets in every compared network. "
            "This is insufficient common coverage, rather than evidence of no regulatory effect."
        )
        tables.update(activities=activities, coverage=coverage)
        return fig, tables
    candidates = db.obs.loc[sorted(common), ["logFC"]].copy()
    candidates["experiment"] = candidates.index.astype(str)
    candidates = candidates.sort_values(["logFC", "experiment"], kind="stable")
    if experiment is None:
        experiment = candidates.index[len(candidates) // 2]
        selection = "Upper median own-TF logFC among experiments supported by every network; activity-blind"
    elif experiment not in common:
        raise ValueError("Choose an experiment supported by every network; available IDs are in the activities table.")
    else:
        selection = "User-selected reference experiment"
    meta = db.obs.loc[[experiment]].copy()
    meta.insert(0, "experiment", experiment)
    meta["selection_rule"] = selection
    meta["n_common_supported_experiments"] = len(common)
    tf = str(meta.TF.iloc[0])
    response = db[[experiment], :].X
    response = response.toarray() if sparse.issparse(response) else np.asarray(response)
    response = pd.Series(response.ravel(), index=db.var_names, name="logFC")
    meta["n_nonzero_response_features_used_by_ulm"] = int((response != 0).sum())
    selected = activities.loc[activities.experiment == experiment].copy()
    targets = []
    for name, graph in result["networks"].items():
        frame = graph.loc[graph.source == tf, ["source", "target", "score"]]
        frame = frame.groupby(["source", "target"], as_index=False).score.mean()
        frame["name"] = name
        frame["logFC"] = frame.target.map(response)
        frame["measured_in_reference"] = frame.logFC.notna()
        frame["used_by_ulm"] = frame.measured_in_reference & (frame.logFC != 0)
        frame["signed_response"] = np.sign(frame.score) * frame.logFC
        targets.append(frame)
    targets = pd.concat(targets, ignore_index=True)

    fig = plt.figure(figsize=(max(11.5, 4.2 * len(names)), 7.0))
    grid = fig.add_gridspec(2, len(names), height_ratios=[2.8, 1.4],
                           left=0.07, right=0.98, bottom=0.17, top=0.80,
                           hspace=0.60, wspace=0.30)
    nonzero = response.loc[response != 0]
    ylo, yhi = float(nonzero.min()), float(nonzero.max())
    padding = max((yhi - ylo) * 0.08, 0.1)
    for i, name in enumerate(names):
        ax = fig.add_subplot(grid[0, i])
        color = COLORS.get(name, f"C{i % 10}")
        frame = targets.loc[(targets.name == name) & targets.used_by_ulm].copy()
        weights = pd.Series(0.0, index=nonzero.index)
        weights.loc[frame.target] = frame.score.to_numpy()
        background = nonzero.loc[~nonzero.index.isin(frame.target)]
        ax.scatter(np.zeros(len(background)), background.to_numpy(), color="#b4bac0",
                   s=5, alpha=0.08, linewidth=0, rasterized=True, zorder=1)
        for positive, marker in ((True, "o"), (False, "s")):
            group = frame.loc[(frame.score >= 0) == positive]
            ax.scatter(group.score, group.logFC, s=48, marker=marker, color=color,
                       edgecolor="white", linewidth=0.7, zorder=3)
        # Draw the actual one-covariate fitted relationship, including the
        # non-regulon genes assigned zero weight by native decoupler ULM.
        x = weights.to_numpy()
        y = nonzero.to_numpy()
        slope = np.dot(x - x.mean(), y - y.mean()) / np.dot(x - x.mean(), x - x.mean())
        intercept = y.mean() - slope * x.mean()
        limits = np.array([-1.08, 1.08]) * max(1.0, float(np.abs(x).max()))
        ax.plot(limits, intercept + slope * limits, color=color, linewidth=1.3, alpha=0.8)
        # A short, deterministic label list: alphabetic target names, never
        # selected by fold change or agreement with the expected response.
        labelled = frame.sort_values("target").head(5)
        for positive in (True, False):
            group = labelled.loc[(labelled.score >= 0) == positive].sort_values("logFC")
            if group.empty:
                continue
            positions_y = group.logFC.to_numpy().copy()
            gap = 0.065 * (yhi - ylo + 2 * padding)
            for j in range(1, len(positions_y)):
                positions_y[j] = max(positions_y[j], positions_y[j - 1] + gap)
            positions_y += float(group.logFC.mean()) - positions_y.mean()
            positions_y += max(0, ylo - positions_y.min())
            positions_y -= max(0, positions_y.max() - yhi)
            for row, label_y in zip(group.itertuples(), positions_y):
                dx = -0.05 if positive else 0.05
                ax.annotate(row.target, (row.score, row.logFC),
                            xytext=(row.score + dx, label_y), textcoords="data",
                            ha="right" if positive else "left", va="center",
                            fontsize=8, color=color,
                            arrowprops={"arrowstyle": "-", "linewidth": 0.5,
                                        "color": color, "alpha": 0.65})
        ax.axhline(0, color="#666666", linewidth=0.7, linestyle=":")
        ax.set(xlim=limits, ylim=(ylo - padding, yhi + padding),
               xlabel="Signed TF–target weight", title=f"{name}: {len(frame)} usable targets")
        if i == 0:
            ax.set_ylabel("Observed gene log fold change\n(external TF knockdown)")
        ax.spines[["top", "right"]].set_visible(False)
    legend = [Line2D([], [], color="#777777", marker="o", linestyle="", label="Positive edge"),
              Line2D([], [], color="#777777", marker="s", linestyle="", label="Negative edge"),
              Line2D([], [], color="#b4bac0", marker="o", linestyle="", markersize=4,
                     label="Other measured genes: weight 0")]
    fig.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.5, 0.88),
               ncol=3, frameon=False, fontsize=9)

    ax = fig.add_subplot(grid[1, :])
    selected = selected.set_index("name").loc[names].reset_index()
    positions = np.arange(len(names))
    ax.barh(positions, selected.activity, color=[COLORS.get(n, f"C{i % 10}") for i, n in enumerate(names)],
            height=0.62)
    ax.axvline(0, color="#444444", linewidth=0.8)
    span = max(float(selected.activity.abs().max()), 1.0)
    ax.set_xlim(min(-0.2 * span, float(selected.activity.min()) - 0.1 * span),
                max(0.2 * span, float(selected.activity.max()) + 0.1 * span) + 0.7 * span)
    text_x = max(0, float(selected.activity.max())) + 0.08 * span
    for i, row in enumerate(selected.itertuples()):
        label = "hit" if row.hit else "not a hit"
        ax.text(text_x, i, f"activity {row.activity:.2f}; q = {row.padj:.3g} ({label})",
                va="center", fontsize=9)
    ax.set(yticks=positions, yticklabels=names,
           xlabel="Native ULM activity (negative = reduced inferred TF activity)")
    ax.invert_yaxis()
    ax.spines[["top", "right"]].set_visible(False)
    specimen = str(meta["Biosample.Name"].iloc[0])
    tissue = str(meta["Tissue.Type"].iloc[0]).replace("_", " ")
    own_lfc = float(meta.logFC.iloc[0])
    fig.suptitle(f"Real external KnockTF experiment: {tf} knockdown in {specimen}\n"
                 f"{experiment} · {tissue} · measured TF logFC = {own_lfc:.2f}",
                 y=0.98, fontsize=13, weight="bold")
    footer = ("Under the signed-edge model, activation edges predict decreases after knockdown; repression edges predict increases. "
              "A hit requires activity < 0 and BH q < 0.05. The q-values use all supported experiments "
              "within each network, not only this example. This is an external perturbation, not an NCI60 intervention.")
    fig.text(0.07, 0.035, textwrap.fill(footer, width=155), ha="left", va="bottom", fontsize=9)
    gene_responses = response.rename_axis("gene").reset_index()
    gene_responses["used_by_ulm"] = gene_responses.logFC != 0
    return fig, dict(experiment=meta.reset_index(drop=True), activities=activities,
                     selected_activity=selected, targets=targets,
                     gene_responses=gene_responses, coverage=coverage)
