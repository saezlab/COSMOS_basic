"""Small, data-backed illustrations of the native GRETA Predictive tasks.

These functions explain individual task mechanics. They do not change the
benchmark, its random control, the score aggregation, or its cached outputs.
"""

import decoupler as dc
import gretapy as gt
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.model_selection import train_test_split
from xgboost import XGBRegressor

from diagnose_gene_sets import infer_pathways, summarize
from rna_benchmark import NATIVE_SPLIT_SEED, NATIVE_TEST_FRACTION


METHOD_COLORS = {"Pearson": "#2A6F97", "CollecTRI": "#D97732", "Random": "#777777"}


def _methods(result):
    networks = result["networks"]
    preferred = [name for name in METHOD_COLORS if name in networks]
    return preferred + sorted(set(networks) - set(preferred))


def _color(name):
    return METHOD_COLORS.get(name, "#5F7B59")


def _pathway_label(name):
    return str(name).removeprefix("HALLMARK_").replace("_", " ").capitalize()


def plot_gene_sets_example(result):
    """Show actual Hallmarks hits, network calls, and native precision/recall.

    The figure includes up to 12 pathway examples chosen to span distinct
    expression/network call patterns. All pathways are returned in
    ``tables['pathway_calls']``; the choice is descriptive, not a score filter.
    The fixed native 20,000-gene ORA background is deliberately retained.
    """
    methods = _methods(result)
    adata = result["adata"].copy()
    db = gt.ds.read_db("hg38", "Hallmarks", verbose=False)
    dc.mt.ulm(adata, db, verbose=False)
    scores = adata.obsm["score_ulm"]
    pvals = 2 * stats.t.sf(np.abs(scores.values), adata.n_vars - 2)
    adjusted = stats.false_discovery_control(pvals.ravel(), method="bh").reshape(pvals.shape)
    fraction = ((adjusted < 0.01) & (scores.values > 0)).mean(axis=0)
    fraction = pd.Series(fraction, index=scores.columns)
    hits = set(fraction.index[fraction > 0.2])
    pathway_targets = {source: set(targets) for source, targets in db.groupby("source").target}
    calls = pd.DataFrame({"pathway": sorted(pathway_targets)})
    calls["label"] = calls.pathway.map(_pathway_label)
    calls["positive_rna_fraction"] = calls.pathway.map(fraction).fillna(0.0)
    calls["expression_hit"] = calls.pathway.isin(hits)
    rows = []
    for name in methods:
        predicted = infer_pathways(result["networks"][name], pathway_targets, 20000, exclude_self=False)
        calls[name] = calls.pathway.isin(predicted)
        rows.append({"name": name, **summarize(predicted, hits)})
    counts = pd.DataFrame(rows)
    counts["counted_unit"] = "Hallmark pathways"

    # Verify that this illustration is the native result, when it was evaluated.
    saved = result.get("results")
    if saved is not None:
        for row in rows:
            match = saved.loc[(saved.name == row["name"]) & (saved.db == "Hallmarks")]
            if len(match) == 1:
                for metric in ("precision", "recall", "f01"):
                    np.testing.assert_allclose(row[metric], match.iloc[0][metric], rtol=1e-10, atol=1e-12)

    # One alphabetically first example per binary call pattern, then a second
    # pass if fewer than 12 patterns occur. Selection never changes the scores.
    pattern_columns = ["expression_hit", *methods]
    ordered = calls.sort_values([*pattern_columns, "pathway"], ascending=[False] * len(pattern_columns) + [True])
    groups = list(ordered.groupby(pattern_columns, sort=False))
    chosen = []
    for position in range(len(calls)):
        for _, group in groups:
            if len(group) > position:
                chosen.append(group.index[position])
                if len(chosen) == min(12, len(calls)):
                    break
        if len(chosen) == min(12, len(calls)):
            break
    examples = calls.loc[chosen].copy()
    calls["shown_in_figure"] = calls.index.isin(chosen)

    figure = plt.figure(figsize=(11.5, 8.3), layout="constrained")
    grid = figure.add_gridspec(2, 2, height_ratios=[1, 2.3], width_ratios=[1.2, 1])
    count_ax = figure.add_subplot(grid[0, 0])
    y = np.arange(len(methods))
    count_ax.barh(y, counts.tp, color="#2A6F97", label="Also an expression hit (TP)")
    count_ax.barh(y, counts.fp, left=counts.tp, color="#C9D4DC", label="Other pathway calls (FP)")
    for i, row in counts.iterrows():
        count_ax.text(row.n_predicted + 0.45, i, f"{row.tp}/{row.n_expression_hits} RNA hits recovered", va="center", fontsize=9)
    count_ax.set(yticks=y, yticklabels=methods, xlabel="Pathways enriched in any TF regulon",
                 title="What the network predicts", xlim=(0, max(counts.n_predicted.max() + 19, 25)))
    count_ax.invert_yaxis()
    count_ax.spines[["top", "right"]].set_visible(False)
    count_ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(0, -0.30))

    metric_ax = figure.add_subplot(grid[0, 1])
    metric_ax.axis("off")
    table_values = [[row["name"], f"{row['precision']:.3f}", f"{row['recall']:.3f}", f"{row['f01']:.3f}"] for row in rows]
    metric_table = metric_ax.table(cellText=table_values, colLabels=["Network", "Precision", "Recall", "F₀.₁"],
                                   cellLoc="center", colWidths=[0.32, 0.23, 0.22, 0.22], loc="center")
    metric_table.auto_set_font_size(False)
    metric_table.set_fontsize(10)
    metric_table.scale(1, 1.6)
    for i, name in enumerate(methods, start=1):
        metric_table[i, 0].get_text().set_color(_color(name))
    metric_ax.set_title("Native score: a strong preference for precision", fontsize=11)

    example_ax = figure.add_subplot(grid[1, :])
    matrix = examples[pattern_columns].to_numpy(dtype=int)
    example_ax.imshow(matrix, cmap=ListedColormap(["#F1F3F5", "#2A6F97"]), vmin=0, vmax=1,
                      aspect="auto", interpolation="nearest")
    labels = [f"{row.label}  ({row.positive_rna_fraction:.0%} of lines)" for row in examples.itertuples()]
    example_ax.set(xticks=np.arange(len(pattern_columns)),
                   xticklabels=["Expression hit", *methods], yticks=np.arange(len(examples)), yticklabels=labels)
    example_ax.tick_params(axis="y", labelsize=9, length=0)
    example_ax.tick_params(axis="x", labelsize=10, length=0)
    example_ax.xaxis.tick_top()
    example_ax.set_xticks(np.arange(-0.5, len(pattern_columns), 1), minor=True)
    example_ax.set_yticks(np.arange(-0.5, len(examples), 1), minor=True)
    example_ax.grid(which="minor", color="white", linewidth=2)
    example_ax.tick_params(which="minor", bottom=False, left=False)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            example_ax.text(j, i, "Yes" if matrix[i, j] else "–", ha="center", va="center",
                            color="white" if matrix[i, j] else "#64717B", fontsize=9)
    example_ax.set_title("Hallmarks examples: RNA enrichment versus pooled target-set enrichment", pad=31, fontsize=11)
    figure.suptitle("Gene Sets: compare RNA-enriched pathways with labels in TF target sets", fontsize=14)
    figure.supxlabel("Expression hit: positive ULM enrichment at q < 0.01 in >20% of all cell lines.\n"
                    "Network call: ORA q < 0.01 in any regulon, using native 20,000-gene background. This does not validate individual edges.",
                    fontsize=9)
    metadata = pd.DataFrame([dict(database="Hallmarks", n_cell_lines=adata.n_obs,
                                  n_rna_genes=adata.n_vars, n_expression_hits=len(hits),
                                  ora_background=20000,
                                  example_selection="First pathway alphabetically per distinct binary call pattern; repeated to show at most 12.",
                                  interpretation="Descriptive pathway concordance; fixed native ORA background remains uncalibrated.")])
    return figure, {"pathway_calls": calls, "summary": counts, "examples": examples, "metadata": metadata}


def _native_predictors(network):
    """Match native uniqueness, absolute-score sorting and top-five selection."""
    edges = network.drop_duplicates(["source", "target"]).copy()
    edges["abs_score"] = edges.score.abs()
    edges = edges.sort_values("abs_score", ascending=False)
    return {target: group.head(5).loc[lambda x: x.source != target].copy()
            for target, group in edges.groupby("target")}


def plot_omics_example(result, target=None):
    """Refit the native model for one target and display held-out predictions.

    Default target selection uses *training RNA only*: among targets with at
    least one incoming TF in every displayed network, select the greatest
    absolute training TF–target Pearson correlation in the Pearson network
    (or first network if Pearson is unavailable), breaking ties by gene name.
    A selected target must still meet native test-size eligibility. The raw
    per-gene Spearman statistics are illustrations, not BH-adjusted discoveries.
    """
    adata = result["adata"]
    methods = _methods(result)
    predictors = {name: _native_predictors(result["networks"][name]) for name in methods}
    train, test = train_test_split(list(adata.obs_names), test_size=NATIVE_TEST_FRACTION,
                                   random_state=NATIVE_SPLIT_SEED, stratify=adata.obs["celltype"])
    expected = result["split"]
    if set(train) != set(expected.loc[expected.partition == "train", "sample"]):
        raise ValueError("The example requires the same training split as the RNA benchmark.")
    common = set.intersection(*(set(table) for table in predictors.values()))
    common = sorted(gene for gene in common if all(len(predictors[name][gene]) for name in methods))
    selection_network = "Pearson" if "Pearson" in methods else methods[0]
    train_values = pd.DataFrame(np.asarray(adata[train].X), index=train, columns=adata.var_names)
    candidates = []
    for gene in common:
        y = train_values[gene].to_numpy()
        nonzero = y != 0.0
        if nonzero.sum() < 10 or np.ptp(y[nonzero]) == 0:
            continue
        tf_names = predictors[selection_network][gene].source.tolist()
        correlations = [np.corrcoef(train_values.loc[nonzero, tf], y[nonzero])[0, 1]
                        for tf in tf_names if np.ptp(train_values.loc[nonzero, tf]) > 0]
        finite = np.asarray(correlations)[np.isfinite(correlations)]
        if len(finite):
            candidates.append(dict(target=gene, max_abs_training_correlation=float(np.abs(finite).max()),
                                   n_nonzero_training=int(nonzero.sum())))
    if not candidates:
        raise ValueError("No common target with variable, nonzero training RNA and eligible predictors.")
    candidates = pd.DataFrame(candidates).sort_values(
        ["max_abs_training_correlation", "target"], ascending=[False, True], kind="stable")
    automatically_selected = target is None
    target = str(candidates.iloc[0].target) if automatically_selected else str(target)
    if target not in common:
        raise ValueError(f"{target}: choose a target with incoming TFs in every displayed network.")
    train_y = np.asarray(adata[train, target].X).ravel()
    test_y = np.asarray(adata[test, target].X).ravel()
    train_keep, test_keep = train_y != 0.0, test_y != 0.0
    if test_keep.sum() < 10:
        raise ValueError(f"{target}: fewer than 10 nonzero held-out responses; native Omics would omit this target.")
    if train_keep.sum() < 2 or np.ptp(train_y[train_keep]) == 0:
        raise ValueError(f"{target}: training response must be variable.")
    test_names = np.asarray(test)[test_keep]
    true_y = test_y[test_keep]
    prediction_rows, score_rows, edge_rows = [], [], []
    for name in methods:
        edges = predictors[name][target]
        sources = edges.source.tolist()
        train_x = np.asarray(adata[train, sources].X)[train_keep]
        test_x = np.asarray(adata[test, sources].X)[test_keep]
        model = XGBRegressor(random_state=0, n_jobs=1).fit(train_x, train_y[train_keep])
        predicted = model.predict(test_x)
        if np.ptp(predicted) > 0 and np.ptp(true_y) > 0:
            coefficient, pvalue = stats.spearmanr(predicted, true_y)
        else:
            coefficient, pvalue = np.nan, np.nan
        score_rows.append(dict(name=name, target=target, n_train=int(train_keep.sum()), n_test=int(test_keep.sum()),
                               n_predictors=len(sources), spearman_rho=coefficient, raw_pvalue=pvalue))
        for sample, observed, estimated in zip(test_names, true_y, predicted):
            prediction_rows.append(dict(name=name, target=target, sample=sample,
                                        tissue=str(adata.obs.loc[sample, "celltype"]),
                                        observed=float(observed), predicted=float(estimated)))
        for rank, edge in enumerate(edges.itertuples(), start=1):
            edge_rows.append(dict(name=name, target=target, rank=rank, source=edge.source, score=float(edge.score)))
    predictions = pd.DataFrame(prediction_rows)
    metrics = pd.DataFrame(score_rows)
    edges = pd.DataFrame(edge_rows)
    tissues = sorted(predictions.tissue.unique())
    tissue_colors = {tissue: plt.get_cmap("tab10")(i % 10) for i, tissue in enumerate(tissues)}
    figure = plt.figure(figsize=(11.5, 6.3), layout="constrained")
    grid = figure.add_gridspec(3, len(methods), height_ratios=[3.7, 1.2, 0.5])
    all_values = predictions[["observed", "predicted"]].to_numpy()
    lo, hi = float(all_values.min()), float(all_values.max())
    padding = max((hi - lo) * 0.07, 0.1)
    limits = (lo - padding, hi + padding)
    for column, name in enumerate(methods):
        ax = figure.add_subplot(grid[0, column])
        frame = predictions.loc[predictions.name == name]
        for tissue in tissues:
            selected = frame.loc[frame.tissue == tissue]
            ax.scatter(selected.observed, selected.predicted, s=43, color=tissue_colors[tissue],
                       edgecolor="white", linewidth=0.5, alpha=0.9)
        ax.plot(limits, limits, linestyle="--", color="#A5ADB4", linewidth=1, zorder=0)
        metric = metrics.loc[metrics.name == name].iloc[0]
        rho_text = f"ρ = {metric.spearman_rho:.2f}" if np.isfinite(metric.spearman_rho) else "ρ undefined (constant values)"
        ax.set(title=f"{name}\n{rho_text}; n = {metric.n_test}", xlabel="Observed RNA", xlim=limits, ylim=limits)
        ax.set_aspect("equal", adjustable="box")
        ax.title.set_color(_color(name))
        if column == 0:
            ax.set_ylabel("Predicted RNA")
        ax.spines[["top", "right"]].set_visible(False)
        table_ax = figure.add_subplot(grid[1, column])
        table_ax.axis("off")
        method_edges = edges.loc[edges.name == name]
        lines = [f"{row.source}: {row.score:+.3f}" for row in method_edges.itertuples()]
        table_ax.text(0.03, 1, "Actual incoming TFs (edge weight)\n" + "\n".join(lines),
                      va="top", fontsize=8.5, wrap=True, transform=table_ax.transAxes)
    legend = [Patch(facecolor=tissue_colors[tissue], label=tissue) for tissue in tissues]
    legend_ax = figure.add_subplot(grid[2, :])
    legend_ax.axis("off")
    legend_ax.legend(handles=legend, title="Tissue of origin", loc="center", ncol=min(9, len(tissues)),
                     frameon=False, fontsize=8, title_fontsize=9)
    figure.suptitle(f"Omics example: can incoming TF expression predict {target} in held-out cell lines?", fontsize=13)
    figure.supxlabel("RNA axes: log₂(FPKM + 1). Dashed line: perfect prediction. Raw held-out Spearman ρ; no per-gene BH claim.\n"
                    "Target chosen from training RNA. XGBoost uses selected TF expression as predictors; zero responses are excluded.", fontsize=9)
    rule = (f"Greatest absolute training TF–target correlation among common eligible targets using {selection_network} predictors; "
            "gene-name tie break; no held-out predictions or correlations used.") if automatically_selected else "User-specified target; no automatic target selection."
    selection = pd.DataFrame([dict(target=target, automatically_selected=automatically_selected,
                                   selection_rule=rule, selection_network=selection_network,
                                   n_candidate_targets=len(candidates), n_train_before_zero_filter=len(train),
                                   n_test_before_zero_filter=len(test), model="XGBRegressor(random_state=0, n_jobs=1)",
                                   note="Edge weights select up to five incoming TFs, but are not regression coefficients. Single-gene raw Spearman is illustrative; native BH adjustment uses all tested genes.")])
    return figure, {"predictions": predictions, "metrics": metrics, "predictors": edges,
                    "selection": selection, "training_candidates": candidates}
