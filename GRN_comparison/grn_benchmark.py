"""Small RNA-only NCI60 GRN comparison using gretapy's public topology tools.

Correlation inference and held-out Ridge scoring here are teaching utilities,
not replacements for gretapy's multimodal inference or benchmark functions.
"""

from dataclasses import asdict, dataclass
from hashlib import sha256
from importlib.metadata import version
import json
from pathlib import Path
import platform
import time

import gretapy as gt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


HERE = Path(__file__).resolve().parent
RNA_PATH = HERE.parent / "data/RNA/RNA_log2_FPKM_cleaned_common.csv"
PRIOR_PATH = HERE / "data/collectri.tsv"


@dataclass(frozen=True)
class Config:
    n_tfs: int = 20
    n_targets: int = 500
    edges_per_tf: int = 10
    test_fraction: float = 0.30
    split_seed: int = 42
    network_seed: int = 42
    alpha: float = 10.0
    max_missing_fraction: float = 0.20

    def validate(self):
        if min(self.n_tfs, self.n_targets, self.edges_per_tf) < 1:
            raise ValueError("TF, target, and edge counts must be positive.")
        if self.n_targets <= self.edges_per_tf:
            raise ValueError("n_targets must exceed edges_per_tf (self-edges are excluded).")
        if not 0 < self.test_fraction < 1 or self.alpha <= 0:
            raise ValueError("test_fraction must be between 0 and 1; alpha must be positive.")
        if not 0 <= self.max_missing_fraction < 1:
            raise ValueError("max_missing_fraction must be in [0, 1).")
        if min(self.split_seed, self.network_seed) < 0:
            raise ValueError("Seeds must be nonnegative.")


def prepare_data(config=Config(), rna_path=RNA_PATH, prior_path=PRIOR_PATH):
    """Split samples first, then choose features and fit predictor medians on train."""
    config.validate()
    rna_path, prior_path = Path(rna_path), Path(prior_path)
    if not rna_path.exists() or not prior_path.exists():
        raise FileNotFoundError(
            "Keep this folder inside COSMOS_basic, including data/collectri.tsv "
            "and ../data/RNA/RNA_log2_FPKM_cleaned_common.csv."
        )
    table = pd.read_csv(rna_path)
    if "gene" not in table or table.gene.isna().any() or table.gene.duplicated().any():
        raise ValueError("RNA input needs a unique, nonmissing 'gene' column.")
    expression = table.set_index("gene").T.apply(pd.to_numeric, errors="raise")
    if expression.index.duplicated().any() or np.isinf(expression.to_numpy()).any():
        raise ValueError("RNA input has duplicate samples or infinite expression values.")
    train_names, test_names = train_test_split(
        expression.index.tolist(), test_size=config.test_fraction, random_state=config.split_seed
    )
    if min(len(train_names), len(test_names)) < 5:
        raise ValueError("The split must contain at least 5 training and 5 test samples.")
    train_raw, test_raw = expression.loc[train_names], expression.loc[test_names]
    variance = train_raw.var()
    valid = (train_raw.isna().mean() <= config.max_missing_fraction) & (variance > 1e-8)
    # Alphabetical order breaks equal-variance ties reproducibly.
    ranked = variance[valid].sort_index().sort_values(ascending=False, kind="stable")
    if len(ranked) < config.n_targets:
        raise ValueError(f"Only {len(ranked)} eligible genes; reduce n_targets.")
    targets = ranked.index[:config.n_targets].tolist()

    prior = pd.read_csv(prior_path, sep="\t")
    if not {"source", "target", "score"}.issubset(prior):
        raise ValueError("Prior needs source, target, and score columns.")
    if prior[["source", "target", "score"]].isna().any().any():
        raise ValueError("Prior contains missing values.")
    prior["score"] = pd.to_numeric(prior.score, errors="raise")
    if not np.isfinite(prior.score).all():
        raise ValueError("Prior scores must be finite.")
    # If the source snapshot has conflicting duplicate signs, their mean is 0.
    prior = prior.groupby(["source", "target"], as_index=False).score.mean()
    reference = prior[
        prior.source.isin(ranked.index) & prior.target.isin(targets)
        & (prior.source != prior.target)
    ].copy()
    counts = reference.groupby("source").size()
    eligible_tfs = counts[counts >= config.edges_per_tf].index
    tfs = ranked[ranked.index.isin(eligible_tfs)].index[:config.n_tfs].tolist()
    if len(tfs) < config.n_tfs:
        raise ValueError(
            f"Only {len(eligible_tfs)} TFs have >= {config.edges_per_tf} prior edges "
            "in the target universe. Reduce n_tfs/edges_per_tf or increase n_targets."
        )
    reference = reference[reference.source.isin(tfs)].reset_index(drop=True)
    genes = sorted(set(tfs) | set(targets))
    medians = train_raw[genes].median()
    train = train_raw[genes].fillna(medians)
    test = test_raw[genes].fillna(medians)
    gene_info = pd.DataFrame({
        "gene": genes,
        "is_tf": [g in tfs for g in genes],
        "is_target": [g in targets for g in genes],
        "train_observed": train_raw[genes].notna().sum().to_numpy(),
        "train_median": medians.to_numpy(),
        "train_variance": variance[genes].to_numpy(),
    })
    return dict(
        train=train, test=test, train_raw=train_raw[genes], test_raw=test_raw[genes],
        tfs=tfs, targets=targets, reference=reference, gene_info=gene_info,
        split=pd.DataFrame({"sample": expression.index,
                            "partition": np.where(expression.index.isin(train_names), "train", "test")}),
        rna_path=rna_path, prior_path=prior_path,
    )


def infer_networks(data, config=Config()):
    """Use identical TFs/targets and exactly k nonself targets per TF in each GRN."""
    networks = {}
    for name, method in [("Pearson", "pearson"), ("Spearman", "spearman")]:
        correlations = data["train"].corr(method=method)
        rows = []
        for tf in data["tfs"]:
            scores = correlations.loc[tf, data["targets"]].drop(labels=tf, errors="ignore")
            scores = scores.dropna().sort_index()
            chosen = scores.abs().sort_values(ascending=False, kind="stable").head(config.edges_per_tf).index
            rows.extend((tf, target, float(scores[target])) for target in chosen)
        networks[name] = pd.DataFrame(rows, columns=["source", "target", "score"])

    prior_seed, random_seed = np.random.SeedSequence(config.network_seed).spawn(2)
    prior_rng, random_rng = np.random.default_rng(prior_seed), np.random.default_rng(random_seed)
    prior_rows, random_rows = [], []
    for tf in data["tfs"]:
        choices = data["reference"].query("source == @tf").sort_values("target")
        selected = choices.iloc[prior_rng.choice(len(choices), config.edges_per_tf, replace=False)]
        prior_rows.append(selected)
        candidates = sorted(set(data["targets"]) - {tf})
        for target in random_rng.choice(candidates, config.edges_per_tf, replace=False):
            random_rows.append((tf, target, 1.0))
    networks["CollecTRI"] = pd.concat(prior_rows, ignore_index=True)
    networks["Random"] = pd.DataFrame(random_rows, columns=["source", "target", "score"])
    for name, network in networks.items():
        if network.duplicated(["source", "target"]).any() or (network.source == network.target).any():
            raise RuntimeError(f"Invalid duplicate/self edge in {name}.")
        if not network.groupby("source").size().eq(config.edges_per_tf).all():
            raise RuntimeError(f"The requested edge budget was not met in {name}.")
        if set(network.source) != set(data["tfs"]):
            raise RuntimeError(f"Not all TFs have valid correlations in {name}.")
    return networks


def prior_agreement(networks, reference):
    """Descriptive reference agreement; unlisted pairs are unknown, not negatives."""
    known = set(zip(reference.source, reference.target))
    rows = []
    for name, network in networks.items():
        edges = set(zip(network.source, network.target))
        shared = len(edges & known)
        rows.append({"name": name, "shared_edges": shared, "reference_edges": len(known),
                     "fraction_in_prior": shared / len(edges),
                     "fraction_of_prior_recovered": shared / len(known),
                     "self_reference": name == "CollecTRI"})
    return pd.DataFrame(rows)


def evaluate_prediction(data, networks, config=Config()):
    """Fit observed train y; score observed test y on the same targets for all GRNs.

    Predictor missing values were filled using train medians. A separate scaler
    is fitted within each target's observed training rows. Edge scores/signs do
    not enter Ridge: this compares graph supports using a common fitted model.
    """
    rows = []
    for name, network in networks.items():
        parents = network.groupby("target").source.apply(list).to_dict()
        for target in data["targets"]:
            observed_train = data["train_raw"][target].notna()
            observed_test = data["test_raw"][target].notna()
            y_train = data["train_raw"].loc[observed_train, target].to_numpy()
            y_test = data["test_raw"].loc[observed_test, target].to_numpy()
            sources = sorted(parents.get(target, []))
            row = dict(name=name, target=target, n_parents=len(sources),
                       n_train_observed=len(y_train), n_test_observed=len(y_test),
                       mse=np.nan, baseline_mse=np.nan, mse_ratio_to_mean=np.nan, r2=np.nan)
            if min(len(y_train), len(y_test)) >= 5:
                baseline = np.full(len(y_test), y_train.mean())
                prediction = baseline.copy()
                if sources:
                    model = make_pipeline(StandardScaler(), Ridge(alpha=config.alpha))
                    model.fit(data["train"].loc[observed_train, sources], y_train)
                    prediction = model.predict(data["test"].loc[observed_test, sources])
                mse = float(np.mean((y_test - prediction) ** 2))
                baseline_mse = float(np.mean((y_test - baseline) ** 2))
                test_variance = float(np.mean((y_test - y_test.mean()) ** 2))
                row.update(mse=mse, baseline_mse=baseline_mse,
                           mse_ratio_to_mean=mse / baseline_mse if baseline_mse > 1e-12 else np.nan,
                           r2=1 - mse / test_variance if test_variance > 1e-12 else np.nan)
            rows.append(row)
    per_target = pd.DataFrame(rows)
    summary = []
    for name, group in per_target.groupby("name", sort=False):
        scored = group[group.mse.notna()]
        summary.append(dict(
            name=name, n_targets=len(group), n_targets_scored=len(scored),
            n_targets_with_parents=int(group.n_parents.gt(0).sum()),
            n_scored_with_parents=int(scored.n_parents.gt(0).sum()),
            n_targets_with_defined_r2=int(scored.r2.notna().sum()),
            n_targets_with_defined_skill=int(scored.mse_ratio_to_mean.notna().sum()),
            median_test_r2=scored.r2.median(),
            mean_skill_vs_train_mean=1 - scored.mse_ratio_to_mean.mean(),
        ))
    return per_target, pd.DataFrame(summary)


def plot_comparison(stats, overlap, agreement, predictions):
    """Small static figure suitable for classroom discussion and export."""
    names = stats.name.tolist()
    colors = ["#3577B3", "#51A898", "#D29D32", "#888888"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
    axes[0, 0].bar(names, stats.n_targets, color=colors)
    axes[0, 0].set(title="Target coverage at equal edge budget", ylabel="Targets with ≥1 parent")
    matrix = pd.DataFrame(np.eye(len(names)), index=names, columns=names)
    for row in overlap.itertuples():
        matrix.loc[row.grn_a, row.grn_b] = row.edge
        matrix.loc[row.grn_b, row.grn_a] = row.edge
    ax = axes[0, 1]
    picture = ax.imshow(matrix, vmin=0, vmax=1, cmap="Blues")
    ax.set(xticks=range(len(names)), yticks=range(len(names)), xticklabels=names,
           yticklabels=names, title="Edge overlap (gretapy)")
    for i in range(len(names)):
        for j in range(len(names)):
            value = matrix.iloc[i, j]
            ax.text(j, i, f"{value:.2f}", ha="center", va="center",
                    color="white" if value > 0.5 else "black")
    fig.colorbar(picture, ax=ax, shrink=0.8)
    agreement = agreement.set_index("name").loc[names]
    axes[1, 0].bar(names, agreement.fraction_in_prior, color=colors)
    axes[1, 0].set(ylim=(0, 1.12), title="Agreement with full eligible CollecTRI",
                   ylabel="Fraction of selected edges listed in prior")
    axes[1, 0].text(2, 1.02, "self-reference", ha="center", fontsize=9)
    predictions = predictions.set_index("name").loc[names]
    axes[1, 1].bar(names, predictions.mean_skill_vs_train_mean, color=colors)
    axes[1, 1].axhline(0, color="black", linewidth=0.8)
    axes[1, 1].set(title="Held-out expression prediction (custom Ridge)",
                   ylabel="Mean skill vs training mean (higher is better)")
    for ax in axes.flat:
        ax.tick_params(axis="x", labelrotation=20)
    fig.suptitle("NCI60: small RNA-only GRN comparison", fontsize=15)
    return fig


def run_benchmark(config=Config(), output_dir=None):
    """Run the practical, save reviewable inputs/results, and return all tables."""
    start = time.perf_counter()
    data = prepare_data(config)
    networks = infer_networks(data, config)
    # These are public gretapy functions for TF–target edge tables, without CREs.
    stats = gt.tl.stats(networks)
    overlap = gt.tl.ocoeff(networks)
    agreement = prior_agreement(networks, data["reference"])
    per_target, summary = evaluate_prediction(data, networks, config)
    fig = plot_comparison(stats, overlap, agreement, summary)
    output = Path(output_dir) if output_dir is not None else HERE / "results"
    output.mkdir(parents=True, exist_ok=True)
    (output / "networks").mkdir(exist_ok=True)
    for name, network in networks.items():
        network.to_csv(output / "networks" / f"{name}.tsv", sep="\t", index=False)
    tables = dict(network_stats=stats, overlap=overlap, prior_agreement=agreement,
                  prediction_per_target=per_target, prediction_summary=summary,
                  split=data["split"], genes=data["gene_info"], collectri_reference=data["reference"])
    for name, table in tables.items():
        table.to_csv(output / f"{name}.tsv", sep="\t", index=False)
    fig.savefig(output / "comparison.png", dpi=160)
    fig.savefig(output / "comparison.pdf")
    info = dict(config=asdict(config), python=platform.python_version(),
                packages={p: version(p) for p in ["gretapy", "mudata", "decoupler", "numpy", "pandas", "scipy", "scikit-learn", "matplotlib"]},
                inputs={key: {"file": str(data[key].relative_to(HERE.parent)),
                              "sha256": sha256(data[key].read_bytes()).hexdigest()}
                        for key in ["rna_path", "prior_path"]},
                n_train=len(data["train"]), n_test=len(data["test"]),
                elapsed_seconds=round(time.perf_counter() - start, 3),
                evaluation="RNA-only topology and custom held-out Ridge; not full GRETA benchmark")
    (output / "run_info.json").write_text(json.dumps(info, indent=2) + "\n")
    return dict(data=data, networks=networks, figure=fig, output_dir=output,
                run_info=info, **tables)
