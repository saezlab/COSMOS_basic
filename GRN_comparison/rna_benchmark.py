"""Full expressed-gene, RNA-only NCI60 benchmark with native GRETA tasks.

Correlation inference is a deliberately simple local baseline. Evaluation and
ranking use gretapy; no ATAC measurements or regulatory regions are invented.
"""

from dataclasses import asdict, dataclass
from hashlib import sha256
from importlib.metadata import version
import json
from pathlib import Path
import platform
import time

import anndata as ad
import gretapy as gt
import gretapy.config as gt_config
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata
from sklearn.model_selection import train_test_split

from greta_benchmark import comparable_results, configure_cache, evaluate_networks


HERE = Path(__file__).resolve().parent
RNA_PATH = HERE / "data/nci60_rna_log2fpkm.tsv.gz"
SAMPLES_PATH = HERE / "data/nci60_samples.tsv"
RNA_METRICS = (
    "Hallmarks", "PROGENy", "Reactome", "Gene ~ TFs", "HPA", "TF-Marker",
    "IntAct", "Europe PMC", "KnockTF (scoring)",
)
RNA_DIAGNOSTICS = ("CollecTRI",)
NATIVE_SPLIT_SEED = 42
NATIVE_TEST_FRACTION = 0.33
PACKAGES = ("gretapy", "decoupler", "anndata", "mudata", "numpy", "pandas",
            "scipy", "scikit-learn", "xgboost", "matplotlib")


@dataclass(frozen=True)
class RnaConfig:
    seed: int = 42
    min_expression: float = 1.0
    min_sample_fraction: float = 0.8
    max_genes: int = 0
    max_tfs: int = 0
    include_spearman: bool = False
    profile: str = "teaching"

    def validate(self):
        if self.profile not in {"quick", "teaching"}:
            raise ValueError("profile must be quick or teaching.")
        if self.seed < 0 or self.max_genes < 0 or self.max_tfs < 0:
            raise ValueError("seed and optional gene/TF limits must be nonnegative.")
        if not np.isfinite(self.min_expression) or self.min_expression < 0:
            raise ValueError("min_expression must be finite and nonnegative.")
        if not 0 < self.min_sample_fraction <= 1:
            raise ValueError("min_sample_fraction must be in (0, 1].")


def _file_hash(path):
    digest = sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _input_record(path):
    path = Path(path).resolve()
    return dict(path=str(path), sha256=_file_hash(path))


def load_rna(config=RnaConfig(), rna_path=RNA_PATH, samples_path=SAMPLES_PATH):
    """Read nonnegative RNA, split first, and filter genes using training RNA only.

    TSV input has a ``gene`` column and one column per sample. Matching sample
    metadata requires ``sample`` and ``tissue_code``. Prepared AnnData can instead
    carry ``obs['celltype']`` or ``obs['group']``; labels must permit stratification.
    """
    config.validate()
    rna_path = Path(rna_path).resolve()
    if not rna_path.exists():
        raise FileNotFoundError(f"RNA input does not exist: {rna_path}")
    samples_used = None
    if rna_path.suffix.lower() == ".h5ad":
        adata = ad.read_h5ad(rna_path)
        if "celltype" not in adata.obs:
            if "group" not in adata.obs:
                raise ValueError("AnnData requires obs['celltype'] or obs['group'] with real sample groups.")
            adata.obs["celltype"] = adata.obs["group"].copy()
    else:
        table = pd.read_csv(rna_path, sep="\t")
        if "gene" not in table or table.gene.isna().any() or table.gene.duplicated().any():
            raise ValueError("RNA TSV requires a unique, nonmissing 'gene' column.")
        matrix = table.set_index("gene").T.apply(pd.to_numeric, errors="raise")
        samples_used = Path(samples_path).resolve()
        metadata = pd.read_csv(samples_used, sep="\t", dtype=str)
        if not {"sample", "tissue_code"}.issubset(metadata):
            raise ValueError("Sample metadata requires sample and tissue_code columns.")
        if metadata[["sample", "tissue_code"]].isna().any().any() or metadata["sample"].duplicated().any():
            raise ValueError("Sample identifiers must be unique and tissue labels nonmissing.")
        metadata = metadata.set_index("sample")
        missing = matrix.index.difference(metadata.index)
        if len(missing):
            raise ValueError(f"Provide --samples metadata for all RNA samples; missing: {list(missing[:5])}")
        obs = metadata.loc[matrix.index].copy()
        obs["celltype"] = obs["tissue_code"]
        adata = ad.AnnData(matrix.to_numpy(dtype=np.float64), obs=obs,
                           var=pd.DataFrame(index=matrix.columns.astype(str)))
    if adata.obs_names.has_duplicates or adata.var_names.has_duplicates:
        raise ValueError("RNA sample and gene identifiers must be unique.")
    if adata.obs_names.isna().any() or adata.var_names.isna().any():
        raise ValueError("RNA sample and gene identifiers must be nonmissing.")
    if adata.n_obs * adata.n_vars > 50_000_000:
        raise ValueError("Subset RNA to at most 50 million sample-by-gene entries for this classroom workflow.")
    if sparse.issparse(adata.X):
        adata.X = adata.X.toarray()
    adata.X = np.asarray(adata.X, dtype=np.float64)
    if not np.isfinite(adata.X).all() or (adata.X < 0).any():
        raise ValueError("Use finite, nonnegative log-normalized RNA; do not supply centered data or missing values.")
    if adata.obs["celltype"].isna().any():
        raise ValueError("Sample groups must be nonmissing.")
    adata.obs["celltype"] = adata.obs["celltype"].astype(str)
    group_sizes = adata.obs["celltype"].value_counts()
    if group_sizes.min() < 2:
        raise ValueError("Each sample group needs at least two samples for the native stratified split.")
    try:
        train, test = train_test_split(
            list(adata.obs_names), test_size=NATIVE_TEST_FRACTION,
            random_state=NATIVE_SPLIT_SEED, stratify=adata.obs["celltype"],
        )
    except ValueError as error:
        raise ValueError(f"Sample groups cannot support GRETA's native stratified split: {error}") from error
    if len(test) < 10:
        raise ValueError("Native Omics evaluation needs at least 10 held-out observations.")
    training = adata[train, :].X
    fraction = (training >= config.min_expression).mean(axis=0)
    variance = training.var(axis=0, ddof=1)
    eligible = (fraction >= config.min_sample_fraction) & (variance > 1e-8)
    feature_info = pd.DataFrame(dict(
        gene=adata.var_names, train_fraction_expressed=fraction,
        train_variance=variance, passes_expression_filter=eligible,
    ))
    ordered = feature_info.loc[eligible].sort_values(
        ["train_variance", "gene"], ascending=[False, True], kind="stable",
    )
    if config.max_genes:
        ordered = ordered.head(config.max_genes)
    genes = sorted(ordered.gene.tolist())
    if len(genes) < 5:
        raise ValueError("Too few training-expressed genes; relax filtering or increase the gene limit.")
    feature_info["retained"] = feature_info.gene.isin(genes)
    adata = adata[:, genes].copy()
    adata.uns["rna_benchmark"] = dict(
        input=_input_record(rna_path),
        samples=_input_record(samples_used) if samples_used else None,
        split_seed=NATIVE_SPLIT_SEED, test_fraction=NATIVE_TEST_FRACTION,
        group_column="celltype", filtering="Training expression and variance only",
    )
    split = pd.DataFrame(dict(
        sample=adata.obs_names,
        celltype=adata.obs["celltype"].to_numpy(),
        partition=["train" if name in train else "test" for name in adata.obs_names],
    ))
    return adata, split, feature_info


def _unit_columns(matrix):
    centered = np.asarray(matrix, dtype=np.float64) - np.mean(matrix, axis=0)
    norm = np.linalg.norm(centered, axis=0)
    if (norm == 0).any():
        raise ValueError("A retained RNA feature is constant in training data.")
    return centered / norm


def build_rna_networks(adata, split, config=RnaConfig()):
    """Keep the complete eligible prior, matching its TF degrees in each baseline.

    Only TF-by-gene correlations are calculated; no gene-by-gene matrix is built.
    TF identities and outgoing edge budgets come from the prior, while inferred
    target identities and correlation signs come only from training expression.
    """
    config.validate()
    prior = gt.ds.read_db("hg38", "CollecTRI", verbose=False).rename(columns={"weight": "score"})
    prior = prior[["source", "target", "score"]].copy()
    prior = prior.loc[prior.source.isin(adata.var_names) & prior.target.isin(adata.var_names)
                      & (prior.source != prior.target)]
    prior = prior.groupby(["source", "target"], as_index=False).score.mean()
    counts = prior.groupby("source").target.nunique()
    tfs = counts[counts >= 3].index.tolist()
    train_names = split.loc[split.partition == "train", "sample"].tolist()
    training = adata[train_names, :].X
    if config.max_tfs:
        tf_variance = pd.Series(training.var(axis=0, ddof=1), index=adata.var_names)
        tfs = (tf_variance.loc[tfs].sort_index()
               .sort_values(ascending=False, kind="stable").head(config.max_tfs).index.tolist())
    tfs = sorted(tfs)
    if not tfs:
        raise ValueError("No eligible CollecTRI TF has at least three measured nonself targets.")
    prior = prior.loc[prior.source.isin(tfs)].sort_values(["source", "target"]).reset_index(drop=True)
    budgets = prior.groupby("source").size().to_dict()
    genes = adata.var_names.to_numpy()
    tf_positions = adata.var_names.get_indexer(tfs)
    networks = {}
    methods = [("Pearson", False)] + ([("Spearman", True)] if config.include_spearman else [])
    for name, use_ranks in methods:
        values = rankdata(training, axis=0) if use_ranks else training
        normalized = _unit_columns(values)
        correlations = np.clip(normalized[:, tf_positions].T @ normalized, -1.0, 1.0)
        rows = []
        for i, tf in enumerate(tfs):
            scores = correlations[i]
            eligible = genes != tf
            candidates = np.flatnonzero(eligible)
            # genes are alphabetic, so stable sorting gives deterministic ties.
            chosen = candidates[np.argsort(-np.abs(scores[candidates]), kind="stable")[:budgets[tf]]]
            rows.extend((tf, str(genes[j]), float(scores[j])) for j in chosen)
        networks[name] = pd.DataFrame(rows, columns=["source", "target", "score"])
    networks["CollecTRI"] = prior.copy()
    random_rng = np.random.default_rng(np.random.SeedSequence(config.seed))
    random_rows = []
    for tf in tfs:
        possible = genes[genes != tf]
        for target in random_rng.choice(possible, size=budgets[tf], replace=False):
            # Random signs make a signed null; they are not regulatory evidence.
            random_rows.append((tf, str(target), float(random_rng.choice([-1, 1]))))
    networks["Random"] = pd.DataFrame(random_rows, columns=["source", "target", "score"])
    for name, network in networks.items():
        if network.duplicated(["source", "target"]).any() or (network.source == network.target).any():
            raise RuntimeError(f"Duplicate/self edges in {name}.")
        if network.groupby("source").size().to_dict() != budgets:
            raise RuntimeError(f"Prior-matched TF degrees were not retained in {name}.")
    return networks


def validate_rna_grn(grn, adata):
    required = ["source", "target", "score"]
    if not set(required).issubset(grn):
        raise ValueError(f"A supplied RNA GRN needs columns {required}.")
    frame = grn[required].copy()
    if frame.isna().any().any():
        raise ValueError("Supplied GRN identifiers and scores must be nonmissing.")
    frame["score"] = pd.to_numeric(frame.score, errors="raise")
    if not np.isfinite(frame.score).all():
        raise ValueError("Supplied GRN scores must be finite.")
    if (set(frame.source) | set(frame.target)) - set(adata.var_names):
        raise ValueError("Every supplied GRN gene must occur in the retained RNA features.")
    frame = frame.loc[frame.source != frame.target]
    frame = frame.groupby(["source", "target"], as_index=False).score.mean()
    if len(frame) < 5:
        raise ValueError("A supplied GRN needs at least five nonself edges.")
    return frame


def _manifest(adata, config, terms, metrics, grn_paths):
    references = {}
    for name in set(metrics) | {"CollecTRI"}:
        filename = gt_config.DATA["hg38"]["dbs"][name]["fname"]
        if filename:
            path = Path(gt_config.PATH_DATA) / filename
            references[name] = dict(file=filename, sha256=_file_hash(path) if path.exists() else None)
    return dict(
        schema=1, config=asdict(config), terms=terms, metrics=list(metrics),
        inputs=adata.uns["rna_benchmark"],
        custom_grns={name: _input_record(path) for name, path in sorted(grn_paths.items())},
        references=references, packages={name: version(name) for name in PACKAGES},
        python=platform.python_version(),
        code={path.name: _file_hash(path) for path in (HERE / "rna_benchmark.py", HERE / "greta_benchmark.py")},
    )


def _save_outputs(output, results, coverage, networks, split, feature_info):
    output.mkdir(parents=True, exist_ok=True)
    (output / "networks").mkdir(exist_ok=True)
    # Remove only obsolete files that this runner previously recorded as its own.
    # In particular, a quick rerun must not retain an old mechanistic panel.
    previous_path = output / "run_info.json"
    try:
        previous = json.loads(previous_path.read_text()).get("output_hashes", {})
    except (OSError, ValueError, TypeError):
        previous = {}
    for relative in previous:
        path = output / relative
        if (path.parent == output / "networks" and path.suffix == ".tsv"
                and path.stem not in networks):
            path.unlink(missing_ok=True)
    for name, network in networks.items():
        network.to_csv(output / "networks" / f"{name}.tsv", sep="\t", index=False)
    ranked, omitted = comparable_results(results)
    tables = dict(
        evaluation_all=results, evaluation_ranked=ranked, coverage=coverage,
        excluded_from_ranking=omitted, inference_split=split, feature_info=feature_info,
        network_stats=gt.tl.stats(networks), overlap=gt.tl.ocoeff(networks),
        ranking=gt.tl.ranking(ranked).reset_index(),
    )
    for name, table in tables.items():
        table.to_csv(output / f"{name}.tsv", sep="\t", index=False)
    panels = [("ranking_task", ranked, "task", (18, 6)),
              ("ranking_class", ranked, "class", (8, 5))]
    for category in ("Predictive", "Literature", "Mechanistic"):
        panel = ranked.loc[ranked["class"] == category]
        if len(panel):
            panels.append((f"category_{category.lower()}", panel, "task",
                           (max(6, panel.db.nunique() * 1.5), 5)))
        else:
            for suffix in ("png", "pdf"):
                (output / f"category_{category.lower()}.{suffix}").unlink(missing_ok=True)
    for name, frame, level, figsize in panels:
        fig = gt.pl.ranking(frame, level=level, return_fig=True, figsize=figsize)
        for suffix in ("png", "pdf"):
            fig.savefig(output / f"{name}.{suffix}", dpi=160, bbox_inches="tight")
        plt.close(fig)
    artifacts = [output / f"{name}.tsv" for name in tables]
    artifacts += [output / "networks" / f"{name}.tsv" for name in networks]
    artifacts += [output / f"{name}.{suffix}" for name, *_ in panels for suffix in ("png", "pdf")]
    return ranked, omitted, artifacts


def _cache_is_valid(output, manifest):
    path = output / "run_info.json"
    try:
        info = json.loads(path.read_text())
        if info.get("manifest") != manifest or not info.get("output_hashes"):
            return False
        return all((output / name).is_file() and _file_hash(output / name) == digest
                   for name, digest in info["output_hashes"].items())
    except (OSError, ValueError, TypeError):
        return False


def run_rna(config=RnaConfig(), output_dir=None, rna_path=RNA_PATH,
            samples_path=SAMPLES_PATH, terms=None, metrics=None, grn_paths=None,
            force=False, verbose=True):
    """Infer full RNA networks and run/reuse a fingerprinted native GRETA benchmark.

    Cached scores are reused only when input, config, code, package and native
    reference hashes match and all saved result/figure hashes remain intact.
    """
    started = time.perf_counter()
    configure_cache()
    terms = {} if terms is None else terms
    if not isinstance(terms, dict):
        raise ValueError("terms must map reference database names to context-term lists.")
    grn_paths = grn_paths or {}
    output = Path(output_dir or HERE / "results/rna").resolve()
    if verbose:
        print("Loading RNA and selecting genes using training samples...", flush=True)
    adata, split, feature_info = load_rna(config, rna_path, samples_path)
    loaded_at = time.perf_counter()
    networks = build_rna_networks(adata, split, config)
    for name, path in grn_paths.items():
        if name in networks or Path(name).name != name or name in {"", ".", ".."}:
            raise ValueError("Use distinct, simple names for supplied GRNs.")
        networks[name] = validate_rna_grn(pd.read_csv(path, sep=None, engine="python"), adata)
    inferred_at = time.perf_counter()
    feature_info["is_tf"] = feature_info.gene.isin(networks["CollecTRI"].source)
    selected = list(metrics or RNA_METRICS)
    if config.profile == "quick":
        selected = [metric for metric in selected if metric != "KnockTF (scoring)"]
    selected = list(dict.fromkeys(selected + list(RNA_DIAGNOSTICS)))
    forbidden = set(selected) - (set(RNA_METRICS) | set(RNA_DIAGNOSTICS))
    if forbidden:
        raise ValueError(f"Unsupported RNA tutorial metrics: {sorted(forbidden)}")
    manifest = _manifest(adata, config, terms, selected, grn_paths)
    if not force and _cache_is_valid(output, manifest):
        if verbose:
            print("Reusing scores and figures from the verified matching cache.", flush=True)
        results = pd.read_csv(output / "evaluation_all.tsv", sep="\t")
        ranked = pd.read_csv(output / "evaluation_ranked.tsv", sep="\t")
        coverage = pd.read_csv(output / "coverage.tsv", sep="\t")
        info = json.loads((output / "run_info.json").read_text())
        elapsed = time.perf_counter() - started
        return dict(adata=adata, networks=networks, split=split, feature_info=feature_info,
                    results=results, res=ranked, coverage=coverage, output_dir=output,
                    timing=dict(cache_load_seconds=elapsed), elapsed_seconds=elapsed,
                    run_info=info, cache_hit=True)
    if verbose:
        print(f"Retained {adata.n_vars:,} genes and {networks['CollecTRI'].source.nunique():,} TFs; "
              f"each baseline has {len(networks['CollecTRI']):,} edges.", flush=True)
    result_frames, coverage_frames, method_times = [], [], {}
    for name, network in networks.items():
        if verbose:
            print(f"Computing native GRETA tasks for {name} ({network.target.nunique():,} targets)...", flush=True)
        method_started = time.perf_counter()
        frame, captured = evaluate_networks(
            adata, {name: network}, terms, selected, dataset_name="NCI60 RNA",
        )
        result_frames.append(frame)
        coverage_frames.append(captured)
        method_times[name] = time.perf_counter() - method_started
        if verbose:
            print(f"  {name} finished in {method_times[name]:.1f} seconds.", flush=True)
    results = pd.concat(result_frames, ignore_index=True)
    coverage = pd.concat(coverage_frames, ignore_index=True)
    evaluated_at = time.perf_counter()
    ranked, omitted, artifacts = _save_outputs(output, results, coverage, networks, split, feature_info)
    timing = dict(
        load_and_filter_seconds=loaded_at - started,
        inference_seconds=inferred_at - loaded_at,
        evaluation_seconds=evaluated_at - inferred_at,
        export_seconds=time.perf_counter() - evaluated_at,
        total_seconds=time.perf_counter() - started,
        evaluation_by_method_seconds=method_times,
    )
    # Populate hashes after native evaluation has downloaded any missing resources.
    initial_code_hashes = manifest["code"]
    manifest = _manifest(adata, config, terms, selected, grn_paths)
    # An edit while a long benchmark is running must invalidate the next cache
    # lookup rather than label old in-memory functions with new source hashes.
    manifest["code"] = initial_code_hashes
    info = dict(
        manifest=manifest, config=asdict(config), packages=manifest["packages"], timing=timing,
        dataset=dict(name="NCI60 RNA", n_samples=adata.n_obs, n_genes=adata.n_vars,
                     n_train=int((split.partition == "train").sum()),
                     n_test=int((split.partition == "test").sum()),
                     groups=adata.obs.celltype.value_counts().to_dict()),
        construction="Full eligible CollecTRI; inferred/random TF degrees matched to that prior; all retained genes are candidate targets",
        inference="Pearson correlation using only native training samples; optional Spearman; correlation signs are association proxies",
        random="Uniform nonself target sampling per TF with random signs, preserving prior outgoing degrees",
        evaluation="Native GRETA RNA-only tasks; no genomic class or ATAC predictors",
        context="Broad mixed-tumor reference screen; no PBMC terms" if not terms else "User-provided reference term filters",
        split="Native Omics split: test_size=0.33, random_state=42, stratified by real tissue/sample groups",
        zero_handling="Native Omics excludes zero-valued target responses; no missing values are imputed",
        custom_grn_note="Supplied networks may already have used evaluation samples; verify their inference provenance" if grn_paths else None,
        ranked_databases=sorted(ranked.db.unique()),
        excluded_databases=sorted(omitted.db.unique()) if len(omitted) else [],
        output_hashes={str(path.relative_to(output)): _file_hash(path) for path in artifacts},
    )
    (output / "run_info.json").write_text(json.dumps(info, indent=2) + "\n")
    return dict(adata=adata, networks=networks, split=split, feature_info=feature_info,
                results=results, res=ranked, coverage=coverage, output_dir=output,
                timing=timing, elapsed_seconds=timing["total_seconds"], run_info=info,
                cache_hit=False)
