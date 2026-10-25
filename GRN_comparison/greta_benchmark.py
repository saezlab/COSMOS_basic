"""Small, genuine GRETA evaluation with task-level ranking.

Metric algorithms and plots come from gretapy 0.0.3. The compatibility context
repairs dispatcher labels, excludes undefined test correlations, and records coverage.
"""

from contextlib import contextmanager
from dataclasses import asdict, dataclass
from hashlib import sha256
from importlib.metadata import version
import json
import inspect
from pathlib import Path
import shutil
import time

import gretapy as gt
import gretapy._utils as gt_utils
import gretapy.config as gt_config
import gretapy.ds._db as gt_db
import gretapy.ds._dts as gt_dts
import gretapy.tl._eval as gt_eval
import gretapy.tl._predictive as gt_predictive
from gretapy.tl._utils import _f_beta_score
import matplotlib.pyplot as plt
import mudata as mu
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.model_selection import train_test_split


HERE = Path(__file__).resolve().parent
QUICK_METRICS = (
    "ENCODE CREs", "ENCODE Blacklist", "HPA", "TF-Marker", "IntAct", "Europe PMC",
    "Hallmarks", "Gene ~ TFs", "Gene ~ CREs", "CRE ~ TFs",
)
TEACHING_METRICS = QUICK_METRICS + ("KnockTF (scoring)",)
DIAGNOSTICS = ("Promoters", "CollecTRI")
OMICS_KEYS = {("source", "target"): "Gene ~ TFs", ("cre", "target"): "Gene ~ CREs",
              ("source", "cre"): "CRE ~ TFs"}


@dataclass(frozen=True)
class GretaConfig:
    n_cells: int = 180
    seed: int = 42
    edges_per_tf: int = 5
    profile: str = "teaching"


def configure_cache(cache_dir=None):
    """Keep upstream downloads in this practical, independent of notebook cwd."""
    path = Path(cache_dir or HERE / "gretapy_data").resolve()
    path.mkdir(parents=True, exist_ok=True)
    for module in (gt_config, gt_utils, gt_db, gt_dts):
        module.PATH_DATA = str(path)
    # This small, attributed snapshot makes the published comparison available offline.
    snapshot = HERE / "data/greta_published_metrics.csv.gz"
    if snapshot.exists() and not (path / "metrics.csv.gz").exists():
        shutil.copyfile(snapshot, path / "metrics.csv.gz")
    return path


@contextmanager
def evaluation_compatibility(coverage):
    """Repair 0.0.3 routing without changing the native metric calculations.

    Restored on exit; intended for sequential notebook/CLI evaluation only.
    """
    if version("gretapy") != "0.0.3":
        raise RuntimeError("Use the tested requirements-lock.txt (gretapy 0.0.3).")
    original_metric = gt_eval._run_metric
    original_fileless = gt_eval._run_fileless_metric
    original_prediction = gt_predictive._test_predictability

    def metric(metric_type, *args, **kwargs):
        aliases = {"Gene Sets": "Gene sets", "CRE to Gene": "CRE to gene links"}
        if metric_type == "TF Scoring":
            # Native ULM needs >=3 targets measured in KnockTF, not merely >=3
            # original edges. Align features before its existing eligibility check.
            bound = inspect.signature(original_metric).bind(metric_type, *args, **kwargs)
            graph, reference = bound.arguments["grn"], bound.arguments["db"]
            bound.arguments["grn"] = graph[graph.target.isin(reference.var_names)]
            eligible, supported, _ = _knocktf_support(
                bound.arguments["grn"], bound.arguments["adata"], reference, bound.arguments["cats"],
            )
            if not supported:
                return np.nan, np.nan, np.nan
            # ULM drops zero-response genes in each experiment. Evaluate only
            # experiments with >=3 usable targets, keeping the original recall universe.
            bound.arguments["db"] = reference[supported, :].copy()
            precision, _, _ = original_metric(*bound.args, **bound.kwargs)
            recall = precision * len(supported) / len(eligible)
            return precision, recall, _f_beta_score(prc=precision, rcl=recall)
        return original_metric(aliases.get(metric_type, metric_type), *args, **kwargs)

    def fileless(metric_type, *args, **kwargs):
        if metric_type == "Steady State Simulation":
            metric_type = "Steady state simulation"
        return original_fileless(metric_type, *args, **kwargs)

    def prediction(*args, **kwargs):
        # Native 0.0.3 passes constant-response Spearman NaNs to BH and crashes.
        # Such features have no defined correlation; omit them before native fitting.
        dataset = kwargs["data"]
        response = dataset.mod[kwargs["mod_target"]] if isinstance(dataset, mu.MuData) else dataset
        constant = []
        for target in kwargs["grn"][kwargs["col_target"]].unique():
            y = response[kwargs["test_obs_names"], target].X.ravel()
            y = y[y != 0]
            if y.size >= 10 and np.ptp(y) == 0:
                constant.append(target)
        safe_kwargs = dict(kwargs)
        safe_kwargs["grn"] = kwargs["grn"][~kwargs["grn"][kwargs["col_target"]].isin(constant)]
        frame = original_prediction(*args, **safe_kwargs)
        key = (kwargs["col_source"], kwargs["col_target"])
        coverage[OMICS_KEYS[key]] = dict(
            n_tested=len(frame),
            n_hits=int(((frame.padj < 0.05) & (frame.coef > 0.05)).sum()),
            n_skipped_constant_response=len(constant),
        )
        return frame

    gt_eval._run_metric, gt_eval._run_fileless_metric = metric, fileless
    gt_predictive._test_predictability = prediction
    try:
        yield
    finally:
        gt_eval._run_metric, gt_eval._run_fileless_metric = original_metric, original_fileless
        gt_predictive._test_predictability = original_prediction


def load_dataset(config=GretaConfig(), dataset_path=None):
    if config.profile not in {"quick", "teaching"} or config.edges_per_tf < 3:
        raise ValueError("Use profile quick/teaching and edges_per_tf >=3.")
    if dataset_path is None:
        if config.n_cells < 90:
            raise ValueError("Use at least 90 toy cells so Omics has enough held-out observations.")
        mdata, _ = gt.ds.toy(n_cells=config.n_cells, seed=config.seed)
        # Discard the hardcoded example GRN. Infer fresh baseline edges below.
    else:
        mdata = mu.read_h5mu(dataset_path)
    if not {"rna", "atac"}.issubset(mdata.mod) or "celltype" not in mdata.obs:
        raise ValueError("MuData requires rna/atac modalities and mdata.obs['celltype'].")
    if mdata.obs_names.has_duplicates or mdata.obs.celltype.isna().any():
        raise ValueError("Cell identifiers must be unique, with nonmissing celltype labels.")
    for name in ("rna", "atac"):
        modality = mdata.mod[name]
        if not modality.obs_names.equals(mdata.obs_names) or modality.var_names.has_duplicates:
            raise ValueError("Modalities must have the same ordered cells and unique feature names.")
        # The pinned inference/predictive helpers expect dense matrices. Bound the conversion.
        if modality.n_obs * modality.n_vars > 50_000_000:
            raise ValueError("For this classroom run, first subset cells/features to a manageable MuData.")
        if sparse.issparse(modality.X):
            modality.X = modality.X.toarray()
        if not np.isfinite(modality.X).all() or (modality.X < 0).any():
            raise ValueError("Use finite, nonnegative normalized RNA/ATAC values, not centered/scaled data.")
    return mdata


def validate_grn(grn, mdata):
    """Validate a student's TF–CRE–target table against the matched measurements."""
    required = ["source", "cre", "target", "score"]
    if not set(required).issubset(grn):
        raise ValueError(f"GRN requires columns {required}.")
    grn = grn[required].copy()
    if grn.isna().any().any():
        raise ValueError("GRN contains missing identifiers or scores.")
    grn["score"] = pd.to_numeric(grn.score, errors="raise")
    if not np.isfinite(grn.score).all():
        raise ValueError("GRN scores must be finite numbers.")
    if not grn.cre.str.fullmatch(r"chr[0-9A-Za-z]+-\d+-\d+").all():
        raise ValueError("CRE names must use chrX-start-end in the correct genome assembly.")
    intervals = grn.cre.str.rsplit("-", n=2, expand=True)
    if not (intervals[2].astype(int) > intervals[1].astype(int)).all():
        raise ValueError("Each CRE end coordinate must exceed its start coordinate.")
    genes = set(mdata.mod["rna"].var_names)
    if (set(grn.source) | set(grn.target)) - genes:
        raise ValueError("Every source/target must occur in the RNA feature names.")
    if set(grn.cre) - set(mdata.mod["atac"].var_names):
        raise ValueError("Every GRN CRE must occur in the ATAC feature names.")
    grn = grn[grn.source != grn.target]
    grn = grn.groupby(["source", "cre", "target"], as_index=False).score.mean()
    if len(grn) < 5:
        raise ValueError("At least five non-self edges are required for evaluation.")
    return grn


def build_networks(mdata, config=GretaConfig()):
    """Public GRETA inference helpers; matched source degrees for readable comparisons."""
    train, test = train_test_split(
        list(mdata.obs_names), test_size=0.33, random_state=42,
        stratify=mdata.obs["celltype"],
    )
    # Match the pinned Omics evaluator's split, so test cells do not select edges.
    training = mdata[train, :].copy()
    prior = gt.mt.lit_grn(mdata, min_targets=0)
    prior = prior[prior.source != prior.target]
    counts = prior.groupby("source").target.nunique()
    tfs = sorted(counts[counts >= 3].index)
    if not tfs:
        raise ValueError("No prior regulators have >=3 eligible targets in this dataset.")
    candidates = {}
    for name, method in [("Pearson", "pearson"), ("Spearman", "spearman")]:
        candidates[name] = gt.mt.correlation(
            training, tfs=tfs, method=method, thr_r=0, min_targets=0,
        )
    # Use a common single accessible promoter per target (an explicit baseline assumption).
    mapping = candidates["Pearson"][["cre", "target"]].sort_values(["target", "cre"]).drop_duplicates("target")
    prior = prior.merge(mapping, on=["target", "cre"])
    counts = prior.groupby("source").target.nunique()
    tfs = sorted(counts[counts >= 3].index)
    budgets = {tf: min(config.edges_per_tf, int(counts[tf])) for tf in tfs}
    if not tfs:
        raise ValueError("No common promoter-linked regulators remain after feature matching.")
    networks = {}
    for name, frame in candidates.items():
        frame = frame.groupby(["source", "target"], as_index=False).score.mean()
        frame["abs_score"] = frame.score.abs()
        selected = []
        for tf in tfs:
            group = frame[frame.source == tf].sort_values(["abs_score", "target"], ascending=[False, True])
            selected.append(group.head(budgets[tf]).drop(columns="abs_score"))
        networks[name] = pd.concat(selected).merge(mapping, on="target")
    prior_rng, random_rng = [np.random.default_rng(s) for s in np.random.SeedSequence(config.seed).spawn(2)]
    selected, random_rows = [], []
    for tf in tfs:
        group = prior[prior.source == tf].sort_values("target").drop_duplicates("target")
        selected.append(group.iloc[prior_rng.choice(len(group), budgets[tf], replace=False)])
        possible = sorted(set(mapping.target) - {tf})
        random_rows.extend((tf, target, 1.0) for target in random_rng.choice(possible, budgets[tf], replace=False))
    networks["CollecTRI"] = pd.concat(selected, ignore_index=True)
    networks["Random"] = pd.DataFrame(random_rows, columns=["source", "target", "score"]).merge(mapping, on="target")
    for name in networks:
        networks[name] = validate_grn(networks[name], mdata)
        actual = networks[name].groupby("source").target.nunique().to_dict()
        if actual != budgets:
            raise RuntimeError(f"Matched TF–target budgets were not met for {name}.")
    split = pd.DataFrame({"cell": list(mdata.obs_names),
                          "partition": ["train" if c in train else "test" for c in mdata.obs_names]})
    return networks, split, mapping


def _knocktf_support(grn, adata, db, cats):
    mask = db.obs.TF.isin(adata.var_names) & (db.obs.logFC < -0.5)
    if cats is not None:
        mask &= db.obs["Tissue.Type"].isin(cats)
    eligible = db.obs_names[mask].tolist()
    targets = grn[grn.target.isin(db.var_names)].groupby("source").target.unique().to_dict()
    supported, insufficient = [], 0
    for experiment in eligible:
        tf = db.obs.loc[experiment, "TF"]
        genes = targets.get(tf, [])
        if len(genes) < 3:
            insufficient += int(tf in targets)
            continue
        values = db[[experiment], list(genes)].X
        values = values.toarray() if sparse.issparse(values) else np.asarray(values)
        if np.count_nonzero(values) >= 3:
            supported.append(experiment)
        else:
            insufficient += 1
    return eligible, supported, insufficient


def tf_scoring_coverage(grn, mdata, terms):
    adata = mdata.mod["rna"] if isinstance(mdata, mu.MuData) else mdata
    db = gt.ds.read_db("hg38", "KnockTF (scoring)")
    eligible, supported, insufficient = _knocktf_support(grn, adata, db, terms.get("KnockTF (scoring)"))
    return dict(n_eligible=len(eligible), n_tested=len(supported),
                n_skipped_insufficient_targets=insufficient)


def evaluate_networks(mdata, networks, terms, metrics=TEACHING_METRICS, dataset_name="Toy PBMC"):
    """Run actual gt.tl.benchmark scores and expose untested results explicitly."""
    metrics = list(dict.fromkeys(metrics))
    frames, records = [], []
    for name, network in networks.items():
        captured = {}
        supported_tfa = tf_scoring_coverage(network, mdata, terms) if "KnockTF (scoring)" in metrics else None
        active_metrics = [m for m in metrics if not (m == "KnockTF (scoring)" and supported_tfa["n_tested"] == 0)]
        with evaluation_compatibility(captured):
            result = gt.tl.benchmark(
                grns={name: {"hg38": {dataset_name: network}}},
                datasets={dataset_name: mdata}, terms={"hg38": {dataset_name: terms}},
                metrics=active_metrics, verbose=False,
            ).rename(columns={"grn": "name"})
        if set(result.db) != set(active_metrics):
            raise RuntimeError(f"GRETA omitted requested tasks for {name}: {set(active_metrics) - set(result.db)}")
        if supported_tfa is not None and supported_tfa["n_tested"] == 0:
            result = pd.concat([result, pd.DataFrame([dict(
                name=name, organism="hg38", dataset=dataset_name, **{"class": "Mechanistic"},
                task="TF Scoring", db="KnockTF (scoring)", precision=np.nan, recall=np.nan, f01=np.nan,
            )])], ignore_index=True)
        result["native_f01"] = result.f01
        for row in result.itertuples():
            info = captured.get(row.db, supported_tfa if row.db == "KnockTF (scoring)" else {})
            info = info or {}
            status = "unsupported" if info.get("n_tested") == 0 or pd.isna(row.f01) else "evaluated"
            if status == "unsupported":
                result.loc[result.db == row.db, ["precision", "recall", "f01"]] = np.nan
            elif row.db == "KnockTF (scoring)":
                info = dict(info, n_hits=int(round(row.precision * info["n_tested"])))
            records.append(dict(name=name, db=row.db, status=status, **info))
        frames.append(result)
    return pd.concat(frames, ignore_index=True), pd.DataFrame(records)


def comparable_results(results):
    """Exclude circular checks and columns without comparable coverage for all methods."""
    names = set(results.name)
    keep, omitted = [], []
    keys = ["organism", "dataset", "class", "task", "db"]
    for _, group in results.groupby(keys, sort=False):
        reason = ""
        if group.db.iloc[0] in DIAGNOSTICS:
            reason = "Construction/self-reference diagnostic; excluded from ranking"
        elif set(group.name) != names or group.f01.isna().any():
            reason = "Incomplete evaluation coverage; excluded for all methods"
        if reason:
            omitted.append(group.assign(reason=reason))
        else:
            keep.append(group)
    if not keep:
        raise ValueError("No comparable evaluation columns remain for ranking.")
    return pd.concat(keep, ignore_index=True), pd.concat(omitted, ignore_index=True) if omitted else pd.DataFrame()


def published_results():
    configure_cache()
    return gt.ds.read_metrics()


def save_evaluation(results, coverage, networks, split, mapping, config, output_dir=None, metadata=None, terms=None):
    output = Path(output_dir or HERE / "results/greta")
    output.mkdir(parents=True, exist_ok=True)
    (output / "networks").mkdir(exist_ok=True)
    for name, network in networks.items():
        network.to_csv(output / "networks" / f"{name}.tsv", sep="\t", index=False)
    ranked, omitted = comparable_results(results)
    tables = dict(evaluation_all=results, evaluation_ranked=ranked, coverage=coverage,
                  excluded_from_ranking=omitted, inference_split=split, candidate_promoters=mapping,
                  network_stats=gt.tl.stats(networks), ranking=gt.tl.ranking(ranked).reset_index())
    for name, frame in tables.items():
        frame.to_csv(output / f"{name}.tsv", sep="\t", index=False)
    for level in ("task", "class"):
        fig = gt.pl.ranking(ranked, level=level, return_fig=True,
                            figsize=(16, 6) if level == "task" else (7, 5))
        for suffix in ("png", "pdf"):
            fig.savefig(output / f"ranking_{level}.{suffix}", dpi=160, bbox_inches="tight")
        plt.close(fig)
    info = dict(config=asdict(config), packages={p: version(p) for p in ("gretapy", "decoupler", "mudata", "numpy", "pandas")},
                compatibility="0.0.3 label routing repaired; undefined constant-response tests excluded; native algorithms/formulas retained",
                dataset="Synthetic toy PBMC" if metadata is None else metadata,
                reference_filter_terms=terms,
                inference_holdout="Omics split: test_size=0.33, random_state=42, stratified by celltype",
                scores="Native GRETA precision/recall/F0.1; unsupported values are NA",
                ranked_databases=sorted(ranked.db.unique()), excluded_databases=sorted(omitted.db.unique()) if len(omitted) else [])
    (output / "run_info.json").write_text(json.dumps(info, indent=2) + "\n")
    return output


def run_greta(config=GretaConfig(), output_dir=None, dataset_path=None, grn_paths=None, terms=None, metrics=None):
    started = time.perf_counter()
    configure_cache()
    mdata = load_dataset(config, dataset_path)
    if dataset_path is not None and terms is None:
        raise ValueError("For custom data provide context-specific terms (JSON or a named GRETA context).")
    terms = gt.get_terms("PBMC") if terms is None else terms
    networks, split, mapping = build_networks(mdata, config)
    for name, path in (grn_paths or {}).items():
        if name in networks:
            raise ValueError(f"Choose a distinct method name for {name}.")
        networks[name] = validate_grn(pd.read_csv(path, sep=None, engine="python"), mdata)
    selected = list(metrics or (TEACHING_METRICS if config.profile == "teaching" else QUICK_METRICS))
    selected += list(DIAGNOSTICS)
    results, coverage = evaluate_networks(mdata, networks, terms, selected,
                                          "Toy PBMC" if dataset_path is None else Path(dataset_path).stem)
    metadata = None if dataset_path is None else dict(path=str(Path(dataset_path).resolve()),
               sha256=sha256(Path(dataset_path).read_bytes()).hexdigest(),
               custom_grn_note="Provided GRNs may have used evaluation cells during inference; verify provenance.")
    output = save_evaluation(results, coverage, networks, split, mapping, config, output_dir, metadata, terms)
    return dict(mdata=mdata, networks=networks, results=results, coverage=coverage,
                output_dir=output, elapsed_seconds=time.perf_counter() - started)
