#!/usr/bin/env python3
"""Compare the supplied SCENIC+ GRN, CollecTRI and a random control with GRETA.

Run from any directory. The default uses cached GRETA resources; pass
--download-resources to fetch missing resources on the first run. Primary
rankings use full networks on a shared TF set. A separate equal-degree analysis
and Gene Sets background sensitivity expose important design dependencies.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from importlib.metadata import version
import json
import os
from pathlib import Path
import sys
import tempfile
import time

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "cosmos-heart-matplotlib"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(Path(tempfile.gettempdir()) / "cosmos-heart-numba"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import anndata as ad
import decoupler as dc
import gretapy as gt
import gretapy.config as gt_config
import numpy as np
import pandas as pd
from scipy import sparse, stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "GRN_comparison"))
from greta_benchmark import configure_cache, evaluate_networks, _knocktf_support
from diagnose_gene_sets import infer_pathways, summarize
from gretapy.tl._prior import _find_pairs

METRICS = ("Hallmarks", "PROGENy", "Reactome", "HPA", "TF-Marker", "IntAct", "Europe PMC", "KnockTF (scoring)")
GENE_SETS = ("Hallmarks", "PROGENy", "Reactome")
MARKERS = ("HPA", "TF-Marker")
METHODS = ("SCENIC+", "CollecTRI", "Random")


def digest(path):
    h = sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_data(path):
    data = ad.read_h5ad(path)
    if data.obs_names.has_duplicates or data.var_names.has_duplicates:
        raise ValueError("Pseudobulk profile and gene names must be unique.")
    if "celltype" not in data.obs:
        if "cell_type" not in data.obs:
            raise ValueError("Pseudobulk requires obs['celltype'] or obs['cell_type'].")
        data.obs["celltype"] = data.obs["cell_type"].copy()
    if data.obs.celltype.isna().any():
        raise ValueError("Cell type labels must be nonmissing.")
    data.X = data.X.toarray() if sparse.issparse(data.X) else np.asarray(data.X)
    data.X = np.asarray(data.X, dtype=float)
    if not np.isfinite(data.X).all() or (data.X < 0).any():
        raise ValueError("Use nonnegative log-normalized pseudobulk RNA, not TF-input Z-scores.")
    if data.n_obs < 2 or data.n_vars < 20:
        raise ValueError("At least two pseudobulk profiles and 20 measured genes are needed.")
    return data


def prepare_networks(data, scenic_path, prior_path, seed=42, min_targets=3):
    """Drop sign conflicts, use shared eligible sources and preserve full degrees.

    Opposing-sign duplicate SCENIC+ pairs are excluded in their entirety. An
    edge's score is its supplied strength times its activator/repressor sign.
    n_topics is retained in the conflict audit but is not a regulatory weight.
    """
    if min_targets < 3 or seed < 0:
        raise ValueError("Use min_targets >=3 and a nonnegative random seed.")
    raw = pd.read_csv(scenic_path, sep="\t")
    required = ["TF", "gene", "type", "score", "n_topics"]
    if not set(required).issubset(raw):
        raise ValueError(f"SCENIC+ TSV needs {required}.")
    if raw[required].isna().any().any():
        raise ValueError("SCENIC+ identifiers, signs, scores and n_topics must be nonmissing.")
    if not raw.type.isin(["activator", "repressor"]).all():
        raise ValueError("SCENIC+ type must be activator or repressor.")
    raw["score"] = pd.to_numeric(raw.score, errors="raise")
    if not np.isfinite(raw.score).all() or (raw.score < 0).any():
        raise ValueError("SCENIC+ strength scores must be finite and nonnegative.")
    signs_per_pair = raw.groupby(["TF", "gene"])["type"].transform("nunique")
    conflicts = raw.loc[signs_per_pair > 1].copy()
    clean = raw.loc[signs_per_pair == 1].copy()
    # Same-sign repeated evidence: keep the strongest supplied strength, with
    # n_topics as a deterministic tie breaker. There are no such duplicates in
    # the provided dataset; opposed signs never get averaged or voted on.
    clean = clean.sort_values(["TF", "gene", "score", "n_topics"], ascending=[True, True, False, False])
    clean = clean.drop_duplicates(["TF", "gene"])
    clean["score"] *= clean.type.map({"activator": 1.0, "repressor": -1.0})
    scenic = clean.rename(columns={"TF": "source", "gene": "target"})[["source", "target", "score"]]
    prior = pd.read_csv(prior_path, sep="\t").rename(columns={"weight": "score", "mor": "score"})
    if not {"source", "target", "score"}.issubset(prior):
        raise ValueError("CollecTRI needs source,target,weight (or score) columns.")
    prior = prior[["source", "target", "score"]].copy()
    if prior.isna().any().any() or not np.isfinite(prior.score.to_numpy(float)).all():
        raise ValueError("CollecTRI edges must be finite and nonmissing.")
    if prior.duplicated(["source", "target"]).any():
        raise ValueError("Resolve duplicate CollecTRI source-target pairs before benchmarking.")
    genes = set(data.var_names)
    frames, filtering = {}, []
    for name, frame in (("SCENIC+", scenic), ("CollecTRI", prior)):
        eligible = frame.source.isin(genes) & frame.target.isin(genes) & frame.source.ne(frame.target)
        retained = frame.loc[eligible].copy()
        filtering.append(dict(name=name, input_rows=len(frame), measured_nonself_edges=len(retained),
                              n_measured_sources=retained.source.nunique()))
        frames[name] = retained
    degrees = pd.DataFrame({name: frame.groupby("source").size() for name, frame in frames.items()}).fillna(0).astype(int)
    degrees["retained_tf"] = (degrees >= min_targets).all(axis=1)
    shared = sorted(degrees.index[degrees.retained_tf])
    if len(shared) < 2:
        raise ValueError("Need at least two shared TFs with enough measured targets in both networks.")
    networks = {name: frame.loc[frame.source.isin(shared)].sort_values(["source", "target"]).reset_index(drop=True)
                for name, frame in frames.items()}
    networks["Random"] = random_network(networks["SCENIC+"], data.var_names, seed)
    matched, rng = {}, np.random.default_rng(np.random.SeedSequence([seed, 2026]))
    budgets = degrees.loc[shared, ["SCENIC+", "CollecTRI"]].min(axis=1).to_dict()
    for name in ("SCENIC+", "CollecTRI"):
        chosen = []
        for tf, group in networks[name].groupby("source", sort=True):
            group = group.copy()
            group["abs_score"] = group.score.abs()
            # Shuffle equal-confidence ties reproducibly before a stable sort.
            # This matters for CollecTRI, whose strengths are all |score| = 1.
            if name == "CollecTRI":
                group = group.iloc[rng.permutation(len(group))]
                group = group.sort_values("abs_score", ascending=False, kind="stable")
            else:
                group = group.sort_values(["abs_score", "target"], ascending=[False, True], kind="stable")
            chosen.append(group.head(budgets[tf]).drop(columns="abs_score"))
        matched[name] = pd.concat(chosen, ignore_index=True)
    matched["Random"] = random_network(matched["SCENIC+"], data.var_names, seed)
    for name, frame in matched.items():
        if frame.groupby("source").size().to_dict() != budgets:
            raise RuntimeError(f"Equal outgoing degrees failed for {name}.")
    audit = dict(conflicts=conflicts, tf_degrees=degrees.reset_index(names="source"), filtering=pd.DataFrame(filtering),
                 raw_rows=len(raw), conflict_pairs=int(conflicts[["TF", "gene"]].drop_duplicates().shape[0]),
                 conflict_rows=len(conflicts), min_targets=min_targets, shared_tfs=len(shared), seed=seed)
    return networks, matched, audit


def random_network(template, gene_names, seed):
    """Preserve each TF's outgoing degree and signed score multiset exactly."""
    rng = np.random.default_rng(np.random.SeedSequence(seed))
    genes = np.asarray(sorted(gene_names), dtype=str)
    rows = []
    for tf, group in template.groupby("source", sort=True):
        possible = genes[genes != tf]
        targets = rng.choice(possible, size=len(group), replace=False)
        scores = rng.permutation(group.score.to_numpy(float))
        rows.extend(zip([tf] * len(group), targets, scores, strict=True))
    return pd.DataFrame(rows, columns=["source", "target", "score"])


def select_references(cache, quick=False, download=False):
    metrics, records = [], []
    for name in METRICS:
        filename = gt_config.DATA["hg38"]["dbs"][name]["fname"]
        path = cache / filename
        skip = quick and name == "KnockTF (scoring)"
        if skip:
            status, reason = "skipped", "--quick excludes the large external KnockTF resource"
        elif path.exists():
            status, reason = "available", "cached resource"
        elif download:
            print(f"Downloading GRETA reference: {name}", flush=True)
            gt.ds.read_db("hg38", name, verbose=True)
            status, reason = "available", "downloaded resource"
        else:
            status, reason = "skipped", "missing resource; rerun with --download-resources"
        if status == "available":
            metrics.append(name)
        records.append(dict(db=name, status=status, reason=reason, file=filename,
                            bytes=path.stat().st_size if path.exists() else 0,
                            sha256=digest(path) if path.exists() else None))
    if not metrics:
        raise ValueError("No GRETA resources are cached. Rerun with --download-resources (KnockTF is approximately 153 MB; --quick omits it).")
    return metrics, pd.DataFrame(records)


def pathway_audit(data, networks, dbname):
    """Reconstruct native pathway calls and repeat ORA with eligible-gene sets.

    Only the network-side ORA universe changes in the sensitivity analysis.
    The RNA-side ULM, global BH, thresholds, and >20% profile rule remain native.
    """
    reference = gt.ds.read_db("hg38", dbname, verbose=False)
    measured = data.copy()
    dc.mt.ulm(measured, reference, verbose=False)
    scores = measured.obsm["score_ulm"]
    raw = 2 * stats.t.sf(np.abs(scores.to_numpy()), data.n_vars - 2)
    q = stats.false_discovery_control(raw.ravel(), method="bh").reshape(raw.shape)
    fraction = pd.Series(((q < .01) & (scores.to_numpy() > 0)).mean(axis=0), index=scores.columns)
    expression_hits = set(fraction.index[fraction > .2])
    pathways = {source: set(targets) for source, targets in reference.groupby("source").target}
    genes = set(data.var_names)
    rows, calls = [], []
    for background, universe, exclude_self in (
        ("native_20000", pathways, False),
        ("eligible_nonself", {pw: targets & genes for pw, targets in pathways.items()}, True),
    ):
        for name, network in networks.items():
            predicted = infer_pathways(network, universe, 20000 if not exclude_self else len(genes), exclude_self=exclude_self)
            rows.append(dict(db=dbname, background=background, name=name, counted_unit="pathway labels",
                             n_rna_tested_pathways=len(scores.columns), n_reference_pathways=len(pathways),
                             **summarize(predicted, expression_hits)))
            for pathway in sorted(pathways):
                calls.append(dict(db=dbname, background=background, name=name, pathway=pathway,
                                  rna_tested=pathway in scores.columns,
                                  positive_rna_fraction=float(fraction.get(pathway, 0)),
                                  expression_hit=pathway in expression_hits, network_call=pathway in predicted))
    return pd.DataFrame(rows), pd.DataFrame(calls)


def marker_counts(networks, genes, terms, dbname):
    reference = gt.ds.read_db("hg38", dbname, verbose=False)
    reference = reference[reference[0].astype(str).isin(genes)]
    cats = terms.get(dbname)
    if cats:
        reference = reference[reference[1].str.contains("|".join(cats))]
    hits = set(reference[0])
    return pd.DataFrame([dict(name=name, db=dbname, counted_unit="source TF identities",
                              **summarize(set(network.source), hits)) for name, network in networks.items()])


def pair_counts(networks, genes, dbname):
    """Expose native overlap-enriched TF pairs and their reference agreements."""
    reference = gt.ds.read_db("hg38", dbname, verbose=False)
    reference = reference[reference[0].astype(str).isin(genes) & reference[1].astype(str).isin(genes)]
    reference_pairs = {"|".join(sorted([a, b])) for a, b in zip(reference[0], reference[1], strict=True)}
    candidate_tfs = set(reference[0]) | set(reference[1])
    rows = []
    for name, network in networks.items():
        predicted = _find_pairs(network[network.source.isin(candidate_tfs)], thr_pval=.01, verbose=False)
        rows.append(dict(name=name, db=dbname, counted_unit="unordered TF pairs",
                         **summarize(predicted, reference_pairs)))
    return pd.DataFrame(rows)


def knocktf_details(data, networks, terms):
    """List eligible experiments, target coverage and the native activity test."""
    reference = gt.ds.read_db("hg38", "KnockTF (scoring)", verbose=False)
    frames, support_sets = [], []
    for name, network in networks.items():
        measured = network[network.target.isin(reference.var_names)].copy()
        eligible, supported, _ = _knocktf_support(measured, data, reference, terms.get("KnockTF (scoring)"))
        support_sets.append(set(supported))
        rows = []
        for experiment in eligible:
            meta = reference.obs.loc[experiment]
            edges = measured[measured.source == meta.TF]
            values = reference[[experiment], edges.target.to_list()].X
            values = values.toarray() if sparse.issparse(values) else np.asarray(values)
            row = dict(name=name, experiment=experiment, TF=meta.TF, tissue=meta["Tissue.Type"],
                       knockdown_logFC=float(meta.logFC), n_measured_targets=len(edges),
                       n_nonzero_targets=int(np.count_nonzero(values)), tested=experiment in supported,
                       activity=np.nan, raw_pvalue=np.nan, adjusted_pvalue=np.nan, recovered=False)
            if experiment in supported:
                score, pvalue = dc.mt.ulm(reference[[experiment], :].to_df(), edges.rename(columns={"score": "weight"}), tmin=3, verbose=False)
                row["activity"] = float(score.iloc[0, 0])
                # One tested TF within each experiment: decoupler's within-row
                # adjustment is the identity. Native GRETA then BH-adjusts
                # across the supported experiments for this network.
                row["raw_pvalue"] = float(pvalue.iloc[0, 0])
            rows.append(row)
        frame = pd.DataFrame(rows)
        if len(frame) and frame.tested.any():
            frame.loc[frame.tested, "adjusted_pvalue"] = stats.false_discovery_control(frame.loc[frame.tested, "raw_pvalue"].to_numpy(), method="bh")
            frame["recovered"] = frame.tested & (frame.activity < 0) & (frame.adjusted_pvalue < .05)
        frames.append(frame)
    common = set.intersection(*support_sets) if support_sets else set()
    details = pd.concat(frames, ignore_index=True)
    if len(details):
        details["tested_in_every_network"] = details.experiment.isin(common)
    return details, len(common)


def comparable(results, coverage=None, min_knocktf_tests=10):
    keep, omitted = [], []
    methods = set(results.name)
    for dbname, group in results.groupby("db", sort=False):
        if dbname in MARKERS:
            reason = "Source-only diagnostic: the compared networks deliberately have the same TF set"
        elif (dbname == "KnockTF (scoring)" and coverage is not None
              and coverage.loc[coverage.db == dbname, "n_common_tested"].min() < min_knocktf_tests):
            n_common = int(coverage.loc[coverage.db == dbname, "n_common_tested"].min())
            reason = (f"Low-coverage diagnostic: {n_common} common tested experiments; "
                      f"minimum {min_knocktf_tests} for composite ranking (pragmatic teaching guard, not a universal statistical threshold)")
        elif set(group.name) != methods or group.f01.isna().any():
            reason = "Incomplete or undefined evaluation coverage; omit this reference for all networks"
        else:
            reason = ""
        if reason:
            omitted.append(group.assign(reason=reason))
        else:
            keep.append(group)
    ranked = pd.concat(keep, ignore_index=True) if keep else results.iloc[:0].copy()
    excluded = pd.concat(omitted, ignore_index=True) if omitted else results.iloc[:0].copy()
    return ranked, excluded


def evaluate(data, networks, metrics, terms, output, min_knocktf_tests=10):
    output.mkdir(parents=True, exist_ok=True)
    (output / "networks").mkdir(exist_ok=True)
    # Remove only named artifacts owned by this runner when a reduced rerun
    # omits a category; cached old diagnostics must not appear current.
    obsolete = []
    if "KnockTF (scoring)" not in metrics:
        obsolete.append("knocktf_experiment_details.tsv")
    if not set(metrics) & set(GENE_SETS):
        obsolete += ["gene_sets_background_sensitivity.tsv", "gene_sets_pathway_calls.tsv",
                     "03_gene_sets_background.png", "03_gene_sets_background.pdf"]
    if not set(metrics) & set(MARKERS):
        obsolete.append("tf_marker_counts.tsv")
    if not set(metrics) & {"IntAct", "Europe PMC"}:
        obsolete.append("tf_pair_counts.tsv")
    for filename in obsolete:
        (output / filename).unlink(missing_ok=True)
    frames, coverage = [], []
    for name, network in networks.items():
        print(f"  {name}: {len(network):,} edges, {network.source.nunique()} TFs", flush=True)
        frame, covered = evaluate_networks(data.copy(), {name: network}, terms, metrics, dataset_name="Heart cell-type pseudobulk")
        frames.append(frame)
        coverage.append(covered)
        network.to_csv(output / "networks" / f"{name}.tsv", sep="\t", index=False)
    results, coverage = pd.concat(frames, ignore_index=True), pd.concat(coverage, ignore_index=True)
    summaries, calls, counts, pair_tables, all_counts = [], [], [], [], []
    for dbname in metrics:
        if dbname in GENE_SETS:
            summary, call = pathway_audit(data, networks, dbname)
            # Audit against native GRETA before exporting any alternate background.
            native = summary.loc[summary.background == "native_20000"]
            for row in native.itertuples():
                actual = results.loc[(results.name == row.name) & (results.db == dbname)].iloc[0]
                for metric in ("precision", "recall", "f01"):
                    np.testing.assert_allclose(getattr(row, metric), actual[metric], rtol=1e-10, atol=1e-12)
                coverage.loc[(coverage.name == row.name) & (coverage.db == dbname), "n_tested"] = row.n_rna_tested_pathways
                coverage.loc[(coverage.name == row.name) & (coverage.db == dbname), "n_reference"] = row.n_expression_hits
                coverage.loc[(coverage.name == row.name) & (coverage.db == dbname), "n_predicted"] = row.n_predicted
                coverage.loc[(coverage.name == row.name) & (coverage.db == dbname), "n_hits"] = row.tp
                if row.n_expression_hits == 0:
                    results.loc[(results.name == row.name) & (results.db == dbname), ["precision", "recall", "f01"]] = np.nan
                    coverage.loc[(coverage.name == row.name) & (coverage.db == dbname), "status"] = "unsupported: no RNA-side pathway hits"
            summaries.append(summary)
            calls.append(call)
            all_counts.append(native.assign(n_reference=native.n_expression_hits))
        elif dbname in MARKERS:
            count = marker_counts(networks, set(data.var_names), terms, dbname)
            counts.append(count)
            all_counts.append(count.assign(n_reference=count.n_expression_hits))
            for row in count.itertuples():
                for metric in ("precision", "recall", "f01"):
                    actual = results.loc[(results.name == row.name) & (results.db == dbname), metric].iloc[0]
                    np.testing.assert_allclose(actual, getattr(row, metric), rtol=1e-10, atol=1e-12)
                mask = (coverage.name == row.name) & (coverage.db == dbname)
                coverage.loc[mask, "n_tested"] = row.n_predicted
                coverage.loc[mask, "n_predicted"] = row.n_predicted
                coverage.loc[mask, "n_reference"] = row.n_expression_hits
                coverage.loc[mask, "n_hits"] = row.tp
            coverage.loc[coverage.db == dbname, "status"] = "diagnostic: invariant source TF set"
        elif dbname in ("IntAct", "Europe PMC"):
            count = pair_counts(networks, set(data.var_names), dbname)
            pair_tables.append(count)
            all_counts.append(count.assign(n_reference=count.n_expression_hits))
            for row in count.itertuples():
                for metric in ("precision", "recall", "f01"):
                    actual = results.loc[(results.name == row.name) & (results.db == dbname), metric].iloc[0]
                    np.testing.assert_allclose(actual, getattr(row, metric), rtol=1e-10, atol=1e-12)
                mask = (coverage.name == row.name) & (coverage.db == dbname)
                candidate_ref = gt.ds.read_db("hg38", dbname, verbose=False)
                candidates = set(candidate_ref[0]) | set(candidate_ref[1])
                n_tf = networks[row.name].loc[networks[row.name].source.isin(candidates), "source"].nunique()
                coverage.loc[mask, "n_tested"] = n_tf * (n_tf - 1) // 2
                coverage.loc[mask, "n_predicted"] = row.n_predicted
                coverage.loc[mask, "n_reference"] = row.n_expression_hits
                coverage.loc[mask, "n_hits"] = row.tp
                if row.n_expression_hits == 0:
                    results.loc[(results.name == row.name) & (results.db == dbname), ["precision", "recall", "f01"]] = np.nan
                    coverage.loc[mask, "status"] = "unsupported: no measured reference TF pairs"
    if "KnockTF (scoring)" in metrics:
        details, n_common = knocktf_details(data, networks, terms)
        details.to_csv(output / "knocktf_experiment_details.tsv", sep="\t", index=False)
        coverage.loc[coverage.db == "KnockTF (scoring)", "n_common_tested"] = n_common
        knock_counts = []
        for row in coverage[coverage.db == "KnockTF (scoring)"].itertuples():
            actual = results.loc[(results.name == row.name) & (results.db == row.db)].iloc[0]
            tested = int(row.n_tested)
            eligible = int(row.n_eligible)
            hits = int(round(actual.precision * tested)) if tested else 0
            recovered = int(details.loc[details.name == row.name, "recovered"].sum()) if len(details) else 0
            if hits != recovered:
                raise RuntimeError(f"KnockTF experiment audit disagrees with native GRETA for {row.name}.")
            knock_counts.append(dict(name=row.name, db=row.db, counted_unit="external TF knockdown experiments",
                                     n_reference=eligible, n_predicted=tested, tp=hits, fp=tested-hits,
                                     fn=eligible-hits, precision=actual.precision, recall=actual.recall, f01=actual.f01))
        all_counts.append(pd.DataFrame(knock_counts))
    ranked, excluded = comparable(results, coverage, min_knocktf_tests)
    tables = dict(evaluation_all=results, evaluation_ranked=ranked, excluded_from_ranking=excluded,
                  coverage=coverage, network_stats=gt.tl.stats(networks), overlap=gt.tl.ocoeff(networks))
    if summaries:
        tables["gene_sets_background_sensitivity"] = pd.concat(summaries, ignore_index=True)
        tables["gene_sets_pathway_calls"] = pd.concat(calls, ignore_index=True)
    if counts:
        tables["tf_marker_counts"] = pd.concat(counts, ignore_index=True)
    if pair_tables:
        tables["tf_pair_counts"] = pd.concat(pair_tables, ignore_index=True)
    if all_counts:
        columns = ["name", "db", "counted_unit", "n_predicted", "n_reference", "tp", "fp", "fn", "precision", "recall", "f01"]
        tables["task_counts"] = pd.concat(all_counts, ignore_index=True)[columns]
    for name, table in tables.items():
        table.to_csv(output / f"{name}.tsv", sep="\t", index=False)
    # The overall native task heatmap is created first, then category details.
    if len(ranked):
        panels = [("01_ranking_task", ranked, (14, 4.5))]
        for category in ("Predictive", "Literature", "Mechanistic"):
            subset = ranked.loc[ranked["class"] == category]
            if len(subset):
                panels.append((f"02_category_{category.lower()}", subset, (max(6, subset.db.nunique() * 2), 4.5)))
        for filename, table, size in panels:
            figure = gt.pl.ranking(table, level="task", return_fig=True, figsize=size)
            for suffix in ("png", "pdf"):
                figure.savefig(output / f"{filename}.{suffix}", dpi=150, bbox_inches="tight")
            plt.close(figure)
    else:
        for suffix in ("png", "pdf"):
            (output / f"01_ranking_task.{suffix}").unlink(missing_ok=True)
    # A reduced rerun must not retain an obsolete mechanistic category figure.
    for category in ("Predictive", "Literature", "Mechanistic"):
        if not len(ranked.loc[ranked["class"] == category]):
            for suffix in ("png", "pdf"):
                (output / f"02_category_{category.lower()}.{suffix}").unlink(missing_ok=True)
    if summaries:
        plot_background_sensitivity(tables["gene_sets_background_sensitivity"], output)
    return results, coverage, ranked


def plot_background_sensitivity(table, output):
    databases = [name for name in GENE_SETS if name in set(table.db)]
    fig, axes = plt.subplots(1, len(databases), figsize=(4 * len(databases), 3.8), squeeze=False, layout="constrained")
    for ax, dbname in zip(axes.ravel(), databases, strict=True):
        subset = table[table.db == dbname]
        x = np.arange(len(METHODS))
        for offset, background, label, color in ((-.18, "native_20000", "Native ORA: 20,000 genes", "#557a95"),
                                                  (.18, "eligible_nonself", "Eligible genes, excluding TF", "#c88d49")):
            scores = subset[subset.background == background].set_index("name").reindex(METHODS).f01
            ax.bar(x + offset, scores, width=.36, label=label, color=color)
        ax.set(title=dbname, xticks=x, xticklabels=METHODS, ylim=(0, 1), ylabel="F₀.₁ (pathway labels)")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0, 0].legend(fontsize=8, loc="upper left", bbox_to_anchor=(0, -.18), frameon=False)
    fig.suptitle("Gene Sets sensitivity: pathway concordance depends on the ORA universe", fontsize=12)
    for suffix in ("png", "pdf"):
        fig.savefig(output / f"03_gene_sets_background.{suffix}", dpi=150, bbox_inches="tight")
    plt.close(fig)


def run(args):
    started = time.perf_counter()
    if args.min_knocktf_tests < 1:
        raise ValueError("--min-knocktf-tests must be at least 1 (default 10; lower explicitly for a tiny-coverage diagnostic ranking).")
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    data = load_data(args.data)
    full, matched, audit = prepare_networks(data, args.scenic, args.collectri, args.seed, args.min_targets)
    audit["conflicts"].to_csv(output / "scenicplus_sign_conflicts.tsv", sep="\t", index=False)
    audit["tf_degrees"].to_csv(output / "candidate_tf_degrees.tsv", sep="\t", index=False)
    audit["filtering"].to_csv(output / "network_filtering.tsv", sep="\t", index=False)
    cache = configure_cache(args.cache_dir)
    metrics, available = select_references(cache, quick=args.quick, download=args.download_resources)
    available.to_csv(output / "reference_resources.tsv", sep="\t", index=False)
    terms = gt.get_terms("Heart")
    print("Primary comparison: full networks on shared eligible source TFs", flush=True)
    results, coverage, ranked = evaluate(data, full, metrics, terms, output / "full", args.min_knocktf_tests)
    print("Sensitivity comparison: equal outgoing degrees per shared TF", flush=True)
    matched_results, _, _ = evaluate(data, matched, metrics, terms, output / "degree_matched", args.min_knocktf_tests)
    omitted = [dict(task="Omics (Gene ~ TFs)", reason="Six pooled cell types, one profile each: no valid native stratified train/test split or ten held-out observations; not scored as zero."),
               dict(task="Genomic/ATAC", reason="The supplied file is a TF-gene projection without regulatory regions or ATAC measurements."),
               dict(task="CollecTRI reference overlap", reason="CollecTRI is one compared network: self-reference is circular and excluded."),
               dict(task="TF Markers ranking", reason="HPA and TF-Marker scores use source identities only; these are identical by construction. Saved as diagnostics."),
               dict(task="Perturbation forecasting / Boolean simulation", reason="Omitted to keep this reduced comparison light; six type averages do not provide replicate-based validation.")]
    pd.DataFrame(omitted).to_csv(output / "omitted_tasks.tsv", sep="\t", index=False)
    info = dict(
        inputs={key: dict(path=str(Path(path).resolve()), sha256=digest(path))
                for key, path in (("pseudobulk", args.data), ("scenicplus", args.scenic), ("collectri", args.collectri))},
        packages={name: version(name) for name in ("gretapy", "decoupler", "numpy", "pandas", "anndata", "scipy")},
        code={str(Path(__file__).name): digest(__file__), "greta_benchmark.py": digest(ROOT / "GRN_comparison/greta_benchmark.py"),
              "diagnose_gene_sets.py": digest(ROOT / "GRN_comparison/diagnose_gene_sets.py")},
        dataset=dict(n_profiles=data.n_obs, n_genes=data.n_vars, cell_types=data.obs.celltype.astype(str).tolist(),
                     values="Nonnegative log1p-normalized selected-gene pseudobulk counts; not gene-wise Z-scores"),
        seed=args.seed, min_targets=args.min_targets,
        min_knocktf_tests=args.min_knocktf_tests,
        knocktf_ranking_gate="At least the specified number of external experiments must be testable in every compared network; otherwise native scores remain a low-coverage diagnostic. Default10 is a pragmatic teaching guard, not a universal statistical threshold.",
        preprocessing={key: value for key, value in audit.items() if not isinstance(value, pd.DataFrame)},
        conflict_policy="Drop every row of an opposing activator/repressor TF-gene pair; do not average or select a winner. Same-sign duplicates retain greatest strength.",
        signed_weights="SCENIC+ activator=+strength, repressor=-strength; n_topics is metadata only. CollecTRI supplied signed weights.",
        primary="All measured nonself edges of shared TFs having >=min_targets targets in both networks. Outgoing degrees differ between SCENIC+ and CollecTRI.",
        random="One illustrative seed; random targets preserve each SCENIC+ TF outgoing degree and its signed score multiset. Does not preserve target in-degree or establish a null distribution.",
        degree_sensitivity="Per TF use min(SCENIC+, CollecTRI) target count; strongest SCENIC+ strengths, seeded tie-sampling CollecTRI. Random exactly matches these outgoing degrees.",
        reference_terms=terms, metrics=metrics, ranked_databases=sorted(ranked.db.unique()),
        interpretation="Descriptive reference agreement, not validated edge accuracy or independent held-out performance. Supplied SCENIC+ may have used a larger original dataset. Shared genes are a preselected subset; excluded genes are unknown, not zero expression.",
        genesets="Do pathways whose genes are relatively highly expressed also appear as enriched labels somewhere among TF target sets? RNA hits: native ULM q<0.01 and positive score in >20% of profiles (>=2 of six). Network calls: ORA q<0.01 in any regulon. Native fixed 20,000-gene ORA retained; eligible-gene restricted memberships/background are a separate sensitivity, not a proven correction for HVG/TF selection.",
        scores="Native GRETA precision/recall/F0.1; F0.1 strongly favors precision. Each task counts different units. Unsupported tasks are NA or explicit omissions, not zero.",
        compatibility="Reuse tested GRETA 0.0.3 dispatcher fixes and explicit KnockTF measured-target coverage handling; native scores/formulas preserved.",
        elapsed_seconds=round(time.perf_counter() - started, 2),
    )
    (output / "run_info.json").write_text(json.dumps(info, indent=2) + "\n")
    print(results[["name", "db", "precision", "recall", "f01"]].to_string(index=False), flush=True)
    print(f"Saved primary plots first at {output / 'full/01_ranking_task.png'}", flush=True)
    print(f"Equal-degree sensitivity: {output / 'degree_matched/01_ranking_task.png'}", flush=True)
    print(f"Completed in {info['elapsed_seconds']:.1f} seconds; all outputs: {output}", flush=True)
    return info


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", type=Path, default=HERE / "data/heart_pseudobulk.h5ad")
    parser.add_argument("--scenic", type=Path, default=HERE / "data/scenicplus_grn.tsv")
    parser.add_argument("--collectri", type=Path, default=ROOT / "GRN_comparison/data/collectri.tsv")
    parser.add_argument("--output-dir", type=Path, default=HERE / "results/grn_comparison")
    parser.add_argument("--cache-dir", type=Path, default=ROOT / "GRN_comparison/gretapy_data")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--min-targets", type=int, default=3, help="Shared-source threshold, distinct from TF activity's 20-target setting")
    parser.add_argument("--min-knocktf-tests", type=int, default=10, help="Minimum common tested KnockTF experiments for ranking; lower explicitly to include tiny-coverage scores")
    parser.add_argument("--quick", action="store_true", help="Skip the external KnockTF task/resource")
    parser.add_argument("--download-resources", action="store_true", help="Download missing GRETA references; KnockTF is approximately 153 MB")
    args = parser.parse_args()
    try:
        run(args)
    except (FileNotFoundError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
