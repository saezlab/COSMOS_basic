"""Audit Gene Sets scores for the saved default NCI60 RNA practical.

Run ``python diagnose_gene_sets.py`` after the default RNA benchmark. This
recomputes only Gene Sets diagnostics, not the full benchmark. It reproduces
native scores, changes the ORA universe in a separate sensitivity analysis,
and evaluates ten additional Hallmarks random controls. Existing benchmark
scores, networks and figures are read without modification.

The expression-side pathway hits, F0.1 formula, per-regulon correction and
union across TFs follow gretapy 0.0.3. The alternate universe excludes genes
that could not have been selected, including each source TF itself. This
adjustment does not turn pathway concordance into validated edge accuracy.
"""

from dataclasses import asdict
import json
from pathlib import Path

import decoupler as dc
import gretapy as gt
import numpy as np
import pandas as pd
from scipy import stats

from greta_benchmark import configure_cache
from rna_benchmark import RnaConfig, load_rna


HERE = Path(__file__).resolve().parent
INPUT = HERE / "results/rna"
OUTPUT = HERE / "results/gene_set_audit"
NETWORK_NAMES = ("Pearson", "CollecTRI", "Random")
DATABASES = ("Hallmarks", "PROGENy", "Reactome")


def validate_network(network, genes, name):
    """Require the nonself, measured-gene networks assumed by this audit."""
    if not {"source", "target", "score"}.issubset(network):
        raise ValueError(f"{name}: expected source, target and score columns.")
    if network[["source", "target", "score"]].isna().any().any():
        raise ValueError(f"{name}: identifiers and scores must be nonmissing.")
    if not np.isfinite(network.score.to_numpy(dtype=float)).all():
        raise ValueError(f"{name}: scores must be finite.")
    if network.duplicated(["source", "target"]).any():
        raise ValueError(f"{name}: duplicate edges are not supported.")
    if (network.source == network.target).any():
        raise ValueError(f"{name}: this audit assumes self edges were excluded.")
    if (set(network.source) | set(network.target)) - genes:
        raise ValueError(f"{name}: sources and targets must be retained NCI60 genes.")


def infer_pathways(network, pathway_targets, n_background, *, exclude_self):
    """Union pathways enriched in any regulon, matching native thresholds.

    The alternate pathway sets must already be restricted to eligible genes.
    With exclude_self=True, remove each TF from its own background and from
    every pathway's size. Validated feature lists already exclude that TF.
    """
    predicted = set()
    for tf, edges in network.groupby("source"):
        features = set(edges.target)
        overlaps = {pw: len(features & targets) for pw, targets in pathway_targets.items()}
        relevant = [pw for pw, overlap in overlaps.items() if overlap > 0]
        if not relevant:
            continue
        actual_background = n_background - int(exclude_self)
        # Hypergeometric upper tails equal one-sided Fisher exact p values.
        pvals = np.array([
            stats.hypergeom.sf(
                overlaps[pw] - 1, actual_background,
                len(pathway_targets[pw]) - int(exclude_self and tf in pathway_targets[pw]),
                len(features),
            )
            for pw in relevant
        ])
        adjusted = stats.false_discovery_control(pvals, method="bh")
        predicted.update(pw for pw, q in zip(relevant, adjusted) if q < 0.01)
    return predicted


def summarize(predicted, hits):
    """Count agreement with expression-side hits, not biological edge truth."""
    tp, fp, fn = len(predicted & hits), len(predicted - hits), len(hits - predicted)
    precision = tp / (tp + fp) if tp else 0.0
    recall = tp / (tp + fn) if tp else 0.0
    f01 = 1.01 * precision * recall / (0.01 * precision + recall) if tp else 0.0
    return dict(n_predicted=len(predicted), n_expression_hits=len(hits),
                tp=tp, fp=fp, fn=fn, precision=precision, recall=recall, f01=f01)


def verify_fisher_equivalence(n_genes):
    checks = []
    for background, size, n_features, overlap in (
        (20000, 200, 30, 3), (n_genes - 1, 87, 34, 4), (n_genes - 1, 50, 3, 1),
    ):
        hypergeom = stats.hypergeom.sf(overlap - 1, background, size, n_features)
        fisher = stats.fisher_exact(
            [[overlap, size - overlap],
             [n_features - overlap, background - size - n_features + overlap]],
            alternative="greater",
        ).pvalue
        np.testing.assert_allclose(hypergeom, fisher, rtol=1e-12, atol=0)
        checks.append(dict(N=background, M=size, n=n_features, a=overlap,
                           hypergeom=float(hypergeom), fisher=float(fisher)))
    return checks


def random_network(genes, degrees, seed):
    """Match each TF's outgoing edge count; sample random targets and signs."""
    rng = np.random.default_rng(np.random.SeedSequence(seed))
    rows = []
    for tf in sorted(degrees):
        possible = genes[genes != tf]
        for target in rng.choice(possible, size=degrees[tf], replace=False):
            rows.append((tf, str(target), float(rng.choice([-1, 1]))))
    return pd.DataFrame(rows, columns=["source", "target", "score"])


def main():
    """Read default NCI60 results and write separate audit JSON and TSV files."""
    config = RnaConfig()
    info_path = INPUT / "run_info.json"
    if not info_path.exists():
        raise FileNotFoundError("Run python run_rna_benchmark.py first using default settings.")
    run_info = json.loads(info_path.read_text())
    if run_info.get("config") != asdict(config) or run_info.get("manifest", {}).get("terms") != {}:
        raise ValueError("This diagnostic supports the default NCI60 configuration and unfiltered references only.")
    configure_cache()
    adata, _, _ = load_rna(config)
    genes = set(adata.var_names)
    networks = {
        name: pd.read_csv(INPUT / "networks" / f"{name}.tsv", sep="\t")
        for name in NETWORK_NAMES
    }
    degrees = networks["CollecTRI"].groupby("source").size().to_dict()
    for name, network in networks.items():
        validate_network(network, genes, name)
        if network.groupby("source").size().to_dict() != degrees:
            raise ValueError(f"{name}: expected the same source TFs and outgoing degrees as CollecTRI.")
    saved = pd.read_csv(INPUT / "evaluation_all.tsv", sep="\t")
    audit = dict(
        scope="Saved default NCI60 RNA practical only; no full benchmark rerun",
        n_genes=len(genes), n_tfs=len(degrees), n_candidate_targets_per_tf=len(genes) - 1,
        background_modes={
            "native20k": "Unfiltered reference pathways and a fixed 20,000-gene background",
            "eligible_nonself": "Pathways restricted to eligible RNA genes; each source TF excluded",
        },
        native_scores_verified=0, databases={},
        fisher_verification=verify_fisher_equivalence(len(genes)),
    )
    rows = []
    for dbname in DATABASES:
        print(f"Auditing {dbname}...", flush=True)
        db = gt.ds.read_db("hg38", dbname, verbose=False)
        dc.mt.ulm(adata, db, verbose=False)
        scores = adata.obsm["score_ulm"]
        pvals = 2 * stats.t.sf(abs(scores.values), adata.n_vars - 2)
        adjusted = stats.false_discovery_control(pvals.ravel(), method="bh").reshape(pvals.shape)
        hits = set(scores.columns[((adjusted < 0.01) & (scores.values > 0)).mean(0) > 0.2])
        pathway_targets = {source: set(targets) for source, targets in db.groupby("source").target}
        backgrounds = {
            "native20k": (20000, pathway_targets, False),
            "eligible_nonself": (len(genes), {s: t & genes for s, t in pathway_targets.items()}, True),
        }
        results = dict(
            n_pathways=len(pathway_targets), n_expression_hits=len(hits),
            reference_total_genes=int(db.target.nunique()),
            reference_expressed_genes=len(set(db.target) & genes), scores={},
        )
        for background, (n_bg, pathways, exclude_self) in backgrounds.items():
            results["scores"][background] = {}
            for name, network in networks.items():
                result = summarize(infer_pathways(network, pathways, n_bg, exclude_self=exclude_self), hits)
                if background == "native20k":
                    expected = saved.loc[(saved.name == name) & (saved.db == dbname)]
                    if len(expected) != 1:
                        raise ValueError(f"Expected one saved native score for {name}/{dbname}.")
                    for metric in ("precision", "recall", "f01"):
                        np.testing.assert_allclose(result[metric], expected.iloc[0][metric], rtol=1e-12, atol=1e-14)
                    audit["native_scores_verified"] += 1
                results["scores"][background][name] = result
                rows.append(dict(db=dbname, background=background, name=name,
                                 seed=config.seed if name == "Random" else "", **result))
        if dbname == "Hallmarks":
            results["random_seeds"] = {}
            for seed in range(10):
                network = random_network(adata.var_names.to_numpy(), degrees, seed)
                results["random_seeds"][seed] = {}
                for background, (n_bg, pathways, exclude_self) in backgrounds.items():
                    result = summarize(infer_pathways(network, pathways, n_bg, exclude_self=exclude_self), hits)
                    results["random_seeds"][seed][background] = result
                    rows.append(dict(db=dbname, background=background, name="Random", seed=seed, **result))
        audit["databases"][dbname] = results
    OUTPUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUTPUT / "scores.tsv", sep="\t", index=False)
    (OUTPUT / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(f"Verified all {audit['native_scores_verified']} native scores; wrote {len(rows)} diagnostic rows to {OUTPUT}")
    return audit


if __name__ == "__main__":
    main()
