"""Prepare six cardiac cell-type pseudobulks and compare TF activities from two GRNs.

Run against the bundled compact h5ad, or regenerate it with --cells rna.h5ad.
Raw counts are summed per cell type, normalized over the supplied genes to
10,000 counts, log1p transformed, then gene-wise Z-scored across cell types.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "cosmos-heart-matplotlib"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(Path(tempfile.gettempdir()) / "cosmos-heart-numba"))

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import pearsonr, spearmanr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "scripts"))
from tf_activity_decoupler import estimate_tf_activities, load_regulons, save_outputs, scale_rna

DEFAULT_H5AD = HERE / "data/heart_pseudobulk.h5ad"
DEFAULT_REGULONS = ROOT / "GRN_comparison/data/collectri.tsv"
DEFAULT_SCENIC = HERE / "data/scenicplus_grn.tsv"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pseudobulk_counts(cells: ad.AnnData, groupby: str = "cell_type", target_sum: float = 10_000) -> ad.AnnData:
    """Sum observed counts, retaining genuine zeros and the supplied gene set.

    These pooled cell-type profiles are descriptive: donors are NOT treated as
    replicates and no differential-expression test is performed. X must contain
    nonnegative integer counts, not already normalized/logged expression.
    """
    if groupby not in cells.obs:
        raise ValueError(f"Missing grouping annotation: {groupby!r}")
    if not cells.obs_names.is_unique or not cells.var_names.is_unique:
        raise ValueError("Cell identifiers and gene names must be unique.")
    if cells.obs[groupby].isna().any() or cells.n_obs == 0 or cells.n_vars == 0:
        raise ValueError("Input must contain cells/genes and nonmissing group labels.")
    if not np.isfinite(target_sum) or target_sum <= 0:
        raise ValueError("target_sum must be positive and finite.")
    x = sparse.csr_matrix(cells.X, dtype=np.float64)
    if not np.isfinite(x.data).all() or (x.data < 0).any() or not np.equal(x.data, np.floor(x.data)).all():
        raise ValueError("Expected finite, nonnegative integer counts in X; do not sum logged expression.")
    labels = cells.obs[groupby].astype(str)
    groups = sorted(labels.unique())
    if len(groups) < 3:
        raise ValueError("At least three cell-type profiles are required for this Z-score tutorial.")
    counts, records, composition = [], [], []
    for group in groups:
        mask = labels.eq(group).to_numpy()
        total = np.asarray(x[mask].sum(axis=0)).ravel()
        if total.sum() <= 0:
            raise ValueError(f"Group {group!r} has no counts.")
        counts.append(total.astype(np.int64))
        record = {
            "cell_type": group,
            "celltype": group,  # GRETA's annotation name; observations are cell-type profiles.
            "n_cells": int(mask.sum()),
            "n_counts_retained_genes": int(total.sum()),
        }
        for annotation, key in (("donor", "n_donors"), ("sample", "n_samples")):
            if annotation in cells.obs:
                record[key] = int(cells.obs.loc[mask, annotation].nunique())
        records.append(record)
        if "donor" in cells.obs:
            for donor, n in cells.obs.loc[mask, "donor"].value_counts().items():
                if n:
                    composition.append({"cell_type": group, "donor": str(donor), "n_cells": int(n)})
    counts = np.stack(counts)
    # Float64 prevents rounding integer sums before normalization.
    normalized = counts.astype(float) / counts.sum(axis=1, keepdims=True) * target_sum
    log_expression = np.log1p(normalized).astype(np.float32)
    obs = pd.DataFrame(records).set_index("cell_type", drop=False)
    obs.index.name = "pseudobulk"
    var = cells.var.copy()
    var["n_cells_detected"] = np.asarray((x > 0).sum(axis=0)).ravel().astype(int)
    var["n_pseudobulks_detected"] = (counts > 0).sum(axis=0)
    result = ad.AnnData(X=log_expression, obs=obs, var=var)
    result.layers["counts"] = counts
    result.uns["pseudobulk"] = {
        "groupby": groupby,
        "aggregation": "sum of provided integer counts across all cells of each cell type",
        "input_n_cells": int(cells.n_obs),
        "input_n_genes": int(cells.n_vars),
        "input_count_total": int(x.sum()),
        "target_sum": float(target_sum),
        "X": "natural log1p of counts normalized to target_sum over the supplied gene subset",
        "normalization_background": "provided genes only; full-transcriptome library totals are unavailable",
        "replication": "one pooled profile per cell type; donor variation is not modeled",
        "organism_tissue": "human heart inferred from source gene IDs and supplied cell annotations",
        "input_assumption": "X contains counts; integer values verified, assay provenance was not supplied",
    }
    if composition:
        result.uns["donor_composition"] = pd.DataFrame(composition)
    # Independently guard aggregation, before any normalization or float32 conversion.
    np.testing.assert_array_equal(counts.sum(axis=0), np.asarray(x.sum(axis=0)).ravel())
    return result


def prepare_file(cells_path: Path, output_path: Path, groupby: str = "cell_type") -> ad.AnnData:
    cells = ad.read_h5ad(cells_path)
    pseudobulk = pseudobulk_counts(cells, groupby=groupby)
    pseudobulk.uns["source"] = {"filename": cells_path.name, "sha256": sha256(cells_path)}
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pseudobulk.write_h5ad(output_path, compression="gzip")
    provenance = {
        "source_filename": cells_path.name,
        "source_sha256": sha256(cells_path),
        "original_shape": list(cells.shape),
        "pseudobulk_file": output_path.name,
        "pseudobulk_shape": list(pseudobulk.shape),
        "pseudobulk_sha256": sha256(output_path),
        "processing": dict(pseudobulk.uns["pseudobulk"]),
        "source_citation": "not supplied; original cell annotations were retained in aggregate form",
    }
    (output_path.parent / "pseudobulk_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Prepared {output_path}: {pseudobulk.n_obs} profiles × {pseudobulk.n_vars:,} genes ({output_path.stat().st_size:,} bytes)")
    return pseudobulk


def validate_prepared(pseudobulk: ad.AnnData) -> None:
    if "counts" not in pseudobulk.layers or "pseudobulk" not in pseudobulk.uns:
        raise ValueError("Expected the prepared pseudobulk with counts layer and processing metadata; regenerate with --cells.")
    counts = np.asarray(pseudobulk.layers["counts"])
    if counts.shape != pseudobulk.shape or not np.isfinite(counts).all() or (counts < 0).any():
        raise ValueError("Invalid pseudobulk counts.")
    if not np.equal(counts, np.floor(counts)).all() or (counts.sum(axis=1) <= 0).any():
        raise ValueError("Pseudobulk counts must be integers with nonzero library totals.")
    target_sum = float(pseudobulk.uns["pseudobulk"]["target_sum"])
    expected = np.log1p(counts / counts.sum(axis=1, keepdims=True) * target_sum)
    np.testing.assert_allclose(pseudobulk.X, expected, rtol=2e-7, atol=2e-7)


def plot_activity(scores: pd.DataFrame, output: Path, network_name: str = "CollecTRI") -> None:
    shown = scores.std(axis=1).sort_values(ascending=False).head(20).index
    view = scores.loc[shown]
    limit = float(np.nanpercentile(np.abs(view), 98))
    cm = plt.get_cmap("RdBu_r").copy()
    cm.set_bad("#dddddd")
    fig, ax = plt.subplots(figsize=(10.2, 7.2), layout="constrained")
    im = ax.imshow(view, cmap=cm, vmin=-limit, vmax=limit, aspect="auto")
    labels = [s.replace(" ", "\n") for s in scores.columns]
    ax.set_xticks(range(len(labels)), labels, fontsize=9)
    ax.set_yticks(range(len(shown)), shown)
    ax.set_title(f"{network_name} TF activities: cell types relative to the six-profile panel", fontsize=12)
    fig.colorbar(im, ax=ax, label="ULM t-statistic", fraction=.035, pad=.02)
    fig.supxlabel("20 most variable TF scores; original t-statistics, without row rescaling. Colors clipped for display.", fontsize=9)
    for extension in ("png", "pdf"):
        fig.savefig(output / f"tf_activity_heatmap.{extension}", dpi=160, bbox_inches="tight")
    plt.close(fig)


def load_scenicplus(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Use the supplied sign and positive strength, resolving pairs explicitly.

    This is the full SCENIC+ regulon, not the smaller common-source network
    used by GRETA. A TF's own RNA is not required to estimate its activity.
    """
    raw = pd.read_csv(path, sep="\t")
    required = ["TF", "gene", "type", "score", "n_topics"]
    if not set(required).issubset(raw) or raw[required].isna().any().any():
        raise ValueError(f"SCENIC+ requires nonmissing {required} columns.")
    if not raw.type.isin(["activator", "repressor"]).all():
        raise ValueError("SCENIC+ type must be activator or repressor.")
    raw["score"] = pd.to_numeric(raw.score, errors="raise")
    raw["n_topics"] = pd.to_numeric(raw.n_topics, errors="raise")
    if not np.isfinite(raw.score).all() or (raw.score <= 0).any():
        raise ValueError("SCENIC+ strengths must be positive and finite.")
    if not np.isfinite(raw.n_topics).all():
        raise ValueError("SCENIC+ n_topics must be finite.")
    opposing = raw.groupby(["TF", "gene"])["type"].transform("nunique") > 1
    conflicts = raw.loc[opposing].sort_values(["TF", "gene", "type"]).reset_index(drop=True)
    clean = raw.loc[~opposing].sort_values(
        ["TF", "gene", "score", "n_topics"], ascending=[True, True, False, False], kind="stable")
    same_sign_duplicates = int(clean.duplicated(["TF", "gene"]).sum())
    clean = clean.drop_duplicates(["TF", "gene"])
    clean["weight"] = clean.score * clean.type.map({"activator": 1.0, "repressor": -1.0})
    regulons = clean.rename(columns={"TF": "source", "gene": "target"})[["source", "target", "weight"]].copy()
    audit = {
        "raw_rows": int(len(raw)), "raw_tfs": int(raw.TF.nunique()),
        "opposing_sign_pairs_excluded": int(conflicts[["TF", "gene"]].drop_duplicates().shape[0]),
        "opposing_sign_rows_excluded": int(len(conflicts)),
        "same_sign_duplicate_rows_excluded": same_sign_duplicates,
        "clean_edges": int(len(regulons)), "clean_tfs": int(regulons.source.nunique()),
        "sign_conflict_policy": "exclude every row of a source-target pair with opposing type labels",
        "same_sign_duplicate_policy": "keep strongest supplied score, then largest n_topics as tie breaker",
        "weight": "supplied positive score multiplied by activator +1 or repressor -1; n_topics is not a weight",
        "source_scope": "all supplied TF sources; a source's own transcript need not be measured",
        "self_edges": "retained as supplied, consistently with CollecTRI TF activity estimation",
        "strength_provenance": "exact statistical meaning and inference settings were not supplied",
    }
    return regulons.reset_index(drop=True), conflicts, audit


def run(pseudobulk_path: Path, regulon_path: Path, output: Path, min_targets: int = 20,
        network_name: str = "CollecTRI", scenicplus: bool = False) -> dict:
    pseudobulk = ad.read_h5ad(pseudobulk_path)
    validate_prepared(pseudobulk)
    # .X is genes in columns; the shared NCI60 helper expects genes in rows.
    rna = pd.DataFrame(np.asarray(pseudobulk.X).T, index=pseudobulk.var_names, columns=pseudobulk.obs_names)
    rna.index.name = "gene"
    scaled = scale_rna(rna, max_missing=0)  # sample SD (ddof=1), across the six cell types.
    output.mkdir(parents=True, exist_ok=True)
    if scenicplus:
        regulons, conflicts, cleaning = load_scenicplus(regulon_path)
        conflicts.to_csv(output / "scenicplus_sign_conflicts.tsv", sep="\t", index=False)
        (output / "regulon_cleaning.json").write_text(json.dumps(cleaning, indent=2) + "\n")
    else:
        regulons = load_regulons(regulon_path)
        cleaning = {"source_scope": "all prior TF sources; source RNA not required", "weight": "supplied signed weights"}
    regulons.to_csv(output / "regulons_cleaned.tsv", sep="\t", index=False)
    result = estimate_tf_activities(scaled, regulons, min_targets=min_targets)
    if result["scores"].empty:
        raise ValueError("No TFs meet the observed-target requirement.")
    result["cell_summary"].index.name = "cell_type"
    result["metadata"]["observation_semantics"] = "pooled cell types, not individual cell lines or donor replicates"
    result["metadata"]["network"] = network_name
    result["metadata"]["regulon_cleaning"] = cleaning
    result["metadata"]["normalization"] = dict(pseudobulk.uns["pseudobulk"])
    # read_h5ad returns NumPy scalar metadata; make the manifest JSON serializable.
    result["metadata"]["normalization"] = json.loads(json.dumps(result["metadata"]["normalization"], default=lambda x: x.item() if isinstance(x, np.generic) else x))
    paths = save_outputs(scaled, result, output, provenance={
        "pseudobulk": {"filename": pseudobulk_path.name, "sha256": sha256(pseudobulk_path)},
        "regulons": {"filename": regulon_path.name, "sha256": sha256(regulon_path)},
        "code": {"script_sha256": sha256(Path(__file__)), "helper_sha256": sha256(ROOT / "scripts/tf_activity_decoupler.py")},
        "cleaned_regulons": {"filename": "regulons_cleaned.tsv", "sha256": sha256(output / "regulons_cleaned.tsv")},
        "reference_workflow": "scripts/prepare_cosmos_inputs.R RNA/TF logic, with count normalization before scaling",
    })
    pseudobulk.obs.to_csv(output / "pseudobulk_summary.csv")
    if "donor_composition" in pseudobulk.uns:
        pseudobulk.uns["donor_composition"].to_csv(output / "donor_composition.csv", index=False)
    pd.DataFrame(pseudobulk.layers["counts"].T, index=pseudobulk.var_names,
                 columns=pseudobulk.obs_names).rename_axis("gene").to_csv(output / "pseudobulk_counts.csv")
    rna.to_csv(output / "rna_log_normalized.csv")
    plot_activity(result["scores"], output, network_name)
    print(f"{network_name}: estimated {int(result['scores'].notna().sum().sum()):,} activities: {result['scores'].shape[0]} TFs × {result['scores'].shape[1]} cell types.")
    print(f"All scores are unthresholded; minimum targets = {min_targets}. Outputs: {paths['run_info'].parent}")
    return result


def comparison_tables(collectri: dict, scenicplus: dict) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    """Align actual tested TFs, retaining unavailable values rather than zero-filling."""
    a, b = collectri["scores"], scenicplus["scores"]
    if set(a.columns) != set(b.columns):
        raise ValueError("TF activities must refer to identical pseudobulk cell types.")
    shared = a.index.intersection(b.index).sort_values()
    if shared.empty:
        raise ValueError("No shared scored TFs are available for comparison.")
    a, b = a.loc[shared].copy(), b.reindex(index=shared, columns=a.columns).copy()
    delta = b - a
    tables = {"collectri_shared_scores": a, "scenicplus_shared_scores": b, "delta_scenicplus_minus_collectri": delta}
    records = []
    for cell in a.columns:
        available = a[cell].notna() & b[cell].notna()
        x, y = a.loc[available, cell], b.loc[available, cell]
        variable = len(x) >= 3 and x.nunique() > 1 and y.nunique() > 1
        nonzero = (x != 0) & (y != 0)
        records.append({
            "cell_type": cell, "n_shared_tested_tfs": int(len(x)),
            "pearson_r": float(pearsonr(x, y).statistic) if variable else np.nan,
            "spearman_rho": float(spearmanr(x, y).statistic) if variable else np.nan,
            "n_nonzero_pairs": int(nonzero.sum()),
            "sign_concordance_nonzero": float((np.sign(x[nonzero]) == np.sign(y[nonzero])).mean()) if nonzero.any() else np.nan,
            "median_delta_scenicplus_minus_collectri": float((y - x).median()),
            "median_absolute_delta": float((y - x).abs().median()),
        })
    return tables, pd.DataFrame(records).set_index("cell_type")


def draw_comparison_heatmap(ax, view: pd.DataFrame, limit: float, title: str):
    cm = plt.get_cmap("RdBu_r").copy()
    cm.set_bad("#dddddd")
    image = ax.imshow(view, cmap=cm, vmin=-limit, vmax=limit, aspect="auto")
    ax.set_xticks(range(view.shape[1]), [s.replace(" ", "\n") for s in view.columns], fontsize=8)
    ax.set_yticks(range(view.shape[0]), view.index, fontsize=9)
    ax.set_title(title, fontsize=12, pad=10)
    return image


def compare_activities(collectri: dict, scenicplus: dict, output: Path,
                       collectri_dir: Path, scenicplus_dir: Path) -> dict:
    """Export matched views and descriptive agreement, without testing method superiority.

    Top rows are selected by mean within-method variance across cell types,
    not by score differences. Different target sets/weights influence t-score
    magnitude; subtraction is descriptive and is not a statistical test.
    """
    output.mkdir(parents=True, exist_ok=True)
    tables, agreement = comparison_tables(collectri, scenicplus)
    a, b, delta = (tables[key] for key in ("collectri_shared_scores", "scenicplus_shared_scores", "delta_scenicplus_minus_collectri"))
    selection = ((a.var(axis=1, ddof=1) + b.var(axis=1, ddof=1)) / 2).rename("mean_within_network_variance")
    selection = selection.sort_values(ascending=False, kind="stable")
    shown = selection.head(20).index
    for name, frame in tables.items():
        frame.to_csv(output / f"{name}.csv", na_rep="NA")
    agreement.to_csv(output / "agreement_by_cell_type.csv", na_rep="NA")
    selection.to_frame().assign(shown=lambda x: x.index.isin(shown)).to_csv(output / "heatmap_tf_selection.csv")
    coverage = []
    for name, result in (("CollecTRI", collectri), ("SCENIC+", scenicplus)):
        frame = result["target_counts"].copy()
        frame.index.name = "TF"
        frame = frame.reset_index().melt(id_vars="TF", var_name="cell_type", value_name="observed_targets")
        frame.insert(0, "network", name)
        frame["score_available"] = [tf in result["scores"].index and pd.notna(result["scores"].at[tf, cell])
                                    for tf, cell in zip(frame.TF, frame.cell_type, strict=True)]
        frame["shared_scored_tf"] = frame.TF.isin(a.index)
        coverage.append(frame)
    pd.concat(coverage, ignore_index=True).to_csv(output / "target_coverage_all_tfs.csv", index=False)
    scope = pd.DataFrame(index=collectri["target_counts"].index.union(scenicplus["target_counts"].index).sort_values())
    scope.index.name = "TF"
    scope["collectri_scored"] = scope.index.isin(collectri["scores"].index)
    scope["scenicplus_scored"] = scope.index.isin(scenicplus["scores"].index)
    scope["shared_scored_tf"] = scope.index.isin(a.index)
    scope.to_csv(output / "tf_scopes.csv")
    views = [a.loc[shown], b.loc[shown], delta.loc[shown]]
    pooled = np.concatenate([views[0].to_numpy().ravel(), views[1].to_numpy().ravel()])
    score_limit = max(float(np.nanpercentile(np.abs(pooled), 98)), 1e-12)
    delta_limit = max(float(np.nanpercentile(np.abs(views[2].to_numpy()), 98)), 1e-12)
    titles = ["CollecTRI", "SCENIC+", "SCENIC+ − CollecTRI"]
    names = ["01_collectri_matched_heatmap", "02_scenicplus_matched_heatmap", "03_delta_heatmap"]
    for i, (view, title, name) in enumerate(zip(views, titles, names, strict=True)):
        fig, ax = plt.subplots(figsize=(10.2, 7.2), layout="constrained")
        im = draw_comparison_heatmap(ax, view, delta_limit if i == 2 else score_limit, title)
        fig.colorbar(im, ax=ax, label="Difference in ULM t-statistics" if i == 2 else "ULM t-statistic", fraction=.035, pad=.02)
        fig.supxlabel("Same 20 shared TFs; selected by mean within-network variance across cell types. No row rescaling.", fontsize=9)
        for extension in ("png", "pdf"):
            fig.savefig(output / f"{name}.{extension}", dpi=160, bbox_inches="tight")
        plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(18, 8.2), layout="constrained")
    for i, (ax, view, title) in enumerate(zip(axes, views, titles, strict=True)):
        im = draw_comparison_heatmap(ax, view, delta_limit if i == 2 else score_limit, title)
        fig.colorbar(im, ax=ax, label="Δ ULM t" if i == 2 else "ULM t", fraction=.04, pad=.02, shrink=.7)
    fig.suptitle("TF activities inferred from the same six RNA pseudobulks", fontsize=15)
    fig.supxlabel("Matched TF rows and cell-type columns. CollecTRI and SCENIC+ share a color scale; Δ has its own scale.\n"
                  "SCENIC+ − CollecTRI is descriptive; a larger difference does not establish better inference or statistical significance.", fontsize=10)
    for extension in ("png", "pdf"):
        fig.savefig(output / f"04_tf_activity_comparison.{extension}", dpi=160, bbox_inches="tight")
    plt.close(fig)
    manifest = {
        "n_collectri_scored_tfs": int(collectri["scores"].shape[0]),
        "n_scenicplus_scored_tfs": int(scenicplus["scores"].shape[0]),
        "n_shared_scored_tfs": int(a.shape[0]), "n_cell_types": int(a.shape[1]),
        "delta": "SCENIC+ ULM t-statistic minus CollecTRI ULM t-statistic; NA remains NA",
        "row_selection": "top 20 shared TFs by mean of each network's score variance across the same six cell types",
        "heatmap_tfs": shown.tolist(), "score_color_limit": score_limit, "delta_color_limit": delta_limit,
        "display": "original unscaled ULM t-statistics; symmetric color ranges clipped at 98th percentile; gray is missing",
        "interpretation": "descriptive agreement only; different target membership, degree, and relative weights influence score magnitude; no method is treated as ground truth",
        "experimental_units": "six pooled cell types, with no donor replicates or independent validation",
        "network_input_to_ulm": "full cleaned prior; source transcript need not be measured; all finite RNA genes form the common regression background",
        "min_targets": collectri["metadata"]["min_targets"],
        "scenicplus_min_targets": scenicplus["metadata"]["min_targets"],
        "provenance": {"code_sha256": sha256(Path(__file__)),
                       "collectri_manifest_sha256": sha256(collectri_dir / "run_info.json"),
                       "scenicplus_manifest_sha256": sha256(scenicplus_dir / "run_info.json"),
                       "collectri_rna_sha256": sha256(collectri_dir / "rna_scaled.csv"),
                       "scenicplus_rna_sha256": sha256(scenicplus_dir / "rna_scaled.csv")},
    }
    if manifest["provenance"]["collectri_rna_sha256"] != manifest["provenance"]["scenicplus_rna_sha256"]:
        raise ValueError("Network comparison must use the identical scaled RNA matrix.")
    (output / "run_info.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Compared {a.shape[0]} shared scored TFs across {a.shape[1]} cell types. Outputs: {output}")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pseudobulk", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--cells", type=Path, help="Optional original count h5ad; regenerate --pseudobulk before scoring")
    parser.add_argument("--groupby", default="cell_type")
    parser.add_argument("--network", choices=["both", "collectri", "scenicplus", "custom"], default=None,
                        help="Default: both networks. An explicit --regulons retains the original custom-network behavior.")
    parser.add_argument("--regulons", type=Path, default=None, help="Custom source/target/score or weight TSV; implies --network custom")
    parser.add_argument("--scenicplus", type=Path, default=DEFAULT_SCENIC)
    parser.add_argument("--output-dir", type=Path, default=HERE / "results/tf_activity")
    parser.add_argument("--scenicplus-output-dir", type=Path, default=HERE / "results/tf_activity_scenicplus")
    parser.add_argument("--comparison-dir", type=Path, default=HERE / "results/tf_activity_comparison")
    parser.add_argument("--min-targets", type=int, default=20)
    parser.add_argument("--prepare-only", action="store_true", help="Only regenerate the compact h5ad; requires --cells")
    args = parser.parse_args()
    if args.prepare_only and args.cells is None:
        parser.error("--prepare-only requires --cells")
    network = args.network or ("custom" if args.regulons else "both")
    if args.regulons and network != "custom":
        parser.error("--regulons is for --network custom; use --scenicplus for a SCENIC+ input")
    if network == "custom" and args.regulons is None:
        parser.error("--network custom requires --regulons")
    if args.cells is not None:
        prepare_file(args.cells, args.pseudobulk, args.groupby)
    if not args.prepare_only:
        if network == "custom":
            run(args.pseudobulk, args.regulons, args.output_dir, args.min_targets, network_name="Custom GRN")
            return
        if network in ("both", "collectri"):
            collectri = run(args.pseudobulk, DEFAULT_REGULONS, args.output_dir, args.min_targets)
        if network in ("both", "scenicplus"):
            scenic = run(args.pseudobulk, args.scenicplus, args.scenicplus_output_dir, args.min_targets,
                         network_name="SCENIC+", scenicplus=True)
        if network == "both":
            compare_activities(collectri, scenic, args.comparison_dir, args.output_dir, args.scenicplus_output_dir)


if __name__ == "__main__":
    main()
