"""Prepare six cardiac cell-type pseudobulks and estimate CollecTRI TF activities.

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

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "scripts"))
from tf_activity_decoupler import estimate_tf_activities, load_regulons, save_outputs, scale_rna

DEFAULT_H5AD = HERE / "data/heart_pseudobulk.h5ad"
DEFAULT_REGULONS = ROOT / "GRN_comparison/data/collectri.tsv"


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


def plot_activity(scores: pd.DataFrame, output: Path) -> None:
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
    ax.set_title("CollecTRI TF activities: cell types relative to the six-profile panel", fontsize=12)
    fig.colorbar(im, ax=ax, label="ULM t-statistic", fraction=.035, pad=.02)
    fig.supxlabel("20 most variable TF scores; original t-statistics, without row rescaling. Colors clipped for display.", fontsize=9)
    for extension in ("png", "pdf"):
        fig.savefig(output / f"tf_activity_heatmap.{extension}", dpi=160, bbox_inches="tight")
    plt.close(fig)


def run(pseudobulk_path: Path, regulon_path: Path, output: Path, min_targets: int = 20) -> dict:
    pseudobulk = ad.read_h5ad(pseudobulk_path)
    validate_prepared(pseudobulk)
    # .X is genes in columns; the shared NCI60 helper expects genes in rows.
    rna = pd.DataFrame(np.asarray(pseudobulk.X).T, index=pseudobulk.var_names, columns=pseudobulk.obs_names)
    rna.index.name = "gene"
    scaled = scale_rna(rna, max_missing=0)  # sample SD (ddof=1), across the six cell types.
    regulons = load_regulons(regulon_path)
    result = estimate_tf_activities(scaled, regulons, min_targets=min_targets)
    if result["scores"].empty:
        raise ValueError("No TFs meet the observed-target requirement.")
    result["cell_summary"].index.name = "cell_type"
    result["metadata"]["observation_semantics"] = "pooled cell types, not individual cell lines or donor replicates"
    result["metadata"]["normalization"] = dict(pseudobulk.uns["pseudobulk"])
    # read_h5ad returns NumPy scalar metadata; make the manifest JSON serializable.
    result["metadata"]["normalization"] = json.loads(json.dumps(result["metadata"]["normalization"], default=lambda x: x.item() if isinstance(x, np.generic) else x))
    paths = save_outputs(scaled, result, output, provenance={
        "pseudobulk": {"filename": pseudobulk_path.name, "sha256": sha256(pseudobulk_path)},
        "regulons": {"filename": regulon_path.name, "sha256": sha256(regulon_path)},
        "reference_workflow": "scripts/prepare_cosmos_inputs.R RNA/TF logic, with count normalization before scaling",
    })
    pseudobulk.obs.to_csv(output / "pseudobulk_summary.csv")
    if "donor_composition" in pseudobulk.uns:
        pseudobulk.uns["donor_composition"].to_csv(output / "donor_composition.csv", index=False)
    pd.DataFrame(pseudobulk.layers["counts"].T, index=pseudobulk.var_names,
                 columns=pseudobulk.obs_names).rename_axis("gene").to_csv(output / "pseudobulk_counts.csv")
    rna.to_csv(output / "rna_log_normalized.csv")
    plot_activity(result["scores"], output)
    print(f"Estimated {int(result['scores'].notna().sum().sum()):,} activities: {result['scores'].shape[0]} TFs × {result['scores'].shape[1]} cell types.")
    print(f"All scores are unthresholded; minimum targets = {min_targets}. Outputs: {paths['run_info'].parent}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pseudobulk", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--cells", type=Path, help="Optional original count h5ad; regenerate --pseudobulk before scoring")
    parser.add_argument("--groupby", default="cell_type")
    parser.add_argument("--regulons", type=Path, default=DEFAULT_REGULONS)
    parser.add_argument("--output-dir", type=Path, default=HERE / "results/tf_activity")
    parser.add_argument("--min-targets", type=int, default=20)
    parser.add_argument("--prepare-only", action="store_true", help="Only regenerate the compact h5ad; requires --cells")
    args = parser.parse_args()
    if args.prepare_only and args.cells is None:
        parser.error("--prepare-only requires --cells")
    if args.cells is not None:
        prepare_file(args.cells, args.pseudobulk, args.groupby)
    if not args.prepare_only:
        run(args.pseudobulk, args.regulons, args.output_dir, args.min_targets)


if __name__ == "__main__":
    main()
