"""Python translation of the RNA/TF section of prepare_cosmos_inputs.R.

Inputs and exported matrices use genes/TFs in rows and cell lines in columns.
The local CollecTRI snapshot is reused; no resource download is needed.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as student_t


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RNA = ROOT / "data/RNA/RNA_log2_FPKM_cleaned_common.csv"
DEFAULT_REGULONS = ROOT / "GRN_comparison/data/collectri.tsv"
DEFAULT_OUTPUT = ROOT / "TF_activity/results"


def _check_matrix(frame: pd.DataFrame) -> None:
    if frame.empty:
        raise ValueError("The RNA matrix is empty.")
    if frame.index.has_duplicates or frame.columns.has_duplicates:
        raise ValueError("Gene names and cell-line names must be unique.")
    if frame.index.isna().any() or frame.columns.isna().any():
        raise ValueError("Gene names and cell-line names must not be missing.")
    if np.isinf(frame.to_numpy(dtype=float)).any():
        raise ValueError("The RNA matrix contains infinite values.")


def load_rna(path: str | Path = DEFAULT_RNA) -> pd.DataFrame:
    """Read the existing gene × cell-line log2(FPKM + 1) CSV, retaining NAs.

    This is the cleaned common-cell-line input used by the R script, rather
    than the uncensored RNA input used in the separate GRN benchmark.
    """
    rna = pd.read_csv(path, index_col=0)
    if rna.index.isna().any():
        raise ValueError("Gene names must not be missing.")
    rna.index = rna.index.astype(str)
    rna.index.name = "gene"
    rna = rna.apply(pd.to_numeric, errors="raise").astype(float)
    _check_matrix(rna)
    return rna


def scale_rna(rna: pd.DataFrame, max_missing: int = 20) -> pd.DataFrame:
    """Filter missingness, then z-score each gene across cell lines.

    Matches the R script's `rowSums(is.na(...)) <= 20`, `mean(na.rm=TRUE)`
    and sample `sd(na.rm=TRUE)` (ddof=1). No missing expression is imputed.
    Genes with zero/undefined SD become NA and are omitted per cell later.
    """
    _check_matrix(rna)
    if not isinstance(max_missing, int) or max_missing < 0:
        raise ValueError("max_missing must be a non-negative integer.")
    keep = rna.isna().sum(axis=1) <= max_missing
    retained = rna.loc[keep].astype(float)
    if retained.empty:
        raise ValueError("No genes remain after filtering missing values.")
    means = retained.mean(axis=1, skipna=True)
    sds = retained.std(axis=1, skipna=True, ddof=1)
    undefined = ~np.isfinite(sds) | (sds == 0)
    scaled = retained.sub(means, axis=0).div(sds.mask(undefined), axis=0)
    scaled.index.name = "gene"
    scaled.attrs["preprocessing"] = {
        "n_input_genes": int(rna.shape[0]),
        "n_removed_missing": int((~keep).sum()),
        "n_retained_genes": int(retained.shape[0]),
        "n_zero_or_undefined_sd": int(undefined.sum()),
        "n_cell_lines": int(rna.shape[1]),
        "max_missing": max_missing,
        "ddof": 1,
        "imputation": "none",
    }
    return scaled


def _check_regulons(regulons: pd.DataFrame) -> None:
    required = {"source", "target", "weight"}
    if not required.issubset(regulons.columns) or regulons.empty:
        raise ValueError("Regulons must contain source, target and weight columns.")
    if regulons[["source", "target", "weight"]].isna().any().any():
        raise ValueError("Regulons must not contain missing names or weights.")
    if regulons.duplicated(["source", "target"]).any():
        raise ValueError(
            "Duplicate source–target pairs must be resolved explicitly; "
            "weights are not averaged or silently discarded."
        )
    if not np.isfinite(regulons["weight"].to_numpy(dtype=float)).all():
        raise ValueError("Regulon weights must be finite.")


def load_regulons(path: str | Path = DEFAULT_REGULONS) -> pd.DataFrame:
    """Load the text export of the exact CollecTRI snapshot used in R.

    Original `mor` annotations are stored as `score` in the bundled TSV;
    decoupler 2.x expects the column name `weight`. Signs are preserved.
    The bundled snapshot has no duplicate source–target pairs.
    """
    regulons = pd.read_csv(path, sep="\t").rename(columns={"score": "weight"})
    _check_regulons(regulons)
    regulons = regulons.loc[:, ["source", "target", "weight"]].copy()
    regulons[["source", "target"]] = regulons[["source", "target"]].astype(str)
    regulons["weight"] = pd.to_numeric(regulons["weight"], errors="raise").astype(float)
    return regulons


def estimate_tf_activities(
    scaled_rna: pd.DataFrame,
    regulons: pd.DataFrame,
    min_targets: int = 20,
    verbose: bool = False,
) -> dict:
    """Run signed ULM independently in each cell line, matching the R loop.

    Each model uses ALL finite RNA z-scores in that cell line, with a zero
    predictor weight for genes outside the TF's regulon. Only TFs with at
    least `min_targets` available targets are tested. `empty=False` keeps
    genuine zero z-scores in both that background and target coverage.

    Returns `scores`, `p_values`, `p_adjusted` (TF × cell-line matrices;
    unavailable TF–cell combinations are NA), `target_counts` (all prior
    TFs × cell lines), `cell_summary`, and `metadata`. Target counts include
    measured zero z-scores and exclude missing values. Decoupler 2.x returns
    BH-adjusted p-values per cell; raw two-sided ULM p-values are recovered
    from the returned t-statistics using df = n_available_genes - 2.
    """
    import decoupler as dc

    _check_matrix(scaled_rna)
    _check_regulons(regulons)
    if not isinstance(min_targets, int) or min_targets < 1:
        raise ValueError("min_targets must be a positive integer.")
    sources = pd.Index(sorted(regulons["source"].unique()), name="TF")
    target_counts = pd.DataFrame(0, index=sources, columns=scaled_rna.columns, dtype=int)
    all_scores, all_pvalues, all_padj, summary = {}, {}, {}, []
    for cell in scaled_rna.columns:
        # Cell-wise omission is essential: missingness differs between lines.
        values = scaled_rna[cell].dropna()
        available = regulons.loc[regulons["target"].isin(values.index)]
        counts = available.groupby("source").size().reindex(sources, fill_value=0)
        target_counts.loc[:, cell] = counts
        eligible = counts.index[counts >= min_targets]
        summary.append({
            "cell_line": cell,
            "n_genes": int(values.size),
            "n_tfs_tested": int(eligible.size),
            "zero_zscore_genes": int((values == 0).sum()),
        })
        if eligible.empty:
            all_scores[cell] = pd.Series(dtype=float)
            all_pvalues[cell] = pd.Series(dtype=float)
            all_padj[cell] = pd.Series(dtype=float)
            continue
        if values.size < 3 or values.std(ddof=1) == 0:
            raise ValueError(f"Cell line {cell!r} has insufficient variation for ULM.")
        cell_frame = values.to_frame(name=cell).T
        scores, padj = dc.mt.ulm(
            data=cell_frame,
            net=regulons,
            tmin=min_targets,
            empty=False,
            verbose=verbose,
        )
        cell_scores = scores.loc[cell]
        if not np.isfinite(cell_scores).all():
            raise ValueError(f"ULM returned undefined scores for cell line {cell!r}.")
        all_scores[cell] = cell_scores
        all_pvalues[cell] = pd.Series(
            2 * student_t.sf(np.abs(cell_scores), df=values.size - 2),
            index=cell_scores.index,
        )
        all_padj[cell] = padj.loc[cell]
    result = {
        "scores": pd.DataFrame(all_scores).sort_index(),
        "p_values": pd.DataFrame(all_pvalues).sort_index(),
        "p_adjusted": pd.DataFrame(all_padj).sort_index(),
        "target_counts": target_counts,
        "cell_summary": pd.DataFrame(summary).set_index("cell_line"),
        "metadata": {
            "decoupler_version": dc.__version__,
            "method": "ULM; t-statistic of signed target-weight coefficient, including intercept",
            "min_targets": min_targets,
            "n_regulon_edges": int(regulons.shape[0]),
            "n_regulon_tfs": int(sources.size),
            "background": "all non-missing RNA z-scores in each cell line, including zeros",
            "degrees_of_freedom": "n_available_genes_in_cell_line - 2",
            "p_values": "raw two-sided Student t p-values reconstructed from ULM scores",
            "p_adjusted": "native decoupler BH correction across tested TFs separately within each cell line",
            "missing_scores": "TF has fewer than min_targets available targets in that cell line",
            "preprocessing": scaled_rna.attrs.get("preprocessing", {}),
        },
    }
    for key in ("scores", "p_values", "p_adjusted"):
        result[key].index.name = "TF"
    return result


def save_outputs(
    scaled_rna: pd.DataFrame,
    result: dict,
    output_dir: str | Path = DEFAULT_OUTPUT,
    provenance: dict | None = None,
) -> dict[str, Path]:
    """Write reusable CSVs and a JSON manifest; never modify the RData inputs."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    frames = {
        "rna_scaled": scaled_rna,
        "tf_scores": result["scores"],
        "tf_pvalues": result["p_values"],
        "tf_padj": result["p_adjusted"],
        "tf_target_counts": result["target_counts"],
        "cell_summary": result["cell_summary"],
    }
    paths = {}
    for name, frame in frames.items():
        path = output_dir / f"{name}.csv"
        frame.to_csv(path, na_rep="NA")
        paths[name] = path
    manifest = dict(result["metadata"])
    manifest["orientation"] = "genes/TFs in rows, cell lines in columns; cell_summary is cell lines in rows"
    manifest["packages"] = {
        name: importlib.metadata.version(name) for name in ("decoupler", "numpy", "pandas", "scipy")
    }
    manifest["provenance"] = provenance or {}
    manifest["outputs"] = {key: path.name for key, path in paths.items()}
    paths["run_info"] = output_dir / "run_info.json"
    paths["run_info"].write_text(json.dumps(manifest, indent=2) + "\n")
    return paths


def _file_info(path: Path) -> dict:
    return {"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rna", type=Path, default=DEFAULT_RNA)
    parser.add_argument("--regulons", type=Path, default=DEFAULT_REGULONS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--max-missing", type=int, default=20)
    parser.add_argument("--min-targets", type=int, default=20)
    args = parser.parse_args()
    rna = load_rna(args.rna)
    scaled = scale_rna(rna, max_missing=args.max_missing)
    regulons = load_regulons(args.regulons)
    result = estimate_tf_activities(scaled, regulons, min_targets=args.min_targets)
    paths = save_outputs(scaled, result, args.output_dir, provenance={
        "rna": _file_info(args.rna),
        "regulons": _file_info(args.regulons),
        "reference_workflow": "scripts/prepare_cosmos_inputs.R, RNA and TF section only",
    })
    print(f"Scored {result['scores'].shape[0]} TFs across {rna.shape[1]} cell lines.")
    print(f"Saved outputs to {paths['run_info'].parent.resolve()}")


if __name__ == "__main__":
    main()
