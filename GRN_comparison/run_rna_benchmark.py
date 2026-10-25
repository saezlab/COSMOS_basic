"""Run the full expressed-gene, RNA-only NCI60 GRETA teaching benchmark."""

import argparse
import json
import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "cosmos-grn-matplotlib"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(Path(tempfile.gettempdir()) / "cosmos-grn-numba"))
import matplotlib
matplotlib.use("Agg")
import gretapy as gt
from rna_benchmark import RNA_PATH, SAMPLES_PATH, RnaConfig, run_rna


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=["teaching", "quick"], default="teaching",
                        help="quick omits the external KnockTF task; both retain full RNA networks")
    parser.add_argument("--rna", type=Path, default=RNA_PATH, help="Gene-by-sample TSV(.gz) or prepared AnnData(.h5ad)")
    parser.add_argument("--samples", type=Path, default=SAMPLES_PATH,
                        help="TSV with sample,tissue_code columns, for TSV RNA inputs")
    parser.add_argument("--terms-json", type=Path, help="JSON mapping GRETA database names to context-term lists")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--include-spearman", action="store_true")
    parser.add_argument("--max-genes", type=int, default=0, help="Optional training-variance cap; 0 keeps all eligible genes")
    parser.add_argument("--max-tfs", type=int, default=0, help="Optional training-variance TF cap; 0 keeps all eligible TFs")
    parser.add_argument("--min-expression", type=float, default=1.0)
    parser.add_argument("--min-sample-fraction", type=float, default=0.8)
    parser.add_argument("--seed", type=int, default=42, help="Random-network seed; native evaluation split remains seed 42")
    parser.add_argument("--grn", action="append", default=[], metavar="NAME=PATH",
                        help="Optional existing TF-target-score network; genes must occur in retained RNA")
    parser.add_argument("--force", action="store_true", help="Recompute even when a matching verified cache exists")
    args = parser.parse_args()
    try:
        terms = json.loads(args.terms_json.read_text()) if args.terms_json else None
        config = RnaConfig(
            seed=args.seed, min_expression=args.min_expression,
            min_sample_fraction=args.min_sample_fraction,
            max_genes=args.max_genes, max_tfs=args.max_tfs,
            include_spearman=args.include_spearman, profile=args.profile,
        )
        grn_paths = dict(item.split("=", 1) for item in args.grn)
        result = run_rna(config, output_dir=args.output_dir, rna_path=args.rna,
                         samples_path=args.samples, terms=terms, grn_paths=grn_paths,
                         force=args.force)
        print(gt.tl.stats(result["networks"]).to_string(index=False))
        print(gt.tl.ranking(result["res"]).to_string())
        mode = "Verified cached scores" if result["cache_hit"] else "Computed native GRETA scores"
        print(f"{mode}; {result['adata'].n_obs} samples and {result['adata'].n_vars} retained genes.")
        print(f"Saved to {result['output_dir']} ({result['elapsed_seconds']:.1f} seconds excluding imports)")
    except (ValueError, FileNotFoundError, json.JSONDecodeError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
