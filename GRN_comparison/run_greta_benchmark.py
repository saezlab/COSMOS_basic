"""Run a reduced GRETA task benchmark and the upstream task-ranking plot."""

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
from greta_benchmark import GretaConfig, configure_cache, published_results, run_greta


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=["quick", "teaching"], default="teaching")
    parser.add_argument("--n-cells", type=int, default=180)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--edges-per-tf", type=int, default=5)
    parser.add_argument("--output-dir")
    parser.add_argument("--data", help="Your prepared human hg38 MuData (.h5mu)")
    parser.add_argument("--grn", action="append", default=[], metavar="NAME=PATH")
    parser.add_argument("--context", help="GRETA terms context, e.g. PBMC; required for custom data unless --terms-json")
    parser.add_argument("--terms-json", help="JSON mapping reference database names to relevant context terms")
    parser.add_argument("--metrics", nargs="+", help="Explicit reference/subtask names, overriding the profile")
    parser.add_argument("--published-only", action="store_true", help="Plot stored published results instead of recomputing metrics")
    args = parser.parse_args()
    try:
        configure_cache()
        if args.published_only:
            output = Path(args.output_dir or Path(__file__).parent / "results/published")
            output.mkdir(parents=True, exist_ok=True)
            res = published_results()
            res.to_csv(output / "published_scores.tsv", sep="\t", index=False)
            fig = gt.pl.ranking(res, level="task", return_fig=True)
            fig.savefig(output / "published_ranking_task.png", dpi=160, bbox_inches="tight")
            fig.savefig(output / "published_ranking_task.pdf", bbox_inches="tight")
            print(f"Stored published results: {res.name.nunique()} methods, {res.task.nunique()} tasks. Saved to {output.resolve()}")
            return
        grns = dict(item.split("=", 1) for item in args.grn)
        if any(Path(name).name != name or name in {"", ".", ".."} for name in grns):
            raise ValueError("Use simple method names in --grn NAME=PATH.")
        terms = json.loads(Path(args.terms_json).read_text()) if args.terms_json else gt.get_terms(args.context) if args.context else None
        config = GretaConfig(args.n_cells, args.seed, args.edges_per_tf, args.profile)
        result = run_greta(config, args.output_dir, args.data, grns, terms, args.metrics)
        print(gt.tl.stats(result["networks"]).to_string(index=False))
        print(result["coverage"].to_string(index=False))
        print(f"Saved to {result['output_dir'].resolve()} ({result['elapsed_seconds']:.1f} seconds excluding imports)")
    except (ValueError, FileNotFoundError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
