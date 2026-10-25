"""Data-backed illustrations for the RNA-only GRETA practical.

The notebook calls each plotting function beside its task explanation. Running
this module also exports every example and its underlying tables, using the
same benchmark result. These illustrations do not change benchmark scores.
"""

from pathlib import Path
import re

import matplotlib.pyplot as plt

from literature_examples import plot_marker_example, plot_tf_pair_example
from mechanistic_examples import plot_knocktf_example
from predictive_examples import plot_gene_sets_example, plot_omics_example


def save_example(result, name, figure, tables):
    """Save an inline example and preserve the data used to draw it."""
    if not re.fullmatch(r"[a-z][a-z0-9_]*", name):
        raise ValueError("Use a short lowercase name for an example.")
    output = Path(result["output_dir"]) / "examples"
    output.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        figure.savefig(output / f"{name}.{suffix}", dpi=160, bbox_inches="tight")
    for table_name, table in tables.items():
        if not re.fullmatch(r"[a-z][a-z0-9_]*", table_name):
            raise ValueError("Use lowercase names for exported example tables.")
        table.to_csv(output / f"{name}_{table_name}.tsv", sep="\t", index=True)
    return output


def main():
    """Prepare the five examples independently of a notebook session."""
    from rna_benchmark import run_rna

    result = run_rna(verbose=False)
    examples = (
        ("gene_sets", plot_gene_sets_example),
        ("omics", plot_omics_example),
        ("tf_markers", plot_marker_example),
        ("tf_pairs", plot_tf_pair_example),
        ("knocktf", plot_knocktf_example),
    )
    for name, plot in examples:
        figure, tables = plot(result)
        output = save_example(result, name, figure, tables)
        plt.close(figure)
        print(f"Saved {name} figure and data to {output}", flush=True)


if __name__ == "__main__":
    main()
