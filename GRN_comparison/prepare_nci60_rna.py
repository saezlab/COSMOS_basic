"""Regenerate the uncensored teaching input from the repository's CellMiner XLS.

Only regeneration requires xlrd==2.0.2. The tutorial reads the bundled TSV/GZ.
"""

from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent


def main():
    source = HERE.parent / "data/RNA/RNA__RNA_seq_composite_expression.xls"
    raw = pd.read_excel(source, skiprows=10, sheet_name=0, engine="xlrd")
    expression = raw.iloc[:, 6:].apply(pd.to_numeric, errors="raise")
    genes = raw.iloc[:, 0].astype(str)
    if genes.duplicated().any() or not np.isfinite(expression.to_numpy()).all():
        raise ValueError("Expected unique genes and finite expression in the source workbook.")
    labels = [str(name).split(":", 1) for name in expression.columns]
    if any(len(parts) != 2 for parts in labels):
        raise ValueError("Expected CellMiner tissue:cell-line column labels.")
    samples = pd.DataFrame({"sample": [p[1] for p in labels], "tissue_code": [p[0] for p in labels]})
    if samples["sample"].duplicated().any():
        raise ValueError("Cell-line names must be unique.")
    expression.columns = samples["sample"].tolist()
    expression.insert(0, "gene", genes)
    output = HERE / "data/nci60_rna_log2fpkm.tsv.gz"
    expression.to_csv(output, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
    sample_file = HERE / "data/nci60_samples.tsv"
    samples.to_csv(sample_file, sep="\t", index=False)
    provenance = dict(
        source="data/RNA/RNA__RNA_seq_composite_expression.xls",
        source_sha256=sha256(source.read_bytes()).hexdigest(),
        expression_sha256=sha256(output.read_bytes()).hexdigest(),
        sample_metadata_sha256=sha256(sample_file.read_bytes()).hexdigest(),
        source_description="Existing CellMiner NCI60 RNA-seq composite workbook; log2(FPKM + 1), as documented in scripts/RNA.R.",
        preparation="First sheet, skip 10 header rows; retain original gene names and all expression values. No censoring, imputation, gene selection or sample split.",
        tissue_codes="Original column prefixes, retained for stratification; bulk tissue labels are not single-cell identities.",
        n_genes=len(expression), n_samples=len(samples),
        regeneration="python GRN_comparison/prepare_nci60_rna.py (requires xlrd==2.0.2)",
    )
    (HERE / "data/nci60_rna_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Prepared {len(expression):,} genes × {len(samples)} cell lines; uncensored values, no missing entries.")


if __name__ == "__main__":
    main()
