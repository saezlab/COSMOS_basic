# TF activity estimation in Python

Open **[01_tf_activity_decoupler.ipynb](01_tf_activity_decoupler.ipynb)** for the Python counterpart of the RNA/TF section of [prepare_cosmos_inputs.R](../scripts/prepare_cosmos_inputs.R).

The question is: **which signed TF-target expression footprints support higher or lower activity in each cell line relative to the NCI60 panel?** The notebook uses Python's `decoupler` package and the same ULM method, existing cleaned RNA, and local CollecTRI snapshot as the R script.

## Run in JupyterLab

If you already installed the GRN practical's environment, no additional packages are needed. From the repository root:

```bash
GRN_comparison/.venv/bin/python -m jupyterlab TF_activity/01_tf_activity_decoupler.ipynb
```

Select **Python (GRN comparison)**. On this teaching Mac, retain the GRN practical's OpenMP setting **before launching Jupyter**:

```bash
export DYLD_FALLBACK_LIBRARY_PATH="$(R RHOME)/lib:${DYLD_FALLBACK_LIBRARY_PATH:-}"
GRN_comparison/.venv/bin/python -m jupyterlab TF_activity/01_tf_activity_decoupler.ipynb
```

ULM does not use XGBoost, but importing decoupler also checks installed optional packages. Reusing an environment that contains XGBoost therefore requires its OpenMP library to load correctly. The separate environment below does not install XGBoost and needs neither R nor this setting. No GPU or reference download is needed.

For a separate environment, use Python 3.12. From the repository root on macOS/Linux:

```bash
cd TF_activity
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt jupyterlab
python -m ipykernel install --sys-prefix --name cosmos-tf --display-name "Python (COSMOS TF activity)"
python -m jupyterlab 01_tf_activity_decoupler.ipynb
```

Choose **Python (COSMOS TF activity)** for that environment. On Windows, use `py -3.12 -m venv .venv` and `.\.venv\Scripts\Activate.ps1`.

## What is reproduced?

| R preprocessing / estimation | Python equivalent |
| --- | --- |
| Read `data/RNA/RNA_log2_FPKM_cleaned_common.csv` | Read the same gene × cell-line matrix |
| Keep genes missing in at most 20 cell lines | The same missing-value filter |
| Subtract gene mean and divide by sample SD across cell lines | Gene-wise Z-scores with `ddof=1`, omitting missing values |
| Drop missing values separately in each cell line | Separate finite-gene input for each ULM call; no imputation |
| `regulons` from `support/collectri_regulon_R.RData` | Its bundled text export, `GRN_comparison/data/collectri.tsv`; `score` becomes `weight` |
| `run_ulm(..., minsize=20)` | `dc.mt.ulm(..., tmin=20, empty=False)` per cell line |

The cleaned input already has low-expression measurements censored to missing values by `scripts/RNA.R`. This tutorial preserves that existing preprocessing and the **56-cell-line** panel. It uses a different input from the GRN comparison practical's uncensored 60-line benchmark.

Setting `empty=False` preserves observed Z-scores equal to zero; these are valid measurements, not missing values. A TF needs at least 20 observed targets in the particular cell line. A missing activity means insufficient coverage, not activity zero. TFs are scored through their targets and do not need their own transcript measured.

The ULM score is the **t-statistic of the fitted relationship** between signed regulatory weights and gene Z-scores. It is neither a log fold change nor a direct measurement of TF protein. Scores depend on the cell-line panel used as the expression reference. The downstream R script's `abs(score) > 2` is shown as a heuristic seed filter; it is not an FDR cutoff. All scores are exported before filtering.

## Outputs and COSMOS handoff

Figures and tables go to `TF_activity/results/`. The notebook exports activities, raw p-values, within-cell BH-adjusted p-values, measured-target counts, scaled RNA, and provenance. Its helper is [scripts/tf_activity_decoupler.py](../scripts/tf_activity_decoupler.py).

This notebook covers RNA preprocessing and TF activity estimation. It does not regenerate metabolomic inputs or overwrite `data/cosmos/cosmos_inputs.RData`. The notebook includes an optional R snippet for reading the exported matrices into the existing COSMOS workflow.

## Numerical validation

The notebook compares its computed activities with a bundled text export of the archived R `TF_scores`, used only after inference. All **10,548 tested TF–cell pairs** match, with maximum score difference below **10⁻¹²** in the tested environment. The score matrix contains **202 TFs × 56 cell lines**; individual cell lines have 149–200 tested TFs because measured-target coverage varies.

Three scientific checks cover sample-SD scaling and missingness, zero-valued measurements and the exact 20-target boundary, the full ULM background, and ambiguous duplicate edges. They can be run from the repository root with `python -m pytest scripts/test_tf_activity_decoupler.py`; installing `pytest` is optional for users who only run the notebook.

The [reference export script](export_r_reference.R) uses base R only and is needed solely to regenerate the bundled validation data. Original R package versions were not recorded in the archive; the export provenance states this explicitly.

## Method references

- [Python decoupler ULM](https://decoupler.scverse.org/en/latest/api/generated/decoupler.mt.ulm.html)
- [R decoupleR ULM](https://saezlab.github.io/decoupleR/reference/run_ulm.html)
- [Local CollecTRI export provenance](../GRN_comparison/data/provenance.json)
