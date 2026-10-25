# Comparing data-inferred and prior-knowledge GRNs

This practical compares small gene regulatory networks (GRNs) built from the same NCI60 RNA expression matrix used elsewhere in this repository. Students construct correlation-based networks, compare them with a literature-derived CollecTRI prior and a random baseline, and test how well each network's chosen regulators predict expression in held-out cell lines.

The session is designed for **60–90 minutes**. It uses local inputs and a small, fixed network budget; running the analysis does not download datasets or evaluation databases. Python package installation may require internet access.

The exercise uses [gretapy](https://github.com/saezlab/gretapy) for network statistics and overlap. Its prediction exercise is a transparent, separate Ridge regression analysis. This is a reduced RNA-only teaching benchmark: it does not run `gt.tl.benchmark` or reproduce the full RNA+ATAC GRETA benchmark. There are no measured chromatin accessibility data here, so no artificial CRE or ATAC features are created.

## Learning objectives

By the end of the session, students should be able to:

- Distinguish a network inferred from expression associations from a curated prior-knowledge network.
- Compare networks on a common TF/gene universe while controlling the number of selected edges.
- Explain why network overlap, reference agreement and held-out prediction answer different questions.
- Recognize circular evaluation, information leakage, missing-data assumptions and limits on causal interpretation.
- Explain how the choice of a GRN prior can matter for the TF activity estimates used in the COSMOS practical.

## Setup and first run

Use Python **3.12**, the tested and locked teaching environment. Gretapy itself supports Python >=3.11, but the lock file here targets 3.12. Work in this subfolder so its environment and outputs are easy to find:

```bash
cd GRN_comparison
uv venv --python 3.12 .venv
source .venv/bin/activate
uv pip install -r requirements-lock.txt
python run_benchmark.py
```

If `uv` is unavailable, create an environment with Python and install the same locked requirements:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-lock.txt
python run_benchmark.py
```

On Windows PowerShell, activate the environment with `.\.venv\Scripts\Activate.ps1` instead of the `source` command.

Open `01_grn_comparison.ipynb` in VS Code or an existing Jupyter environment and select the `.venv` interpreter/kernel for the guided version. The locked requirements include kernel support but do not install a Jupyter server. To start a local server if needed, install JupyterLab with `python -m pip install jupyterlab`, then run `jupyter lab 01_grn_comparison.ipynb`. The notebook and command-line runner use the same analysis helpers in `grn_benchmark.py`.

The runner also works from the repository root:

```bash
python GRN_comparison/run_benchmark.py
```

Use `python run_benchmark.py --help` for the available settings. `requirements.txt` describes the direct dependencies; `requirements-lock.txt` records the environment to use for the practical.

The package environment is larger than the analysis: gretapy imports its multiome and XGBoost dependencies even though this exercise uses only network statistics. Install it before class. Gretapy 0.0.3 also requires the `mudata<0.4` constraint included here; an unconstrained installation currently fails on import.

On macOS, if importing gretapy reports that `libomp.dylib` cannot be loaded, install OpenMP with `brew install libomp`. On this repository's teaching Mac, the existing R installation already provides it; the following shell setting also allowed the workflow to run without installing another copy:

```bash
export DYLD_FALLBACK_LIBRARY_PATH="$(R RHOME)/lib:${DYLD_FALLBACK_LIBRARY_PATH:-}"
```

Set this before starting Python or Jupyter. On Windows, XGBoost may require the Visual C++ Redistributable; see its [installation guide](https://xgboost.readthedocs.io/en/stable/install.html). Python 3.12 on macOS ARM64 is locally tested; the lock includes platform markers, but other operating systems still need a pre-course installation check.

## Inputs and provenance

The expression input is `../data/RNA/RNA_log2_FPKM_cleaned_common.csv`: genes by cell lines, containing the existing log2(FPKM + 1) measurements. The analysis uses expression values directly, rather than the per-cell-line z-scores or CollecTRI-derived TF activity scores prepared for COSMOS.

The local `data/collectri.tsv` contains 43,178 rows exported from `../support/collectri_regulon_R.RData`. `export_collectri.R` regenerates that export and its `data/provenance.json` checksums using base R. The original snapshot's release was not recorded; it is not assumed to be current. The benchmark does not fetch a newer prior during a run. CollecTRI is a literature-curated resource, not a complete or NCI60-specific catalogue of true interactions; its absence of an edge does not establish that the interaction is false.

Here, “TFs” is shorthand for CollecTRI source transcriptional regulators; some entries are cofactors rather than sequence-specific DNA-binding TFs.

There is an upstream preprocessing limitation. The repository's `scripts/RNA.R` replaced low expression values with `NA` and filtered genes across the full panel before this practical's split. Missingness therefore partly represents censored low expression. Selection and fitting **within this module** use training cell lines only, but this is not a pristine holdout starting from untouched raw data. Prediction scores describe the observed expression entries retained in this input.

## Worked session

### 1. Inspect the data and define a fair comparison — 10–15 minutes

Inspect the matrix dimensions, missing values and cell-line names. Rows represent genes and columns represent different cancer cell lines; these are not repeated measurements from one cell type.

The default run uses a random 70:30 training/test split with `--split-seed 42`. Genes must be observed in at least 80% of training cell lines and have nonzero training variance. It selects up to 500 variable target genes using training observations and 20 TFs using training expression variance. Candidate TFs must also be CollecTRI sources and have enough eligible prior edges to fill the selected edge budget. All networks use this same eligible TF/target universe and exclude self-edges.

This is a deliberately **prior-guided TF universe**. It makes the comparison manageable and ensures the curated network has comparable coverage, while limiting conclusions to those TFs and genes. It does not test inference for arbitrary genes or previously unannotated TFs.

Record the split, selected genes and parameters before inspecting prediction scores. Do not choose parameters because they produce the best test result.

### 2. Construct four small networks — 15–20 minutes

The default budget is 10 outgoing edges per TF, giving 200 edges in each selected network when all eligibility requirements are met:

| Network | How its TF-to-target edges are selected |
| --- | --- |
| Pearson | Largest absolute TF–target Pearson correlations in training expression |
| Spearman | Largest absolute TF–target Spearman correlations in training expression |
| CollecTRI | A reproducible random sample of eligible curated edges for each TF |
| Random | Uniform random eligible targets for each TF, using the same outgoing-edge budget |

The prior is subsampled because its signed weights do not provide a comparable continuous confidence ranking. The full eligible prior is retained as a descriptive reference, not added as a fifth competitor with a larger edge budget. Network sampling uses `--network-seed 42`, independently of the sample split seed.

Correlation estimation uses training-median-filled expression for both TFs and targets. This is a simple teaching approximation that can alter correlations, especially when low expression is censored. Prediction also fills missing TF predictors with training medians, while preserving missing target responses for exclusion from fitting/scoring. Correlation signs describe positive or negative associations in this dataset. CollecTRI signs describe curated regulatory effects. These meanings differ, and a negative correlation is not evidence of direct repression. Direction in the expression-derived networks comes from designating the source as a TF; correlation itself is symmetric.

The random baseline matches the selected TFs and their outgoing degree. It does not preserve each target's incoming degree, so target coverage and regulator allocation may differ between networks.

### 3. Compare network structure and prior agreement — 10–15 minutes

Use `gt.tl.stats` to compare the numbers of TFs, target genes and unique TF–target edges, together with mean regulon size. Equal edge counts need not give equal target coverage.

Use `gt.tl.ocoeff` to compare overlap. For edge sets A and B, its overlap coefficient is:

```text
overlap(A, B) = number of shared edges / min(number of edges in A, number in B)
```

The edge statistic counts unique TF–target pairs and ignores edge signs. An overlap of 1 means the smaller network is contained in the larger one; it establishes agreement, not regulatory correctness. With equal-sized networks, it is the fraction of selected edges shared by the two methods. CRE comparisons are unavailable for these RNA-only networks.

The separate prior-agreement table compares selected edges with the full eligible CollecTRI prior. **CollecTRI's agreement with its own source is circular by construction.** Treat it as a check on selection and coverage, not independent evidence that CollecTRI outperforms the inferred networks. Agreement of other networks with this reference is also not a complete measure of biological accuracy.

### 4. Test held-out expression prediction — 15–20 minutes

For each selected target gene, use the network's incoming TFs as predictors in a Ridge regression. Each method learns regression coefficients using training cell lines only, with standardized predictors and the same fixed regularization strength (`--alpha 10`). Edge scores and curated signs select or describe the support; they are not imposed as the fitted regression coefficients.

Fit targets only on observed training responses and evaluate them only on originally observed test responses. A minimum of five observed responses is required in each partition. Missing test responses are not replaced by medians and scored as if measured. All methods are assessed on the same evaluable target genes, including targets with no selected regulators: those receive the observed training-target mean as their prediction.

Compare mean squared error (MSE), improvement relative to the training-mean baseline, target coverage and per-target test R². Inspect both the summary and per-target tables:

- Negative test R² is meaningful and should be retained. Standard R² compares error with variation around the **test** mean, whereas the practical's explicit baseline predicts the **training** mean.
- Constant test responses have undefined R²; use the reported eligibility and observation counts to interpret omissions.
- Good prediction can reflect shared tissue programs or indirect associations. It does not establish a direct TF-to-gene regulatory mechanism.

The holdout is protected from this module's correlation selection and coefficient fitting. Nevertheless, a single small random split estimates prediction within this panel, and does not establish performance in a new tissue, experiment or perturbation.

### 5. Interpret the comparison — 10–15 minutes

Write a short interpretation that addresses all three views: structure, prior agreement and prediction. Give the observed values from your own run, and identify where the evidence is insufficient to choose a method.

Useful discussion prompts:

- Do Pearson and Spearman agree on edges? Does similar overlap imply similar prediction?
- Which methods spread edges over more targets? How does the training-mean fallback affect their overall prediction scores?
- Can a network predict well while agreeing poorly with CollecTRI? What biological and statistical explanations are plausible?
- Which comparison is circular, and what independent evidence would strengthen the conclusion?
- Would the preferred network necessarily be the best prior for TF activity estimation or COSMOS? What task-specific evidence is still missing?

## Outputs

The selected output directory contains the following records. Inspect `run_info.json` before comparing runs so differences in the universe, split or parameters are not mistaken for method effects.

| File | Purpose |
| --- | --- |
| `networks/*.tsv` | Four selected edge lists |
| `collectri_reference.tsv` | Full eligible prior, retained as a descriptive reference |
| `network_stats.tsv` | gretapy network size and coverage statistics |
| `overlap.tsv` | Pairwise overlap of sources, targets and TF–target edges |
| `prior_agreement.tsv` | Descriptive agreement with the full eligible prior |
| `prediction_per_target.tsv` | Held-out results and observation counts for individual targets |
| `prediction_summary.tsv` | Prediction and coverage summaries across common evaluable targets |
| `split.tsv` | Training/test cell-line assignments |
| `genes.tsv` | Selected gene universe |
| `run_info.json` | Parameters, versions and run metadata |
| `comparison.png`, `comparison.pdf` | Comparison figure for discussion or export |

## Exercises

Keep outputs from each run in a different directory. The examples below run from `GRN_comparison/`.

1. **Change network density.** Compare outgoing-edge budgets of 5 and 10 while holding the seeds fixed:

   ```bash
   python run_benchmark.py --edges-per-tf 5 --output-dir results/k5
   python run_benchmark.py --edges-per-tf 10 --output-dir results/k10
   ```

   Check whether the eligibility rule changed the TF or target universe before interpreting the difference. A larger budget can reduce the number of TFs with enough eligible prior interactions.

2. **Change the data split.** Repeat with several `--split-seed` values while keeping the network seed fixed:

   ```bash
   python run_benchmark.py --split-seed 43 --network-seed 42 --output-dir results/split43
   ```

   Compare score variability and rank changes. This changes training observations, so variable-gene selection and correlations can change too. Repeated splits of the same panel are sensitivity analyses, not independent biological replications.

3. **Change prior/random edge sampling.** Keep the split fixed and change only the network seed:

   ```bash
   python run_benchmark.py --split-seed 42 --network-seed 43 --output-dir results/network43
   ```

   Compare the stability of the subsampled prior and random baseline. Distinguish sensitivity to edge subsampling from sensitivity to which cell lines are available for training.

4. **Change the modeling assumptions.** Compare a few prespecified values of `--alpha`, or change `--n-tfs` and `--n-targets`. Explain how those choices affect coverage and model complexity. Selecting the best setting using these test scores turns the test set into a validation set; a new untouched test set would then be needed for a final performance claim.

5. **Design a stronger evaluation.** Propose a tissue-aware split, a target-degree-preserving randomization, or evaluation against a genuine independent TF perturbation dataset. Identify what extra data and controls are required. Binding evidence, perturbation responses, expression predictability and literature agreement each support different claims; none alone proves every selected edge is direct and causal.

## Instructor validation

The default command-line workflow and all 22 notebook cells were executed with the locked Python 3.12 environment on macOS ARM64. The default split contains 39 training and 17 test cell lines. It creates four 200-edge networks and scores all 500 candidate targets. The analysis took approximately two seconds after imports on the teaching machine; installation and first-import caches take longer.

The default prediction gains over the training-mean baseline are small. Use this as an interpretation exercise; the single run does not establish a clear biological winner. Three focused safeguards also passed: changing held-out values cannot change selected features or networks, changing the edge-sampling seed leaves correlation networks unchanged, and missing target responses are excluded while targets without parents remain in the evaluation.

To rerun those checks from the repository root, with the environment active and any required OpenMP setting applied:

```bash
python -m unittest discover -s GRN_comparison -v
```

## Relation to GRETA and COSMOS

The [official gretapy comparison tutorial](https://gretapy.readthedocs.io/en/latest/notebooks/usage.html) demonstrates literature, correlation and random baselines with TF–CRE–gene networks. This practical retains the TF–gene comparison and uses gretapy's statistics tools without constructing unsupported chromatin features. Its [statistics implementation](https://github.com/saezlab/gretapy/blob/main/src/gretapy/tl/_stats.py) supports these TF–target edge lists directly.

The [official evaluation tutorial](https://gretapy.readthedocs.io/en/latest/notebooks/eval.html) covers additional genomic, literature, predictive and mechanistic tasks. Those tasks require appropriate reference resources and, for some metrics, paired RNA+ATAC data. A future extension can add a small real multiome dataset or independent perturbation experiment, with explicit checks for overlap between evidence used to build and evaluate the networks.

In the existing COSMOS workflow, CollecTRI targets are used to estimate TF activities from RNA expression. Comparing a GRN against TF activities computed from that same GRN would introduce another form of circularity. The held-out exercise here instead predicts target expression from TF expression. It evaluates that prediction task, and does not directly rank priors for footprint-based activity estimation or causal network reconstruction.
