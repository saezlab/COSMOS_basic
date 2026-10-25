> Retained preview for the later RNA + ATAC session. Start with [the RNA-only practical](README.md) for the current course.

# Evaluating data-driven and prior-knowledge GRNs with GRETA

This practical compares small gene regulatory networks (GRNs) using gretapy's **native evaluation tasks** and the task-level ranking figure from the [official evaluation tutorial](https://gretapy.readthedocs.io/en/latest/notebooks/eval.html):

```python
gt.pl.ranking(res, level="task")
```

Start with **[03_rna_atac_preview.ipynb](03_rna_atac_preview.ipynb)** and keep **[RNA_ATAC_TASKS.md](RNA_ATAC_TASKS.md)** open as the student guide. The emphasis is on what each evaluation task measures, what a score supports, and why different tasks can prefer different networks. Network inference is a short introductory step; students can subsequently bring their own GRNs inferred from single-cell RNA + ATAC.

The main session uses gretapy's **synthetic, PBMC-like paired RNA + ATAC data**. This makes genuine multimodal evaluation practical on a laptop. It is a teaching example, not a biological comparison of inference methods. The original NCI60 exercise remains available in [NCI60_EXTENSION.md](NCI60_EXTENSION.md) and [01_grn_comparison.ipynb](01_grn_comparison.ipynb), with its original runner, `run_benchmark.py`.

## What the session does

The practical has two separate result panels:

- **Live reduced benchmark:** generate small Pearson, Spearman, CollecTRI, and random networks, run selected native GRETA tasks, and draw `gt.pl.ranking(..., level="task")`.
- **Full reference panel:** load the authors' stored manuscript scores for 20 methods and all 11 task types, then draw the full task-level figure. This uses existing scores rather than recalculating the full benchmark or adding an imaginary method.

Keep the panels separate: toy results and manuscript results have different datasets, network sizes, and evaluation coverage. They cannot support a combined method ranking.

By the end, students should be able to explain the four evaluation classes, distinguish evidence about binding, prediction, prior agreement, and perturbation response, read F₀.₁ colors and within-column ranks, and identify circularity and missing evaluation coverage.

## Environment and first run

Reuse the existing Python **3.12** environment and `requirements-lock.txt`; no extra scientific dependencies are needed for this update. From this folder, a fresh persistent setup is:

```bash
cd /Users/dugourd/Documents/COSMOS_basic/GRN_comparison
uv venv --python 3.12 .venv
source .venv/bin/activate
uv pip install -r requirements-lock.txt
```

If `uv` is unavailable, use `python3.12 -m venv .venv`, activate it, and run `python -m pip install -r requirements-lock.txt`. On Windows PowerShell, activate with `.\.venv\Scripts\Activate.ps1`. Use the locked environment: it pins gretapy 0.0.3 and a compatible MuData version.

On this teaching Mac, R supplies the OpenMP library required by XGBoost. Set this **before starting Python or Jupyter**:

```bash
export DYLD_FALLBACK_LIBRARY_PATH="$(R RHOME)/lib:${DYLD_FALLBACK_LIBRARY_PATH:-}"
```

On another Mac without that library, `brew install libomp` is an alternative. Other operating systems need a pre-course installation check; see the [XGBoost installation guide](https://xgboost.readthedocs.io/en/stable/install.html) for platform prerequisites.

**JupyterLab is recommended** for teaching: students can run cells, inspect tables, and discuss the figure as it develops. The locked requirements include a Python kernel; install the notebook interface separately if needed:

```bash
uv pip install --python .venv/bin/python jupyterlab
python -m jupyterlab 03_rna_atac_preview.ipynb
```

Choose the `.venv` kernel and run cells from the top with **Shift+Enter**. Keep the Terminal session running while JupyterLab is open. VS Code with the same Python interpreter is also suitable.

To run the same live analysis from Terminal:

```bash
python run_greta_benchmark.py
```

The default writes to `results/greta/`. For a smaller first run that omits the external KnockTF perturbation resource:

```bash
python run_greta_benchmark.py --profile quick --output-dir results/quick
```

To draw only the full manuscript-style reference panel:

```bash
python run_greta_benchmark.py --published-only
```

This last command uses the bundled score snapshot and does not need the live evaluation databases.

## Downloads and preparation before class

The default teaching profile downloads approximately **84 MB** of reference files on a clean run. Most is the approximately 74 MB compressed KnockTF resource. The final reference cache occupies roughly **170 MB**, because KnockTF is stored decompressed. The quick profile omits KnockTF and downloads approximately 10 MB. Package installation is separate and larger than these reference downloads.

References are cached in `GRN_comparison/gretapy_data/`, regardless of the working directory. The cache, environments, and results are ignored by Git. Run the intended profile before class to populate the cache, test the installation, and avoid relying on classroom connectivity. No command here automatically starts the full multi-dataset benchmark described upstream.

The small bundled `data/greta_published_metrics.csv.gz` snapshot is about **152 kB**. It is copied into the cache and read through `gt.ds.read_metrics()`. Attribution, source record, license (**CC BY 4.0**), resource sizes, and checksums are recorded in [data/greta_resources.json](data/greta_resources.json). The source is the [GRETA Zenodo record](https://zenodo.org/records/19600656); see the [GRETA manuscript](https://doi.org/10.1101/2024.12.20.629764).

## Classroom flow: about 60–90 minutes

### 1. Inspect the example and understand network construction — 10–15 minutes

The default generates 180 synthetic cells using `gt.ds.toy(n_cells=180)`, with B-cell, T-cell, and monocyte labels. Its hardcoded example GRN is discarded. We infer fresh correlation networks using the public `gt.mt.correlation` helper and obtain the prior baseline using `gt.mt.lit_grn`.

At a high level, correlation inference connects candidate TFs to genes with strong expression associations and attaches accessible promoter-overlapping regions. CollecTRI supplies curated TF–target links, which are filtered against measured genes and promoter-overlapping peaks. A common accessible promoter is retained per target. Random edges provide a simple baseline.

All four networks share eligible source TFs and the same outgoing-edge budget for each TF. By default this is capped at five distinct targets per TF, with at least three eligible prior targets required. The tested default gives **7 TFs and 31 TF–target edges per network**. This is a prior-guided, promoter-constrained comparison; it does not test arbitrary regulators or enhancer discovery.

Correlation estimation uses **120 training cells**, leaving **60 test cells** untouched for edge selection. This matches the pinned native Omics evaluator's split exactly: `test_size=0.33`, `random_state=42`, stratified by `celltype`. The configurable `--seed` controls toy generation and prior/random sampling; it does not change that native evaluation split. A TF-to-target direction is imposed by designating TF sources; correlation alone does not establish regulatory direction or signs.

### 2. Understand the tasks before calculating scores — 15–20 minutes

Use [RNA_ATAC_TASKS.md](RNA_ATAC_TASKS.md) for explanations of all eleven task types. The default **teaching** profile requests eleven database/subtask columns spanning six task types and all four classes:

| Class | Live task | Databases or prediction directions | Main question |
|---|---|---|---|
| Genomic | CREs | ENCODE CREs; ENCODE Blacklist | Do selected regions overlap regulatory annotations and avoid problematic regions? |
| Literature | TF Markers | HPA; TF-Marker | Are the selected source TFs supported as context-relevant markers? |
| Literature | TF Pairs | IntAct; Europe PMC | Do TFs sharing targets agree with reported TF associations? |
| Predictive | Gene Sets | Hallmarks | Does network functional coverage agree with pathways active in expression? |
| Predictive | Omics | `Gene ~ TFs`; `Gene ~ CREs`; `CRE ~ TFs` | Can the network's chosen predictors explain held-out molecular variation? |
| Mechanistic | TF Scoring | KnockTF (scoring) | Does a signed regulon identify decreased TF activity in external knockdown profiles? |

The **quick** profile uses the same setup without KnockTF, covering three classes. Tissue-specific references use the PBMC context terms for this toy example.

Two additional diagnostic columns, **Promoters** and **CollecTRI reference agreement**, are calculated but excluded from the comparison ranking. Promoter-based construction makes promoter agreement partly predetermined. A CollecTRI-derived network's agreement with CollecTRI is self-reference. Keep these checks visible without treating them as independent validation.

The full reference panel additionally illustrates TF Binding, CRE to Gene, Perturbation Forecasting, and Steady State Simulation. These tasks have different evidence and coverage requirements; they are explained in the task guide rather than all being recomputed in the live subset.

### 3. Run the benchmark and inspect eligibility — 15–20 minutes

The notebook calls the native `gt.tl.benchmark` on the selected networks and references. Inspect `coverage.tsv` alongside precision, recall, and F₀.₁. For Omics, the helper records the number of tested features, significant hits, and constant-response features excluded because their correlation is undefined. For TF Scoring, it records eligible and testable perturbation experiments.

**Zero and unavailable are different outcomes.** No eligible tests are recorded as missing values. A task with eligible tests and no successful hits can legitimately score zero. On the tested default, three KnockTF profiles are eligible and two have usable regulons; all four methods score zero in this column. This limited evidence does not establish that their biological regulation is absent.

Before ranking, the helper excludes any database/subtask column missing a score for any compared method. This keeps the ranked methods on a common set of available columns. Inspect `excluded_from_ranking.tsv` to understand what was removed. Equal availability does not imply equal numbers of tested features, so retain the coverage discussion.

### 4. Read the task-level figure — 15–20 minutes

The ranked table is passed directly to:

```python
gt.pl.ranking(res, level="task")
```

Color shows mean **F₀.₁**; the printed number shows the rank within that database/subtask column, with 1 best and average ranks for ties. The score strongly emphasizes precision. The row ordering averages databases within tasks, tasks within classes, combines the available classes with their weights per dataset, then averages datasets. The guide explains the formula and how missingness affects aggregation.

Explain differences across columns before naming a preferred network. In particular:

- Identical TF source sets make TF Marker scores tie even when their edges differ.
- A common single promoter per target means `Gene ~ CREs` largely compares which genes are included and their accessibility–expression association, rather than TF choice.
- Correlation-selected predictors can explain cell-type differences without demonstrating direct regulation.
- The synthetic biological structure and tiny networks limit generalization. Low scores and ties are useful discussion material.

Change class weights as an interpretation exercise. The underlying task scores stay the same; the scientific preference used to order methods changes.

### 5. Read the full reference panel and design an evaluation — 10–15 minutes

Load the authors' stored scores and draw their full task-level panel. Identify which extra columns require binding collections, region–gene evidence, perturbation models, or signed dynamical assumptions. These are results from the original benchmark, not additional observations from our toy run.

Finish by choosing tasks for a student's own GRN. Ask which evidence was already used during inference, what independent evidence remains, and whether the network has the necessary regions, signs, cell types, and reference coverage. [RNA_ATAC_TASKS.md](RNA_ATAC_TASKS.md) includes six discussion exercises.

## Bringing an existing RNA + ATAC GRN

The current runner supports a prepared **human hg38** MuData file and one or more existing GRN tables. It also constructs the four teaching baselines on the supplied data for comparison.

The MuData must contain:

- `rna` and `atac` modalities with the same ordered, unique cell identifiers.
- A nonmissing `mdata.obs["celltype"]` annotation, with enough observations per class for the stratified split.
- Finite, nonnegative normalized RNA and ATAC values, not centered/scaled matrices.
- Unique RNA gene symbols compatible with the human reference resources and ATAC regions matching the network. Regions use `chrX-start-end` coordinates in hg38; map Ensembl identifiers to compatible symbols before this exercise.

Each GRN needs **`source`, `cre`, `target`, `score`** columns. Sources and targets must occur in RNA feature names; CREs must occur in ATAC feature names. Scores must be finite; biological signs are needed for signed mechanistic interpretations. Duplicate triples are averaged, self-edges removed, and at least five remaining edges are required.

The pinned helpers use dense matrices. This teaching wrapper converts sparse inputs only when each modality has at most **50 million entries**; first subset larger data to a suitable teaching size. Custom data must also support the promoter-linked CollecTRI baseline. This is a bounded classroom workflow, not a large-scale production benchmark.

Choose the reference-filtering context explicitly:

```bash
python run_greta_benchmark.py \
  --data my_data.h5mu \
  --grn MyGRN=my_network.tsv \
  --context PBMC \
  --output-dir results/my_grn
```

Use `PBMC` only for appropriate blood-cell data. For another context, choose a suitable GRETA context or supply `--terms-json my_terms.json`: a JSON mapping database names to relevant reference terms. Inspect `gt.show_terms(organism="hg38")` when selecting terms. Passing `--grn` repeatedly adds multiple student methods.

The same entry point is available interactively:

```python
import gretapy as gt
from greta_benchmark import configure_cache, run_greta

configure_cache()
result = run_greta(
    dataset_path="my_data.h5mu",
    grn_paths={"MyGRN": "my_network.tsv"},
    terms=gt.get_terms("PBMC"),  # Choose the context for your data.
    output_dir="results/my_grn",
)
```

The wrapper cannot undo leakage in an existing GRN inferred from all supplied cells. Record its inference data and preprocessing before interpreting native held-out Omics scores. For transfer claims, evaluate on genuinely independent cells, donors, or samples. Also compare network sizes and universes: student networks are validated but are not automatically subsampled to the teaching baselines' edge budget.

To override the profile with selected database/subtask names:

```bash
python run_greta_benchmark.py \
  --metrics "ENCODE CREs" "Hallmarks" "Gene ~ TFs" \
  --output-dir results/selected_tasks
```

The two construction/reference diagnostics remain separately calculated and excluded from ranking. Requesting additional resources can substantially increase downloads or computation. Check `python run_greta_benchmark.py --help` for the supported options.

## Outputs and implementation notes

Live runs save the following under the selected output directory:

| Output | Purpose |
|---|---|
| `networks/*.tsv` | Evaluated TF–CRE–target edge tables |
| `inference_split.tsv`, `candidate_promoters.tsv` | Training/test allocation and the common target–promoter map |
| `network_stats.tsv` | Network size and coverage |
| `evaluation_all.tsv` | All requested scores, including diagnostics and unsupported outcomes |
| `coverage.tsv` | Tested-feature/hit counts for Omics and eligibility counts for TF Scoring |
| `evaluation_ranked.tsv`, `excluded_from_ranking.tsv` | Common comparable columns and explicit exclusion reasons |
| `ranking.tsv` | Native composite ranking values |
| `ranking_task.png/.pdf`, `ranking_class.png/.pdf` | Native GRETA figures |
| `run_info.json` | Configuration, package versions, selected columns, and evaluation notes |

The published-only run writes `published_scores.tsv` and `published_ranking_task.png/.pdf` to `results/published/` by default.

[greta_benchmark.py](greta_benchmark.py) contains the teaching helpers; [run_greta_benchmark.py](run_greta_benchmark.py) provides the CLI. A version-checked local compatibility context repairs three label-routing inconsistencies in gretapy **0.0.3** so requested tasks reach their native implementations. It also excludes constant held-out responses before native prediction, because their undefined Spearman correlations otherwise cause an error, and records testing coverage. Native model algorithms, score formulas and plotting functions are retained, and the original functions are restored on exit. The benchmark's `grn` method column is renamed to the `name` column expected by the ranking API. Keep the pinned version when reproducing this session.

The default live CLI was executed on the teaching Mac: four networks, 52 evaluation rows including two diagnostic columns per method, and 44 ranked rows covering eleven database/subtask columns, six task types, and all four classes. The analysis took approximately **5 seconds after imports** with the reference cache populated. Installation, imports, and initial downloads take longer.

All **28 cells** in the new notebook were executed successfully, including both native ranking figures. A custom MuData/GRN import was also exercised. Eight checks passed across the original and revised workflows, including holdout isolation, requested-task completeness, missing-coverage exclusions, constant-response handling, and matched feature validation. With the environment active, run them from the repository root using:

```bash
python -m unittest discover -s GRN_comparison -v
```

The cached integration checks skip themselves if the teaching reference files have not yet been downloaded. Exact reference-filter terms are saved in `run_info.json` for live runs.

For continuity with the COSMOS/NCI60 material, use the preserved [NCI60 guide](NCI60_EXTENSION.md) and [original notebook](01_grn_comparison.ipynb). That optional extension uses local NCI60 RNA, gretapy topology tools, and a separate held-out Ridge exercise; its scores should be interpreted within that RNA-only exercise.
