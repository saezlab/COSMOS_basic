# RNA-only GRN comparison with GRETA

Compare a large expression-derived gene regulatory network (GRN) with **CollecTRI**, using the existing **NCI60 RNA-seq data** and gretapy's native evaluation tasks. Start with **[02_greta_evaluation.ipynb](02_greta_evaluation.ipynb)**; use **[TASKS.md](TASKS.md)** to explain the results.

The notebook shows the overall heatmap first:

```python
gt.pl.ranking(res, level="task")
```

It then briefly explains how the example networks were constructed and examines **Predictive**, **Literature**, and **Mechanistic** evaluation separately. The aim is to interpret a benchmark: inference provides networks to compare, rather than becoming a second practical. This is an RNA-only application of the [GRETA evaluation framework](https://gretapy.readthedocs.io/en/latest/notebooks/eval.html); region-based tasks need ATAC or other genomic evidence and are outside this session.

## Environment and first run

Reuse your Python **3.12** teaching environment and the existing locked dependencies. To create a persistent environment, run:

```bash
cd /Users/dugourd/Documents/COSMOS_basic/GRN_comparison
uv venv --python 3.12 .venv
source .venv/bin/activate
uv pip install --python .venv/bin/python -r requirements-lock.txt jupyterlab
python -m ipykernel install --sys-prefix --name cosmos-grn --display-name "Python (GRN comparison)"
```

If the environment already exists, activate it; run the kernel-registration command once if the named kernel is not listed. Without `uv`, use `python3.12 -m venv .venv` and `python -m pip install` after activation. On Windows PowerShell, activate with `.\.venv\Scripts\Activate.ps1`.

On this teaching Mac, make R's OpenMP library available **before launching Python or Jupyter**:

```bash
export DYLD_FALLBACK_LIBRARY_PATH="$(R RHOME)/lib:${DYLD_FALLBACK_LIBRARY_PATH:-}"
python -m jupyterlab 02_greta_evaluation.ipynb
```

Choose **Python (GRN comparison)** and run cells from the top with **Shift+Enter**. The named kernel is installed inside `.venv`; it does not change your global Python setup. JupyterLab is recommended for the interactive class. A Mac without R's library can use `brew install libomp`; consult the [XGBoost installation guide](https://xgboost.readthedocs.io/en/stable/install.html) for platform prerequisites.

Prepare the full benchmark before class:

```bash
python run_rna_benchmark.py
```

Results go to `results/rna/`. The tested full run took **62 seconds after imports**, with loading and inference together under one second. Native held-out evaluation fits thousands of small prediction models; other machines or inputs may take longer. Verified cached reuse took **0.9 seconds**. Reuse requires matching input/reference hashes, settings, context terms, code, package versions and saved output hashes. Use `--force` to recompute. A shorter profile omits the KnockTF resource and mechanistic column:

```bash
python run_rna_benchmark.py --profile quick --output-dir results/rna_quick
```

Reference files are downloaded into `gretapy_data/` on their first use. KnockTF is the largest individual download, approximately 74 MB compressed. Prepopulate the cache and run the intended profile before teaching. The bundled expression snapshot needs no extra dependency; regenerating it from the original XLS requires `xlrd==2.0.2`.

### Recovering from an incomplete temporary environment

An error such as `Failed to read metadata ... /private/tmp/cosmos-grn-venv/.../METADATA: No such file or directory` means installed package metadata is missing. That temporary environment may also be missing other files, so reinstalling JupyterLab alone can be insufficient. Use the persistent `.venv` in this folder with the setup commands above. Temporary directories are unsuitable for an environment you want to keep between teaching sessions.

Once `.venv` is installed, launch it explicitly from a fresh Terminal session:

```bash
cd /Users/dugourd/Documents/COSMOS_basic/GRN_comparison
export DYLD_FALLBACK_LIBRARY_PATH="$(R RHOME)/lib:${DYLD_FALLBACK_LIBRARY_PATH:-}"
.venv/bin/python -m jupyterlab 02_greta_evaluation.ipynb
```

Open the URL printed by this new server and select **Python (GRN comparison)**. An older browser tab can still be connected to the temporary environment's server or kernel. Inside the notebook, `import sys; print(sys.executable)` should show this folder's `.venv/bin/python`.

## Real RNA data and larger networks

The input is an uncensored local export of CellMiner's **23,808 genes × 60 cancer cell lines**, in `log2(FPKM + 1)` units. The source workbook, checksums, and extraction procedure are recorded in [data/nci60_rna_provenance.json](data/nci60_rna_provenance.json). This export preserves genuine zero measurements; it does not reuse the censored/imputed matrix from the earlier NCI60 exercise.

The code first reserves **20 test cell lines**, stratified by their original tissue codes, using the pinned evaluator's `test_size=0.33`, `random_state=42` split. It retains genes expressed at least at level 1 in at least 80% of the **40 training lines**, and removes genes constant in training. There is **no default gene or TF cap**. This retains **7,530 genes and 195 source TFs** with at least three eligible prior targets. The full eligible CollecTRI network has **6,871 edges** covering **2,002 distinct targets**; Pearson and Random each have the same 6,871 edges, with potentially different target coverage.

NCI60 consists of bulk cell-line measurements. Tissue labels are saved in `obs["celltype"]` to satisfy GRETA's interface; they are original tissue strata, not single-cell identities. With only 60 lines and several tumour origins, this is a practical comparison on a heterogeneous panel, not a definitive assessment of general GRN accuracy.

Three methods share the same TF identities and outgoing edge count for each TF:

- **Pearson:** correlate every eligible TF with every eligible gene across training lines; retain each TF's strongest absolute correlations up to its CollecTRI degree.
- **CollecTRI:** retain all eligible prior links and their curated signs, without randomly reducing the prior to a tiny network.
- **Random:** draw random target links with the same per-TF edge counts, as a simple control.

The Random control assigns random positive/negative signs; it matches outgoing degrees, not sign distributions. Here “TFs” includes CollecTRI transcriptional regulators and cofactors, rather than only sequence-specific DNA-binding proteins.

The candidate TF list and edge budgets therefore use prior knowledge; Pearson's target choices use expression. This is a **prior-guided comparison**, not fully prior-free discovery. The TF→gene direction comes from specifying TF sources. Correlation signs describe association and are only proxies for activating/repressing effects. Add a rank-correlation baseline with `--include-spearman`.

## Teaching flow: overview first, then explain each category

Allow about **60–90 minutes**, with computation completed beforehand or reused from a valid cache.

1. **Overall heatmap — 10 minutes.** Read database labels, class/task headers, colors, within-column rank numbers, and the composite row order. Identify disagreements before discussing a preferred network. Missing scores are unavailable evidence, not zeros.
2. **Where the networks came from — 10 minutes.** Use the short correlation explanation above. Inspect network sizes, the training/test allocation, and the common regulator/degree budgets. Students with their own GRNs can largely skip this construction stage.
3. **Predictive category — 20 minutes.** For Gene Sets, ask: **do pathways whose genes are relatively highly expressed also appear as enriched labels somewhere among the TF target sets?** Follow the RNA-side and TF-target-side pathway lists, then compare their overlap. Contrast this functional concordance with `Gene ~ TFs` held-out RNA prediction. Examine precision, recall, and tested-feature coverage.
4. **Literature category — 15 minutes.** Discuss TF Markers (HPA and TF-Marker) and TF Pairs (IntAct and Europe PMC). Shared source identities make marker scores tie by construction; pair scores examine shared-target organization. Inspect CollecTRI overlap separately as a construction diagnostic.
5. **Mechanistic category — 15 minutes.** Inspect KnockTF activity scoring, the number of eligible/testable perturbations, and the dependence on signs and external experimental context. This tests recovery of reduced TF activity, not a full perturbation forecast.
6. **Return to the heatmap — 10 minutes.** Explain which evidence matters for the intended application, change class weights, and describe what would count as an independent evaluation of a student's network.

The teaching profile requests **nine columns spanning five task types**. CollecTRI reference agreement is also calculated but excluded from ranking because the CollecTRI baseline reuses it. Columns unavailable for any compared method are excluded from the common ranking and recorded explicitly. Read [TASKS.md](TASKS.md) for thresholds, denominators, coverage, and interpretation.

The default uses broad reference resources (`terms={}`): this mixed tumour panel has no single PBMC-like context. Broad marker or KnockTF evidence may therefore come from unrelated tissues. For a context-specific dataset, select appropriate terms explicitly with `--terms-json`; record this choice when comparing runs.

## Worked visual examples from this run

Each category now includes a figure that exposes the observations behind its scores, followed by a discussion question:

| Task | Data shown | What to interpret |
| --- | --- | --- |
| Gene Sets | Expression-side Hallmarks calls beside each network's pathway calls, with recovered and missed pathways | Why high precision can coexist with low pathway recovery |
| Omics | Observed versus predicted RNA for one real gene in held-out NCI60 lines, colored by tissue, plus its actual incoming TFs | Predictability, calibration, tissue structure and the limits of a single example |
| TF Markers | Actual marker-reference counts and example TF memberships | Why matched regulator identities necessarily produce tied scores |
| TF Pairs | One reported TF pair's target memberships, overlap counts and shared-target enrichment in each network | How target sharing becomes a proposed TF association |
| TF Scoring | Target responses and inferred activity for a real external KnockTF experiment | How signs, perturbation responses, significance and coverage affect a hit |

Examples are selected by documented rules. The Omics gene uses training information; its displayed test correlations do not select the example. The TF pair is explicitly selected as an illustration of reference-supported target sharing. The perturbation example uses knockdown magnitude and common coverage rather than activity-score contrast. Inspect the exported selection tables and try an explicit target, TF pair or experiment in the notebook.

The notebook includes saved inline figures for immediate viewing; rerun the cells to explore other examples. It also saves PNG/PDF copies with their underlying tables in `results/rna/examples/`. To prepare them without opening Jupyter, run:

```bash
python task_examples.py
```

This reuses the benchmark when its cache is valid and calculates the additional examples. It retains the existing scores and references. KnockTF examples are labeled as external perturbations; they do not represent perturbations measured in NCI60.

## Adapting the exercise

Use `python run_rna_benchmark.py --help` for all settings. `--max-genes`/`--max-tfs` can limit exploratory runs; **0 means all eligible features**. Gene filtering and correlation inference use training observations only. Add an existing RNA network with:

```bash
python run_rna_benchmark.py --grn MyGRN=my_network.tsv --output-dir results/my_grn
```

Its columns are `source`, `target`, `score`; identifiers must occur in the retained RNA genes, weights must be finite, and at least five nonself edges must remain. Supplied networks keep their own sizes rather than automatically matching baseline budgets.

For matching custom RNA, pass `--rna my_rna.h5ad`: use unique HGNC symbols, finite nonnegative log-normalized RNA, and real sample groups in `obs["celltype"]` (or `obs["group"]`) that support stratification. A gene-by-sample TSV with a `gene` column is also supported, with `--samples samples.tsv` providing `sample` and `tissue_code` columns. The default expression threshold is in NCI60 log2(FPKM+1) units; reconsider `--min-expression` for other normalization schemes. The teaching wrapper accepts at most 50 million sample-by-gene entries.

Students should record the expression preprocessing, gene identifiers, candidate TFs, network signs, and priors used during inference. A network inferred using all observations has already seen the later test set: fitting a prediction model on a training split cannot undo that leakage. New donors or independently collected samples are preferable when the question is transfer across datasets. A network's size and available reference coverage should accompany its score.

[rna_benchmark.py](rna_benchmark.py) contains the RNA-only helpers; [run_rna_benchmark.py](run_rna_benchmark.py) is the command-line entry point. Key outputs are `networks/*.tsv`, `inference_split.tsv`, `feature_info.tsv`, `evaluation_all.tsv`, `evaluation_ranked.tsv`, `coverage.tsv`, `excluded_from_ranking.tsv`, `ranking_task.png/.pdf`, `category_*.png/.pdf`, and `run_info.json`. Keep gretapy **0.0.3** pinned: a local compatibility context repairs its task-label routing and excludes undefined constant-response correlations. It also excludes KnockTF experiments with fewer than three nonzero measured regulon targets, as required by native ULM after zero-response filtering. Native ULM and F₀.₁ formulas are retained; precision uses actually tested experiments and recall retains the original eligible-experiment universe.

## Validation and instructor notes

The complete RNA benchmark and all **34 notebook cells** were executed successfully, including nine inline figures. The worked-example pathway/marker/pair counts and perturbation hit counts were reconciled with the saved native scores; the Omics example was checked against native predictions and the held-out split. **Twelve benchmark checks passed**, including full-prior retention, matching degree budgets, holdout isolation, cache reuse/tamper detection, and eligibility/recall handling. From the repository root, with the environment active:

```bash
python -m unittest discover -s GRN_comparison -v
```

Reference-dependent integration checks skip themselves until the relevant resources are cached.

The default result needs a calibration warning: Pearson leads held-out expression prediction and CollecTRI leads TF-loss activity scoring, but Random leads **native gene-set F₀.₁ despite much lower pathway recall**. This is not better functional coverage. For Hallmarks, Random recovers 10 of 24 expression-supported pathways, while CollecTRI recovers all 24; F₀.₁ rewards Random's higher precision. The fixed 20,000-gene enrichment background also differs from the 7,530 retained candidate genes. Restricting pathway memberships and the background to retained genes materially changes the scores, but repeated random controls show that this adjustment alone does not calibrate the task. See [TASKS.md](TASKS.md) for interpretation.

TF Marker columns are regulator-selection diagnostics and tie because the baselines share their source TFs. They provide no evidence about differences in wiring. The native overview retains these tasks for reproduction; its composite ordering should not be presented as a calibrated ranking of GRN quality. A comparison focused on wiring should keep marker checks separate and assess gene-set scores against repeated appropriate null networks before using them in an aggregate. Category panels recompute their row order and can rescale colors; read their labels and colorbars before comparing them.

The original introductory [01_grn_comparison.ipynb](01_grn_comparison.ipynb) and [NCI60_EXTENSION.md](NCI60_EXTENSION.md) remain available. The previous synthetic RNA+ATAC example is preserved as [03_rna_atac_preview.ipynb](03_rna_atac_preview.ipynb), with its original [preview guide](RNA_ATAC_PREVIEW.md) and [task guide](RNA_ATAC_TASKS.md). A dedicated RNA+ATAC teaching session can build on that preview afterward.
