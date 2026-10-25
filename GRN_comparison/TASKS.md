# Interpreting the RNA-only GRETA benchmark

Start with the overall heatmap in [02_greta_evaluation.ipynb](02_greta_evaluation.ipynb), then use this guide as you examine each evaluation category. Pearson, CollecTRI, and Random use NCI60 RNA with identical TF identities and outgoing edge budgets. Predictive, literature, and perturbation tests ask different questions; success in one need not imply success in another.

## Read the overall figure first

```python
gt.pl.ranking(res, level="task")
```

The class/task headers organize the figure; the lower labels identify each database or prediction subtask. One task can have several database columns.

- **Color is F₀.₁.** A brighter color means a higher score under that task's definition. Inspect absolute scores as well as relative ranks.
- **The printed number is the rank within a column.** Rank 1 is best; ties receive average ranks. First place among three weak results remains a weak result.
- **Row order uses a composite score.** GRETA averages databases within tasks, tasks within classes, then combines available classes with their weights. With several datasets it finally averages dataset scores. This prevents database count determining class weight.

Available classes have equal default weight; the absent Genomic class contributes none. Changing weights changes row ordering, not task results. Missing classes are omitted and weights renormalized. The wrapper ranks only columns available for every compared network and reports exclusions.

**Calibration limitation in this example:** the native overview is a reproduction of GRETA's aggregation, not an established ranking of wiring quality. TF Markers cannot discriminate networks with identical TFs, and native Gene Sets can assign high scores to random wiring. Keep the raw results visible, inspect precision/recall, and assess appropriate repeated null controls before selecting tasks for a biological comparison. Choose inclusion rules from the scientific question, rather than from which method wins.

**Always inspect coverage.** Unavailable tests are missing values; eligible tests with no successful predictions can legitimately score zero. Equal availability does not imply equal numbers of tested genes or perturbations.

## What does F₀.₁ mean?

For reference overlap, precision is supported predicted items / predicted items; recall is recovered reference items / reference items. Prediction and activity tasks count genes or experiments instead. Identify the counted unit first.

\[
F_{0.1}=\frac{1.01PR}{0.01P+R}.
\]

Higher is better under each task's definition. This harmonic combination strongly favors precision (reciprocal weights 100:1, because β is squared), while penalizing extremely low recall. F₀.₁ is not the fraction of biologically correct edges.

## Briefly: where did the data-driven GRN come from?

Across the **40 training cell lines**, correlate each eligible TF's expression with every eligible gene. Retain the strongest absolute Pearson correlations, with each TF receiving as many targets as it has in the eligible CollecTRI network. CollecTRI retains its full eligible links; Random uses the same regulator/degree budgets. Spearman correlation is optional.

Gene filtering and edge selection exclude the **20 held-out lines**. Training expression ≥1 in ≥80% of lines and nonconstant expression retain 7,530 genes; 195 TFs produce 6,871 edges per baseline. There is no default gene/TF cap. Candidate TFs and degrees use CollecTRI, so inference is prior-guided.

Correlation is association; TF→gene direction is imposed by choosing TF sources, and signs are not established activation/repression effects. Tissue identity and shared regulators can generate correlations. NCI60 contains bulk cancer cell lines: `celltype` holds authentic tissue-code strata.

## Predictive: does the network explain RNA or functional patterns?

### Gene Sets — Hallmarks, PROGENy, Reactome

**Question: do pathways whose genes are relatively highly expressed also appear as enriched labels somewhere among the TF target sets?**

GRETA performs two different enrichment analyses and compares their **lists of pathway names**. For the unweighted Hallmarks example:

1. **RNA → List A.** For each cell line, ULM tests whether pathway genes have higher average log-expression than other measured genes. Its score is the t-statistic of that relationship. GRETA adjusts p-values across all cell-line × pathway tests, requires a positive score and adjusted p<0.01, then retains pathways passing in more than 20% of lines. At least 13 of our 60 lines must pass. This produces **24 RNA-enriched pathways**, using all 60 lines. Relative expression within a line does not establish activation compared with a control condition.
2. **TF target sets → List B.** For each TF, a one-sided Fisher test asks whether its targets contain more pathway genes than expected by chance. BH correction is applied within that TF's pathways with nonzero overlap. Pool pathways with adjusted p<0.01 in **any TF regulon**, ignoring GRN edge signs and strengths. Count each pathway once, regardless of how many TFs support it. Pearson's List B contains **32 pathways**.
3. **Compare the names.** Shared pathways are TP; GRN-only pathways are FP; RNA-only pathways are FN. Pearson has **18 TP, 14 FP and 6 FN**, giving precision 18/32, recall 18/24 and F₀.₁ ≈0.564. The counted items are pathways, rather than genes, TFs or regulatory edges.

**Actual example — Coagulation:** it passes the RNA-side test in **30/60 lines**, so it enters List A. CollecTRI's HIF1A regulon contains **6 Coagulation genes among 190 targets**; under the native background its adjusted enrichment p-value is **0.00324**. Coagulation therefore enters CollecTRI's List B and counts as a TP. No Pearson regulon passes for it, so it counts as a Pearson FN.

Hallmarks and Reactome use unweighted memberships. PROGENy follows the same two-list comparison, but its RNA-side ULM uses signed pathway weights rather than the simple member-versus-other-gene contrast. See the [ULM documentation](https://decoupler.scverse.org/en/latest/api/generated/decoupler.mt.ulm.html).

A high score supports **functional concordance**: RNA-enriched pathway labels also appear somewhere among the TF target sets. It does not identify which TF correctly controls a pathway or validate individual edges. Distinct networks can recover the same pathway labels. CollecTRI and pathway collections may also share studies, limiting independence.

The pinned native enrichment uses a fixed **20,000-gene background**, not this dataset's measured/eligible universe. Our random targets can only come from the 7,530 retained genes (excluding the source TF). A matched enrichment null must restrict each pathway's membership as well as the background size to the eligible targets. Changing only the background size is insufficient. The rationale is the set of genes that could enter the tested list, consistent with [Gene Ontology's enrichment guidance](https://www.geneontology.org/docs/go-enrichment-analysis/). Expression-side activity uses the supplied RNA dataset, so this is a functional concordance task, not a held-out prediction test.

**Why can Random win?** The task pools pathways significant in any of 195 regulons, then scores agreement with expression-supported pathways. It does not check which TF should regulate which pathway, nor regulatory signs. Multiple regulons provide repeated opportunities for chance enrichment; native BH adjustment is within each regulon, without an additional correction for the final union. F₀.₁ strongly rewards precision even when only a few expression-supported pathways are recovered. In the saved Hallmarks result:

| Network | Supported pathways recovered | Precision | Recall | F₀.₁ |
| --- | ---: | ---: | ---: | ---: |
| Pearson | 18 / 24 | 0.563 | 0.750 | 0.564 |
| CollecTRI | 24 / 24 | 0.511 | 1.000 | 0.513 |
| Random | 10 / 24 | 0.667 | 0.417 | 0.663 |

A diagnostic restricting pathway memberships and the background to eligible genes lowers the saved Random Hallmarks score to zero. The exact per-TF universe has 7,529 targets because self-edges are excluded. However, across ten additional random seeds with the same per-TF degree budgets, five still exceed Pearson's Hallmarks F₀.₁ under that adjusted background, while recovering at most six of the 24 expression-supported pathways. The background matters, but correcting it alone does not validate this score for comparing wiring. These are sensitivity checks, not replacements for the saved native scores or proof of general method performance. Report repeated null distributions and counts, and treat this task as a functional-concordance diagnostic until calibrated for the intended comparison.

To reproduce this diagnostic after running the default NCI60 benchmark, use `python diagnose_gene_sets.py`. The [audit script](diagnose_gene_sets.py) verifies all nine native network/database scores, compares the exact eligible background, and evaluates ten additional Hallmarks random networks. It writes counts and scores to `results/gene_set_audit/`, separately from the benchmark results, without refitting the Omics prediction models.

**Inspect:** precision, recall, and agreement across gene-set collections.

**Worked figure:** labeled Hallmarks pathways span different expression/network call patterns. Follow one expression-supported pathway across the columns, then relate the recovered, additional and missed pathways to precision and recall. Counts and scores use all 50 pathways, whose full call table is exported; the displayed subset illustrates the mechanism rather than sampling method performance. A binary pathway call does not reveal the number of correctly assigned TF–target edges.

### Omics — `Gene ~ TFs`

**Question:** do the TF predictors chosen by the network predict target RNA in held-out cell lines?

For each target, native GRETA selects at most five incoming TFs by absolute edge score, fits an XGBoost regressor on training observations, and predicts test expression. A successful feature has positive held-out Spearman correlation >0.05 and Benjamini–Hochberg adjusted p<0.05. Precision is successful genes / tested genes; recall is successful genes / the full supplied RNA feature universe.

The split is stratified by original tissue code, with `test_size=0.33`, `random_state=42`; inference uses that same training partition. Only 20 test lines limit precision of the result.

The pinned evaluator **removes observations whose target expression is exactly zero**, although zero can be a genuine measurement here, and requires at least ten remaining test observations. Constant predictions are not tested; our compatibility guard also excludes constant held-out responses with undefined correlations. Read tested-gene and exclusion counts alongside the score. High precision with low recall can mean that only a restricted subset is predictable.

**Interpretation:** useful predictive structure can reflect tumour/tissue differences or shared upstream signals. The model refits coefficients and does not validate individual causal TF→target edges. Same TF outdegrees do not guarantee the same targets or incoming predictor counts.

**Worked figure:** observed versus predicted target RNA in the held-out NCI60 lines, colored by tissue, alongside the actual incoming TF predictors. The gene is selected using training correlations among targets shared by the networks. Read the selection table: this is an illustrative training-selected association. Distance from the identity line indicates prediction calibration; the displayed Spearman coefficient measures ordering. Its raw correlation does not reproduce the full benchmark's BH-adjusted success decision. Try another gene with `plot_omics_example(result, target="GENE_SYMBOL")`.

## Literature: does the network recover reported biology?

### TF Markers — HPA, TF-Marker

**Question:** are selected source TFs represented in the marker reference?

The counted items are TF identities. Precision asks what fraction of network TFs are markers; recall asks what fraction of eligible reference markers were recovered. This task evaluates **regulator selection**, not target assignment.

All teaching baselines have identical TF sources, so their marker scores tie by construction. These columns cannot distinguish differences confined to edges. They are useful regulator-selection/coverage diagnostics, but should be separate from a primary comparison conditioned on a fixed TF list. Keeping a constant task in the native aggregate also dilutes the contribution of TF Pairs within the Literature class. The default broad references include many tissues; context-specific conclusions require suitable explicit term filters.

**Worked figure:** actual source-TF counts listed and absent in HPA/TF-Marker, plus example memberships for named TFs. All methods should have the same counts. Absence from the resource is an annotation result, not a confirmed biological negative.

### TF Pairs — IntAct, Europe PMC

**Question:** do TFs with significantly overlapping target sets correspond to reported TF associations?

GRETA first proposes pairs through shared-target enrichment (adjusted p<0.01), then compares these pairs with the resource. It filters to resource-covered TFs; the overlap background depends on the network's target universe. Precision is supported proposed pairs / proposed pairs; recall is recovered eligible reference pairs / eligible reference pairs.

This does **not** simply check direct TF→TF edges. It tests organization of target sets against reported TF associations. Physical interaction, published association, and co-regulation are related but different observations. Coverage and publication bias matter; a zero can mean no inferred significant pairs or no supported pairs.

**Worked figure:** one reference-supported pair followed across all networks, with actual target memberships, overlap counts, and enrichment q-values. The pair is selected to illustrate the mechanism and is not representative of all pairs. Displayed targets can be a subset for readability; the statistics use the complete target sets and native multiple-testing family. Ask why a large overlap might still be unsurprising for two large regulons. Inspect a specific pair with `plot_tf_pair_example(result, pair=("TF_A", "TF_B"))`.

### CollecTRI reference agreement — diagnostic only

Directed TF–target overlap measures recovery of prior knowledge. Missing reference edges are not established biological negatives. A CollecTRI-derived baseline scoring against CollecTRI reuses its construction evidence, so this column is shown separately and **excluded from the main ranking**. This also applies to student networks using CollecTRI during inference.

## Mechanistic: does a signed regulon recover known TF loss?

### TF Scoring — KnockTF (scoring)

**Question:** after an experimentally measured TF knockdown, does its regulon indicate decreased TF activity?

GRETA selects external KnockTF experiments with the perturbed TF measured in the RNA universe and knockdown log fold change <−0.5. A regulon needs at least three targets with nonzero measured responses in each perturbation profile, because native ULM removes zero-response genes. The compatibility helper excludes untestable experiments before native ULM. A hit has negative activity and adjusted p<0.05. Precision counts successful / actually tested experiments; recall retains the original eligible-experiment universe, including experiments not testable by that network.

Inspect eligible, testable, and successful experiment counts. Broad KnockTF coverage includes contexts different from NCI60: this is **external perturbation evidence**, with possible context mismatch.

A high score supports recovery of reduced activity, not forecasting of all downstream expression changes. CollecTRI has curated regulatory signs; Pearson has association signs, and the Random control is arbitrary. Results reflect sign quality as well as target selection.

**Worked figure:** a real KnockTF experiment's target log fold changes, network signs, and inferred TF activity. The default selection uses the median own-TF knockdown among experiments supported by all networks, independently of activity. The adjusted p-values retain native BH correction across supported experiments within each network. A positive edge with a decreasing target and a negative edge with an increasing target both agree with a TF-loss footprint; ULM uses the broader measured background as well. This is an external experimental context, not an NCI60 knockdown. Use `plot_knocktf_example(result, experiment="EXPERIMENT_ID")` to inspect another eligible experiment.

## What is outside this RNA-only lesson?

No genomic regions are fabricated. TF Binding, CRE annotation, and CRE-to-gene evaluation require suitable genomic evidence and coordinate conventions; they belong in the later RNA+ATAC session. Perturbation forecasting and Boolean steady-state simulation can use RNA-based networks, but are omitted here to keep this teaching benchmark manageable and avoid introducing additional dynamical assumptions. Their omission is a scope choice, not proof that they require ATAC.

## Discussion exercises

1. Find a disagreement between columns. State both biological questions before explaining the difference.
2. Explain why TF Marker scores tie and why equal outgoing degrees do not ensure equal Omics coverage.
3. Inspect a zero or missing score. Distinguish failed tests, narrow coverage, and unavailable evidence.
4. Explain why CollecTRI self-agreement stays outside the ranking. Name another potential source of reused evidence.
5. Change class weights and describe the scientific preference this encodes.
6. Propose an evaluation for your own network: document construction data/priors and choose an independent test. A GRN learned on all supplied observations has already seen the nominal Omics test set; an evaluation split cannot undo that leakage.

This guide follows the [official evaluation tutorial](https://gretapy.readthedocs.io/en/latest/notebooks/eval.html) and pinned **gretapy 0.0.3**. Primary implementation details are in its [predictive tasks](https://github.com/saezlab/gretapy/blob/d5137ee35585885fb6e0a029d808d83cc4a52f0c/src/gretapy/tl/_predictive.py), [literature tasks](https://github.com/saezlab/gretapy/blob/d5137ee35585885fb6e0a029d808d83cc4a52f0c/src/gretapy/tl/_prior.py), [mechanistic tasks](https://github.com/saezlab/gretapy/blob/d5137ee35585885fb6e0a029d808d83cc4a52f0c/src/gretapy/tl/_mechanistic.py), and [ranking aggregation](https://github.com/saezlab/gretapy/blob/d5137ee35585885fb6e0a029d808d83cc4a52f0c/src/gretapy/tl/_ranking.py). The practical preserves the native algorithms while repairing version-specific routing/undefined-correlation behavior and reporting coverage.
