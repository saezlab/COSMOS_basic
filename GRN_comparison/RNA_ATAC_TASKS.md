# What does a GRN benchmark tell us?

This practical asks whether a gene regulatory network (GRN) is useful under several tests. A network can agree with published interactions, explain expression variation, and still predict a perturbation poorly. No single task establishes that a network is biologically correct.

The live exercise uses gretapy's small, **synthetic PBMC-like dataset**, with paired RNA and chromatin accessibility measurements and B-cell, T-cell, and monocyte labels. It runs selected evaluation tasks. A separate panel uses the authors' **precomputed published benchmark scores** to illustrate the wider framework. These are not recalculated in class. Do not combine toy and published results as though they came from the same experiment. A toy ranking is not evidence that one method generally outperforms another.

## Where do the networks come from?

For a simple data-driven baseline, start with candidate transcription factors (TFs) and ask which genes have expression that varies with each TF across the training cells. Pearson correlation measures approximately linear association; Spearman correlation compares ranks. Retain sufficiently strong associations, remove self-edges, and attach accessible regions overlapping the target genes' annotated promoters. This gives a table of `source` TF, `cre` region, `target` gene, and `score`.

This is a convenient network to evaluate. Correlation does not establish the direction of regulation: the TF-to-gene direction comes from our choice of candidate regulators. A positive or negative correlation is also not proof of activation or repression. Cell identity, shared upstream regulators, and technical effects can produce correlations. The attached promoter region is an annotation-based assignment, not evidence that this TF binds there.

The CollecTRI baseline starts from curated TF–target interactions compatible with the measured genes and promoter-overlapping regions. Both baselines produce the same table format, while using different evidence. Students bringing an existing single-cell RNA + ATAC GRN can replace this stage and concentrate on evaluation.

## Precision, recall, and the score

Each task defines its own successful predictions. Precision asks, “Of the predictions evaluated, how many passed this test?” Recall asks, “How much of the reference or measured universe did the successful predictions cover?” For edge-overlap tasks, the conventional definitions are

\[
P=\frac{TP}{TP+FP},\qquad R=\frac{TP}{TP+FN}.
\]

For prediction and simulation tasks, the counted units can instead be genes, perturbation experiments, or steady states. Always identify the unit before interpreting a score.

Gretapy combines precision and recall using

\[
F_{0.1}=\frac{(1+0.1^2)PR}{0.1^2P+R}
       =\frac{1.01PR}{0.01P+R}.
\]

Higher is better under the task's definition. This score strongly emphasizes precision while still penalizing very low coverage. The evaluation documentation describes this as weighting precision “ten times” more; the implemented formula uses **β squared**, so β is not a simple tenfold weight. Use the formula and inspect precision and recall separately. In reciprocal harmonic form, their relative weights are 100:1.

## The eleven tasks

The four classes group related questions. Within a task, several reference databases or prediction subtasks can provide separate columns in the final figure.

### Predictive: does the network explain molecular or functional patterns?

| Task | Inputs and what counts as a successful hit | How to interpret it |
|---|---|---|
| **Omics** | Network-selected predictors and paired measurements. `Gene ~ TFs` predicts target RNA from parent TF RNA; `Gene ~ CREs` predicts RNA from linked region accessibility; `CRE ~ TFs` predicts accessibility from TF RNA. A hit is a molecular feature whose predicted and observed held-out values have a significant positive correlation. | Tests useful predictive structure. Predictability can arise from cell identity or shared regulation; it does not demonstrate a direct causal edge. The three prediction directions require different inputs. |
| **Gene Sets** | Expression, pathway gene sets such as Hallmarks or Reactome, and each TF's target set. A hit is a pathway both enriched among network targets and classified as active in the expression data. | Tests agreement in broad functional coverage. It does not test individual edges or predict a pathway's response to perturbation. Many different networks can recover the same broad pathways. |

In the Omics implementation, gretapy fits a common prediction model using at most five incoming predictors per feature, ordered by absolute edge score. It judges predictions with a held-out Spearman correlation greater than 0.05 and a multiple-testing-adjusted p-value below 0.05. Precision is the fraction of tested features that succeed; recall is the fraction of all measured features that succeed. The pinned implementation excludes zero-valued responses. Features with too few usable test observations are not tested. A high precision can therefore coexist with narrow coverage.

Our baselines share one promoter per target. Consequently, `Gene ~ CREs` largely reflects which genes are included and their promoter accessibility–RNA association; it cannot distinguish the quality of TF assignments on its own.

Gene Sets compares expression-derived pathway activity with enrichment among TF target sets. The expression side requires positive activity with adjusted p-value below 0.01 in more than 20% of samples; network enrichment also uses an adjusted p-value below 0.01. Small regulons and a tiny measured gene universe can give few hits. A zero is an opportunity to inspect coverage and thresholds, not a reason to adjust them until a preferred method wins.

### Genomic: are the proposed regions and links supported by annotations?

| Task | Inputs and what counts as a successful hit | How to interpret it |
|---|---|---|
| **TF Binding** | TF–CRE pairs, genomic coordinates, and binding collections such as ChIP-Atlas, ReMap 2022, or UniBind. A hit overlaps a reference binding region for the corresponding TF. | Supports physical binding plausibility. Binding alone does not establish functional regulation of the linked gene. Results depend on which TFs and biological contexts have reference coverage. |
| **CREs** | Network regions and annotations of regulatory elements, promoters, conservation, or other genomic features. Usually a hit is a network region overlapping an annotation. For **ENCODE Blacklist**, the desired outcome is avoiding blacklisted regions. | Asks whether selected regions have the property represented by that annotation. Promoter overlap, conservation, and avoidance of technical artefacts are different properties. |
| **CRE to Gene** | Region–target assignments and an external region–gene association resource, such as eQTL Catalogue. A hit overlaps a reference region associated with the same target gene. | Supports the region–gene link. It does not establish the TF–region link, and association need not imply direct regulation. |

These tasks need CRE coordinates in the correct genome assembly and matching measured ATAC peaks. A TF–target list alone is insufficient. Reference peaks from unrelated tissues may be poor evidence for the cell types being studied. Unrepresented TFs or regions can be outside a resource's evaluation coverage rather than experimentally disproven.

**Construction can predetermine a score.** Our baselines attach promoter-overlapping regions. Their performance against a promoter annotation partly reflects that construction rule. It is useful for checking the pipeline, but is not independent validation of enhancer discovery. Record all resources used to infer, filter, or annotate a network before choosing evaluation resources.

### Literature: does the network recover previously reported biology?

| Task | Inputs and what counts as a successful hit | How to interpret it |
|---|---|---|
| **TF Markers** | Source TF identities and tissue- or cell-type-relevant marker collections, such as HPA or TF-Marker. A hit is a network TF present in the selected marker collection. | Tests regulator selection, not target assignment. Networks with identical source TFs necessarily tie, even when their edges differ completely. |
| **TF Pairs** | TF target sets and reference TF interaction pairs, such as Europe PMC or IntAct. The network first proposes TF pairs from significantly overlapping target sets; a hit is a proposed pair present in the reference. | Tests whether shared-target organization agrees with known TF associations. It does not simply look for direct TF→TF edges. A physical interaction and shared target regulation are related but distinct properties. |
| **Reference GRN** | Directed TF–target pairs and a curated reference, here CollecTRI. A hit is the same directed pair in both. | Measures agreement with existing knowledge. Missing reference edges are not established biological negatives, and well-studied TFs have more opportunities for support. |

**CollecTRI against itself is circular.** A CollecTRI-derived network should recover CollecTRI edges. A high score in this column demonstrates retention of the input prior, not superior independent accuracy. Keep it as an explicitly labelled diagnostic or remove it from an independent comparison. The same concern applies when a student's method used CollecTRI during inference. Other literature resources may also share underlying publications, so different database names do not guarantee independent evidence.

### Mechanistic: can the network account for functional outcomes?

| Task | Inputs and what counts as a successful hit | How to interpret it |
|---|---|---|
| **TF Scoring** | Signed TF regulons and external TF knockdown expression profiles, such as KnockTF. A hit is an eligible knockdown for which the perturbed TF's inferred activity decreases significantly. | Tests whether targets and signs detect a known loss of TF activity. It does not require forecasting every target's expression change. Regulons need sufficient measured target overlap. |
| **Perturbation Forecasting** | Baseline RNA, network topology, and external knockdown responses. The method fits regulatory coefficients, simulates TF loss, and counts perturbations with significant positive agreement between predicted and observed responses. | Tests response prediction under a particular fitted model. The coefficients are re-estimated: this is not a direct test of the original edge weights. Context mismatch between baseline and perturbation experiments can limit interpretation. |
| **Steady State Simulation** | A signed TF-to-TF subnetwork and cell-type marker TFs. Boolean rules generate steady states. A successful state matches exactly one cell type; recall asks how many cell types have at least one matching state. | Tests consistency with a simplified model of stable cell identities. Matching a marker pattern does not prove that the biological system follows the assumed Boolean dynamics. |

These tasks can need substantially more coverage than our small networks provide. TF Scoring needs at least three usable regulon targets. Forecasting requires sufficient incoming regulation of the perturbed TF and sufficient shared response genes. Boolean evaluation needs a suitable signed TF subnetwork and distinguishable marker sets. An unavailable task must be reported as unavailable, rather than presented as a failed biological prediction.

## Reading `gt.pl.ranking(res, level="task")`

Read the figure in three steps:

1. **Identify the column.** Class and task headers group the columns; the lower labels identify the database or prediction subtask. Several columns can belong to one task.
2. **Read color and number separately.** Color represents mean F₀.₁, with higher values better. The printed number is the rank **within that column**: 1 is best, ties receive average ranks, and missing values appear as `-`. Rank 1 can still correspond to a low absolute score.
3. **Inspect the ordering and coverage.** Rows are ordered by a composite score, not by how often a method is ranked first. Gretapy averages databases within tasks, tasks within classes, combines classes using their weights within each dataset, and finally averages datasets. The default gives each available class equal weight. This hierarchy stops a class with many databases automatically dominating the score.

Missing classes are excluded and the remaining weights are renormalized. Consequently, two methods with different available tasks may be compared on different evidence. Check the coverage table and compare a common set of tasks; do not replace missing results with zero. Changing `metric_weights` changes row ordering, not the underlying task scores. At task level, a column's color averages its available dataset scores, while the overall ordering uses the hierarchy above.

## Keeping a comparison interpretable

Split cells **before** learning correlations, selecting variable genes, choosing thresholds, or fitting preprocessing that uses expression. A held-out prediction test is optimistic if the test cells already helped select the network's edges. Choosing features or tuning the network on all cells can leak information even when the final predictor is fitted only on training cells. For real single-cell work, donor or sample splits may be more informative than random cell splits when the intended use is transfer to new samples.

Keep network size, regulator and target universes, usable observations, reference coverage, and score signs visible. Equal edge budgets help, but do not make all networks equally predictive or equally testable. A positive importance score is not automatically an activating effect; tasks requiring biological signs need an appropriate signed network.

For a future evaluation of your own GRN, supply the corresponding expression/accessibility data, matching gene identifiers and genome build, cell-type labels, and a clear account of which priors influenced inference. Decide the biological question and applicable tasks before interpreting the final ranking.

## Discussion exercises

1. Find a method that performs differently across two tasks. Explain the biological questions behind those columns before proposing a reason for the difference.
2. Remove the CollecTRI reference column. Does the ranking change, and which comparison becomes less circular?
3. Explain why a perfect promoter score would be unsurprising for a promoter-constructed baseline.
4. Find a tied column or a zero. Use precision, recall, and coverage to distinguish identical predictions, few tested features, and unsuccessful tests.
5. Change the class weights. Describe the scientific preference this encodes; do not describe the resulting order as a newly discovered universal winner.
6. Propose an evaluation for your own RNA + ATAC network, naming one independent resource and one resource that would be circular to reuse.

## Sources and implementation scope

This guide follows the [official evaluation tutorial](https://gretapy.readthedocs.io/en/latest/notebooks/eval.html) and the implementation tested here, **gretapy 0.0.3**, at [source commit d5137ee](https://github.com/saezlab/gretapy/tree/d5137ee35585885fb6e0a029d808d83cc4a52f0c). Relevant primary implementation references are the [prediction tasks](https://github.com/saezlab/gretapy/blob/d5137ee35585885fb6e0a029d808d83cc4a52f0c/src/gretapy/tl/_predictive.py), [genomic tasks](https://github.com/saezlab/gretapy/blob/d5137ee35585885fb6e0a029d808d83cc4a52f0c/src/gretapy/tl/_genomic.py), [literature tasks](https://github.com/saezlab/gretapy/blob/d5137ee35585885fb6e0a029d808d83cc4a52f0c/src/gretapy/tl/_prior.py), [mechanistic tasks](https://github.com/saezlab/gretapy/blob/d5137ee35585885fb6e0a029d808d83cc4a52f0c/src/gretapy/tl/_mechanistic.py), [score formula](https://github.com/saezlab/gretapy/blob/d5137ee35585885fb6e0a029d808d83cc4a52f0c/src/gretapy/tl/_utils.py), and [ranking aggregation](https://github.com/saezlab/gretapy/blob/d5137ee35585885fb6e0a029d808d83cc4a52f0c/src/gretapy/tl/_ranking.py). The [current main branch](https://github.com/saezlab/gretapy) may evolve; use the pinned version and this practical's coverage report when reproducing the session.
