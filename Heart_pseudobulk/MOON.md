# RNA-only MOON: upstream hypotheses from TF activities

`03_moon.R` asks **which upstream signaling activities are compatible with a
cell type's inferred TF activity pattern?** It adapts the NCI60 MOON practical
to a case with no measured upstream layer and no metabolomics.

The inputs are the unthresholded `tf_scores.csv` and `rna_scaled.csv` written
by `01_tf_activity.py` in `results/tf_activity/` (CollecTRI) and
`results/tf_activity_scenicplus/` (SCENIC+). Both use exactly the same RNA
preprocessing and gene-wise Z-scores across the six cell types.

| Quantity | Role in this run |
| --- | --- |
| Signed CollecTRI or SCENIC+ ULM TF scores | Downstream activity footprints; full-prior runs use `abs(score)>2`, while the controlled comparison uses all finite shared TF scores |
| Upstream measurements | Unspecified: `upstream_input=NULL` |
| Gene-wise RNA Z scores | Node annotations only; do not substitute RNA abundance for protein activity |
| SCENIC+ network | Supplies one alternative TF-target footprint definition; the protein/complex signaling PKN is unchanged |
| Identifier domain | Original gene symbols and named signaling complexes from the packaged PKN |

These TF scores describe **relative transcriptional activity across the six
pooled cardiac cell types**. They are not donor-level differential statistics,
and a positive activity does not mean activation relative to a healthy or
unstimulated control. Cell types contain differing donor mixtures. The RNA
baseline and pseudobulk construction remain those documented for script 01.

## Install the R environment once

An existing environment that runs the repository's NCI60 MOON script can run
this script. R is separate from the Python virtual environment. In a fresh R
session, install dependencies explicitly:

```r
install.packages(c("BiocManager", "remotes", "jsonlite"))
BiocManager::install(c("CARNIVAL", "decoupleR", "GSEABase"),
                     ask = FALSE, update = FALSE)
remotes::install_github(
  "saezlab/cosmosR@b2ed195738479c4ad204cace37fd11ae2fd02bb8",
  dependencies = NA, upgrade = "never")
```

That commit is the tested cosmosR 1.21.2 version. The script records the actual
installed version and Git commit in every run. cosmosR imports CARNIVAL as an R
dependency, but this MOON workflow uses no CARNIVAL optimization and needs no
CPLEX, CBC, or other external optimization solver. The analysis never installs
or upgrades packages while running.

## Run from the repository root

After script 01:

```bash
Rscript Heart_pseudobulk/03_moon.R
```

This analyzes **Ventricular Cardiomyocyte** under both priors, including the
shared-TF control. To choose a different column in the TF score matrix or run
all six cell types and generate comparison figures:

```bash
Rscript Heart_pseudobulk/03_moon.R --cell-type "Myeloid"
Rscript Heart_pseudobulk/03_moon.R --all
python Heart_pseudobulk/04_compare_moon.py
```

Single-prior runs remain available with `--mode collectri` and
`--mode scenicplus`. `--help` lists input/output directories, layer count,
thresholds, and compression options. A custom `--input-dir` retains the previous
single-run behavior unless `--mode compare` is specified explicitly.

Full-prior outputs are `results/moon/<cell_type>/` and
`results/moon_scenicplus/<cell_type>/`. Controlled runs are saved separately in
`results/moon_comparison/controlled/{collectri,scenicplus}/<cell_type>/`.
Comparison tables and ledgers are in `results/moon_comparison/`, and script 04
saves figures and annotated tables in its `figures/` subfolder. Re-running
updates the selected examples. Rerunning a single cell type refreshes comparison
tables for that requested type; use `--all` to regenerate the six-type figures.

## Compare TF activity sources while holding coverage constant

CollecTRI scores 180 TFs and SCENIC+ scores 102, with 40 shared. Their full-prior
`abs(score)>2` inputs also differ between cell types. Comparing those results
alone mixes the footprint estimates with the number and identity of TF inputs,
the pruned PKN and potentially the scoring levels. Those full-prior runs are
retained as a supplementary view.

The main controlled analysis uses TFs finite in **both methods in all supplied
profiles** and represented in the signaling PKN. It includes every such score,
without an activity cutoff, and preserves genuine zero values. Forty eligible
TFs become **39 observable inputs** because KLF2 has no usable incoming
signaling support after the same pruning. The final **12,170-edge scoring
network**, compression mapping, input identities and scored-node levels are
identical for both methods and all six profiles; the script asserts those
invariants before interpreting differences. Each profile has 3,368 inferred
original-node IDs, plus the 39 observed TFs.

`moon(..., downstream_cutoff=-Inf)` explicitly preserves observed level-0 inputs
whose value equals zero. The package's default cutoff would otherwise omit
them from its returned table. The controlled graph is prepared once from input
identities; each run supplies the actual method-specific numeric activities.

Both comparison CSVs retain the union of original node IDs. A score difference
is defined only when the node is **inferred in both runs at the same level**.
Observed level-0 TFs, absent scores and changed-level full-prior results have
missing differences and an explicit status. Missing results are never zeroed.

| Column | Meaning |
| --- | --- |
| `delta_score` | SCENIC+ score minus CollecTRI score; a shift toward a more positive or negative inferred activity |
| `delta_abs_score` | `abs(SCENIC+) − abs(CollecTRI)`; positive means larger absolute alignment under SCENIC+-derived footprints |
| `rank_abs_collectri`, `rank_abs_scenicplus` | Absolute-score ranks within the same cell type and scoring level, among comparable original node IDs |
| `rank_gain_scenicplus` | CollecTRI rank minus SCENIC+ rank; positive means a higher rank with SCENIC+ |
| `source_scoring_collectri`, `source_scoring_scenicplus` | Exact scoring/compression identities used to trace expanded original IDs |

Ranks over original IDs can include tied members of the same inference group.
The main figures show one representative per indistinguishable footprint, as
described below. These are descriptive scores and ranks, without a significance
test for their difference. Different regulon membership, signs, target counts
and relative weights are part of the comparison; fixed TF identities do not
remove those differences or establish which GRN is more accurate.

## What the script does

1. Loads the bundled `cosmosR::meta_network`, preserving a raw copy and package
   provenance. Keeps its protein/complex signaling layer, excluding metabolite
   nodes (`Metab__...`) and reaction-specific enzyme copies (`Gene<digits>__...`).
   It cleans self-loops, duplicate pairs and conflicting signs using the package
   cleaner, with before/after counts in `pkn_preparation.csv`.
2. Chooses covered downstream TF activities (the full-prior cutoff or the
   globally shared unthresholded identity set),
   then keeps nodes up to four directed steps upstream of those observations.
   **It does not prune by reachability from an empty upstream vector.** It also
   does not expression-filter the PKN: the 4,070 supplied genes were already
   selected as highly variable genes and TFs. An absent gene is unmeasured in
   this subset, rather than evidence that its signaling protein is absent.
3. Removes outgoing edges from the observed TFs, matching the treatment inside
   `moon()`, and safely compresses parents with identical signed children.
   Protected observed TF identifiers must survive unchanged. It saves the full
   original-to-scoring identifier mapping and uses the uncompressed PKN if this
   protection check fails.
4. Runs `cosmosR::moon(..., upstream_input=NULL, statistic="ulm", n_layers=4)`.
   At level 1, a candidate regulator's signed effects on observed TFs are tested
   against their activity scores. Subsequent levels repeat this procedure using
   the previous level's inferred scores. The maximum requested depth is four;
   the actual number of realized levels is recorded. No upstream regulator is
   supplied or independently validated.
5. Expands scores to original identifiers and makes a compact explanatory
   network. It keeps nodes with `abs(score)>1.5`, retains edges whose signs agree
   with both endpoint scores, and requires consecutive levels (`k → k−1`) with
   a directed route ending at an observed level-0 TF. An independent audit
   verifies signs, layer descent, and downstream reachability. It repeats the
   display extraction at cutoffs 1, 1.5 and 2 to show threshold sensitivity.

**Why use a custom display extraction?** The current package reducers require
upstream seed nodes; the single-threshold reducer returns an empty network with
`NULL` upstream input. Inventing upstream measurements would change the question.
The downstream-anchored extraction here is an explicit teaching adaptation. It
shows a subset of the scoring relationships and excludes same-level edges,
feedback loops and longer layer-skipping edges from the display. Those display
rules are not an additional causal inference algorithm.

**Why no RNA coherence loop?** This run adds no separate TF-to-transcript
regulatory layer to the signaling PKN. RNA already supplies the TF footprints,
and a target's transcript abundance is not its signaling protein activity.
The NCI60 transcript-sign filtering step is therefore not repeated on these
signaling edges. A future run with explicit transcriptional or protein evidence
can define that additional consistency rule separately. The ledger records a
single fixed-PKN scoring pass and that coherence iteration is not applicable.

## Read the outputs

| Output | Meaning |
| --- | --- |
| `tf_input_coverage.csv` | Every available TF score, threshold decision and PKN coverage; inspect losses before interpretation |
| `moon_scores.csv` | All scores in original IDs, including observed level-0 TFs and every inferred layer |
| `top_regulators.csv` | Inferred nodes ordered by absolute score, with RNA Z score as an annotation and display membership |
| `SIF.csv`, `ATT.csv` | Signed explanatory edges and node scores/levels; import into Cytoscape if desired |
| `network_audit.csv` | Every displayed node has a route to an observed TF; isolated observed TFs are explicitly flagged |
| `threshold_sensitivity.csv` | Changes in network size under alternative display cutoffs |
| `pkn_uncompressed.csv`, `pkn_scoring.csv`, `compression_mapping.csv`, `compression.rds` | Exact calculation topology, identifier mapping and full compression return value |
| `moon_diagnostics.pdf` / `.png` | Inferred score distribution and the largest inferred scores |
| `run_info.json`, `sessionInfo.txt` | Input checksums, mapping, settings, package versions, coverage, actual levels and audit results |

The comparison directory additionally contains `controlled_shared_tf_comparison.csv`,
`native_full_tf_comparison.csv`, `controlled_topology_checks.csv`,
`input_coverage.csv` and `shared_tf_ids.csv`. The controlled manifest states all
fixed inputs and verified graph/level invariants.

## Read the comparison figures

Run `python Heart_pseudobulk/04_compare_moon.py` after script 03. It begins with
the direct level-1 CollecTRI/SCENIC+ score heatmaps and their absolute-score
contrast, using matching rows and color limits for the two source panels. It
then separates levels 1–4 and shows named examples in the selected cell type.
Red in a **support difference** plot means a larger absolute score under
SCENIC+; it does not necessarily mean positive inferred activity. Signed scores
and sign flips are reported alongside the bars. The earlier **TF activity delta**
heatmap instead uses signed `SCENIC+ − CollecTRI`; its red/blue meaning is different.

Illustrative rows have large mean absolute support differences within their
level, with `abs(score)>1.5` in at least one method/profile. The named example
uses up to five groups favoring each source at level 1. This display selection
explains differences; it is not an unbiased estimate of overall performance.

Several different PKN regulators can have the **same signed effects on the
observed TFs**, even when their full-PKN children differ. Script 04 therefore
groups level-1 nodes by this effective measured-input footprint: 377 original
level-1 nodes reduce to 138 signed patterns. For example, HDAC4, HDAC5 and
CCND1_CDK4 all have MEF2A− and MEF2C− as their observed children. Their controlled
ventricular scores are all −3.63 with CollecTRI and −0.62 with SCENIC+. The plot
shows `HDAC4 (+2)` and exports all three names. At later levels, the plots collapse
full-PKN compression groups. Original-ID scores remain available in full.

Inspect `figures/direct_tf_footprint_groups.csv`,
`figures/indistinguishable_regulators.csv` and the selected example table before
interpreting a representative name. The plotted `nTF` counts direct observed
TF children. Current MOON scoring allows a **one-target** footprint; for
example, MED1's direct support here is GATA6 alone. Such predictions should not
be treated as uniquely identified upstream regulators. A larger absolute score
can reflect a coherent but small and prior-dependent footprint.

For example, a positive score for an upstream kinase means its signed downstream
pattern agrees with the inferred TF activities under this PKN. It does **not**
show that the kinase's RNA is high, that phosphorylation was measured, or that
the kinase caused the cell-type difference. Several alternative regulators can
explain overlapping footprints, including compressed groups with the same
children. Thresholds 2 and 1.5 are teaching heuristics, not FDR or universal
significance cutoffs. With no measured upstream anchor, results are especially
useful for generating hypotheses and inspecting prior-dependent explanations.

The small scientific checks exercise the key adaptations:

```bash
Rscript Heart_pseudobulk/test_moon.R
```

They reject disconnected branches, feedback and sign-incoherent display edges,
verify preservation of isolated and zero-valued observed TFs, protect TF IDs
during compression, check complete score expansion without compressed groups,
and verify identical-input reproduction and missing/changed-level exclusion.

Package source for the scoring and reducer behavior:
[cosmosR `moon()` and extraction helpers](https://github.com/saezlab/cosmosR/blob/b2ed195738479c4ad204cace37fd11ae2fd02bb8/R/decoupleRnival.R).
