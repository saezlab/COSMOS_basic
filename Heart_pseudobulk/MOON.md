# RNA-only MOON: upstream hypotheses from TF activities

`03_moon.R` asks **which upstream signaling activities are compatible with a
cell type's inferred TF activity pattern?** It adapts the NCI60 MOON practical
to a case with no measured upstream layer and no metabolomics.

The inputs are the unthresholded `tf_scores.csv` and `rna_scaled.csv` written
by `01_tf_activity.py` in `results/tf_activity/`.

| Quantity | Role in this run |
| --- | --- |
| Signed CollecTRI ULM TF scores | Downstream protein-activity readouts; retain finite scores with `abs(score)>2` by default |
| Upstream measurements | Unspecified: `upstream_input=NULL` |
| Gene-wise RNA Z scores | Node annotations only; do not substitute RNA abundance for protein activity |
| SCENIC+ network | Evaluated separately in script 02; not used to constrain this MOON run |
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

This analyzes **Ventricular Cardiomyocyte**. To choose a different column in the
TF score matrix or run all six cell types:

```bash
Rscript Heart_pseudobulk/03_moon.R --cell-type "Myeloid"
Rscript Heart_pseudobulk/03_moon.R --all
```

`--help` lists options for input/output directories, layer count, thresholds,
and disabling compression. Default output directories are
`Heart_pseudobulk/results/moon/<cell_type>/`. Re-running a cell type updates that
cell type's outputs in this new practical; the NCI60 results are untouched.

## What the script does

1. Loads the bundled `cosmosR::meta_network`, preserving a raw copy and package
   provenance. Keeps its protein/complex signaling layer, excluding metabolite
   nodes (`Metab__...`) and reaction-specific enzyme copies (`Gene<digits>__...`).
   It cleans self-loops, duplicate pairs and conflicting signs using the package
   cleaner, with before/after counts in `pkn_preparation.csv`.
2. Retains sufficiently large TF activities that occur in the signaling PKN,
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
verify preservation of an isolated observed TF, protect TF IDs during
compression and check complete score expansion even without compressed groups.

Package source for the scoring and reducer behavior:
[cosmosR `moon()` and extraction helpers](https://github.com/saezlab/cosmosR/blob/b2ed195738479c4ad204cace37fd11ae2fd02bb8/R/decoupleRnival.R).
