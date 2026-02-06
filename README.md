# COSMOS tutorial with the NCI-60 cell line panel

This repository provides a step-by-step tutorial for running [COSMOS](https://github.com/saezlab/cosmosR) (Causal Oriented Search of Multi-Omic Space) on multi-omics data from the NCI-60 cancer cell line collection. It demonstrates how transcription factor activities, metabolomic measurements, and RNA expression data can be connected through a prior knowledge network (PKN) using the **MOON** (Meta-fOOtprint aNalysis) scoring algorithm.

**Publication:** Dugourd *et al.*, "Modeling causal signal propagation in multi-omic factor space with COSMOS", *bioRxiv* (2024). [doi:10.1101/2024.07.15.603538](https://doi.org/10.1101/2024.07.15.603538)

## What does the main script do?

The main analysis script (`scripts/net_compr_MOON.R`) takes pre-processed inputs for a single NCI-60 cell line and builds a mechanistic network connecting upstream transcription factor (TF) activities to downstream metabolomic measurements through a signed, directed prior knowledge network. The pipeline has five stages:

1. **PKN filtering.** The generic COSMOS prior knowledge network (~100k interactions spanning signaling, gene regulation, and metabolism) is pruned by removing genes not expressed in the cell line of interest. Then, nodes that cannot be reached from the upstream TF inputs (*controllable*) or that have no path to the downstream metabolomic inputs (*observable*) within a given number of steps are removed. This focuses the network on the relevant biology.

2. **Network compression.** Parent nodes that regulate the exact same set of children are merged into a single representative node. This reduces the network size (and MOON runtime) without discarding any information; the mapping is stored so that results can be decompressed afterwards.

3. **MOON scoring.** MOON iteratively scores every node of the compressed network using decoupleR's Univariate Linear Model (ULM). Scoring starts at the downstream layer (metabolomic measurements) and propagates upstream, layer by layer, toward the TF activity scores. At each layer, a node's score reflects how consistently its direct downstream targets are up- or down-regulated (accounting for interaction signs). After each full scoring pass, TF-to-target edges that are incoherent with the observed RNA expression are removed, and the scoring is repeated until no more edges are pruned (convergence).

4. **Decompression.** The compressed node groups are expanded back to individual nodes, recovering the full-resolution scores.

5. **Solution network extraction.** A score threshold is applied to select the most significant nodes, and the corresponding subnetwork is extracted from the original (uncompressed) PKN. Metabolite HMDB identifiers are translated to human-readable names. The output is a **SIF** (edge list) and **ATT** (node attribute table) pair that can be visualized in Cytoscape or similar tools.

## Key parameters to explore

The script exposes several parameters that control the scope, resolution, and stringency of the analysis. Experimenting with them is a good way to build intuition for what MOON captures and where its limits are.

### Network depth (`n_steps`)

Controls how many steps away from the input nodes the PKN is allowed to extend during the controllability/observability pruning. A smaller value (e.g. 4) yields a compact network focused on direct regulatory mechanisms, while a larger value (e.g. 8) allows more indirect paths and longer signaling cascades, at the cost of a larger search space and potentially more noise. The default of **6** is a reasonable trade-off for the NCI-60 data. This parameter also sets the number of MOON scoring layers (`n_layers`), meaning it determines how far upstream the iterative footprint propagation can reach.

### TF activity significance threshold

The line `tf_input <- tf_input[abs(tf_input) > 2]` filters the upstream TF input to only keep TFs whose estimated activity score exceeds 2 in absolute value (roughly a t-value > 2 significance level). Lowering this threshold includes more TFs as upstream seeds (broader but noisier), while raising it focuses on the most confidently active TFs. This directly affects the set of upstream nodes MOON tries to connect to.

### Network compression

Compression (`compress_same_children`) is not just a performance optimization. Because MOON scores nodes through iterative ULM runs, redundant parent nodes (those sharing the exact same children) would each receive their own score and then jointly contribute to subsequent upstream ULM rounds as if they were independent observations, when in reality they carry the same information. Merging them into a single representative node prevents this inflated influence on upstream scores. Comparing results with and without compression (by passing the uncompressed network to `moon()` directly) can reveal whether redundant paths are biasing certain scores.

### Dual score threshold for solution extraction

The function `reduce_solution_network_double_thresh()` uses two thresholds to extract the final subnetwork:

- **Primary threshold** (e.g. 1.7): selects the high-confidence core nodes whose absolute MOON score exceeds this value. These serve as seeds for the solution network.
- **Secondary threshold** (e.g. 1.0): allows neighbors of the core nodes to be included if their score exceeds this softer cutoff, filling in the connecting paths.

This two-tier approach gives richer, more connected networks than a single hard cutoff, while still anchoring the result around the most significant nodes. Raising the primary threshold produces smaller, higher-confidence networks; lowering the secondary threshold includes more of the surrounding context. The single-threshold alternative (`reduce_solution_network()` with a single `cutoff` parameter) is also available for simpler exploration.

### Metabolite compartment assignment

`prepare_metab_inputs(metab_input, c("c", "m"))` duplicates each metabolite measurement across the specified cellular compartments (here cytosol and mitochondria). Adding or removing compartments (e.g. `"r"` for endoplasmic reticulum, `"e"` for extracellular) changes which metabolic reactions in the PKN can be reached, and therefore which metabolic mechanisms appear in the solution network.

## Installation

```r
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("saezlab/cosmosR")
```

Other dependencies (`reshape2`, `readr`, `igraph`, `decoupleR`) will be pulled in automatically.

## Repository structure

```
scripts/
  net_compr_MOON.R            # Main MOON analysis pipeline (start here)
  prepare_cosmos_inputs.R     # Preprocessing: RNA, metabolomics, TF activities -> per-cell-line inputs
  RNA.R                       # RNA-seq preprocessing from raw Excel
  metabolomic.R               # Metabolomics preprocessing and VSN normalization
data/
  cosmos/                     # Pre-processed COSMOS inputs (RData, per cell line)
  metabolomic/                # Raw and processed metabolomics data
  RNA/                        # Raw Excel and processed RNA expression data
support/                      # TF regulons, pathway GMTs, metabolite matching tables
results/                      # Output SIF/ATT network files and MOON scores
```

## Running the tutorial

1. Install `cosmosR` (see above).
2. Open the project in RStudio (`NCI60.Rproj`).
3. Run `scripts/net_compr_MOON.R`. By default it analyzes the **786-0** cell line; change the `cell_line` variable to any of the 56 available NCI-60 cell lines.
4. Outputs are written to `results/`: a SIF edge list and an ATT node attribute table, ready for visualization.

To regenerate the input data from scratch, run the preprocessing scripts in order: `RNA.R` and `metabolomic.R` first, then `prepare_cosmos_inputs.R`.
