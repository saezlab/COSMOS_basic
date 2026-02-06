# CLAUDE.md

## Project Overview

This is a **COSMOS (Causal Oriented Search of Multi-Omic Space)** tutorial/analysis project using the NCI60 cancer cell line panel. It demonstrates multi-omics data integration by connecting transcription factor activities (upstream) to metabolomic measurements (downstream) through a prior knowledge network (PKN), using the `cosmosR` R package and the MOON (Multi-Omics Optimized Network) scoring algorithm.

## Language & Framework

- **Language**: R
- **Core package**: `cosmosR` (installed from GitHub: `saezlab/cosmosR`)
- **Key dependencies**: `decoupleR`, `reshape2`, `readr`, `ggplot2`, `GSEABase`, `igraph`
- **Project type**: RStudio project (`NCI60.Rproj`)
- **Main analysis**: R Markdown (`net_compr_MOON.Rmd`)

## Directory Structure

```
├── data/                  # Input data
│   ├── cosmos/            # Pre-processed COSMOS inputs (RData, per cell line)
│   ├── metabolomic/       # Raw metabolomic data
│   └── RNA/               # RNA expression data
├── scripts/               # R scripts
│   ├── prepare_cosmos_inputs.R   # Data preprocessing
│   ├── net_compr_MOON.R          # Main MOON analysis (script version)
│   ├── net_compr_carnival.R      # CARNIVAL-based analysis (compressed)
│   ├── net_no_compr_carnival.R   # CARNIVAL-based analysis (uncompressed)
│   ├── RNA.R / metabolomic.R     # Data processing scripts
│   ├── support_compression.R     # Network compression utilities
│   ├── support_decompression.R   # Network decompression utilities
│   └── support_decoupleRnival.R  # decoupleR/CARNIVAL support functions
├── support/               # Support files (regulons, pathway GMTs, helper functions)
├── results/               # Output results (SIF/ATT network files, MOON scores)
│   └── moon/              # MOON scoring results per cell line
├── EBI_practical/         # Presentation materials
├── cbc/                   # CBC tools
├── net_compr_MOON.Rmd     # Main R Markdown tutorial
└── net_compr_MOON_files/  # Rendered figures from Rmd
```

## Key Concepts & Workflow

1. **Load data**: TF activity scores, metabolomics, RNA expression for a given NCI60 cell line
2. **Prepare inputs**: Assign metabolite compartments, filter significant TF inputs
3. **Filter PKN**: Remove unexpressed genes, prune network to reachable/observable nodes within `n_steps` (default: 6)
4. **Compress network**: `compress_same_children()` to reduce redundant parent nodes
5. **Run MOON**: Iterative scoring with `moon()` + `filter_incohrent_TF_target()` until convergence
6. **Decompress & extract**: Decompress results, apply cutoff, extract solution subnetwork (SIF + ATT format)
7. **Pathway enrichment**: Optional Reactome pathway scoring via `decoupleR::run_ulm()`

## Output Formats

- **SIF**: Source-Interaction-Target edge list (columns: `source`, `target`, `sign`)
- **ATT**: Node attribute table (columns: `Nodes`, `score`, plus activity values)
- Results are written as CSVs to `results/`

## Conventions

- Cell line names follow NCI60 naming (e.g., `"786-0"`, `"A549/ATCC"`, `"MCF7"`)
- Metabolite IDs use HMDB format with compartment suffix (e.g., `Metab__HMDB0011747_c`)
- Internal `cosmosR` functions are accessed via `cosmosR:::` (triple colon)
- Network operations use data frames with `source`/`target`/`interaction` columns
- DoRothEA regulons are loaded from `support/dorothea_reg.RData`
