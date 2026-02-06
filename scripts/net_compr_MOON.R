# We advise to instal from github to get the latest version of the tool.
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("saezlab/cosmosR")


library(cosmosR)
library(reshape2)
library(readr)
library(igraph)

# --- Load prior knowledge network (PKN) --------------------------------------

data("meta_network")

meta_network <- meta_network_cleanup(meta_network)

# --- Load pre-processed cell line data ----------------------------------------

load("data/cosmos/cosmos_inputs.RData")

names(cosmos_inputs)

cell_line <- "786-0"

tf_input    <- cosmos_inputs[[cell_line]]$TF_scores
metab_input <- cosmos_inputs[[cell_line]]$metabolomic
rna_input   <- cosmos_inputs[[cell_line]]$RNA

# --- Prepare MOON inputs ------------------------------------------------------

# Assign metabolites to cellular compartments (c = cytosol, m = mitochondria)
metab_input <- prepare_metab_inputs(metab_input, c("c", "m"))

# Keep only TFs with significant activity (|t-value| > 2)
tf_input <- tf_input[abs(tf_input) > 2]

# Remove genes not expressed in this cell line from the PKN
meta_network <- cosmosR:::filter_pkn_expressed_genes_fast(
  names(rna_input), meta_pkn = meta_network)

# --- Prune PKN to reachable subnetwork ----------------------------------------

# Iteratively remove:
#   - input nodes not present in the PKN
#   - PKN nodes not reachable from upstream inputs within n_steps (controllable)
#   - PKN nodes that cannot reach downstream inputs within n_steps (observable)
# Repeat until no more nodes are removed.

n_steps <- 6

n_inputs_before <- c(length(tf_input), length(metab_input))
n_inputs_after  <- 0

while (sum(n_inputs_before == n_inputs_after) != 2) {
  n_inputs_before <- c(length(tf_input), length(metab_input))
  tf_input    <- cosmosR:::filter_input_nodes_not_in_pkn(tf_input, meta_network)
  meta_network <- cosmosR:::keep_controllable_neighbours(meta_network, n_steps, names(tf_input))
  metab_input <- cosmosR:::filter_input_nodes_not_in_pkn(metab_input, meta_network)
  meta_network <- cosmosR:::keep_observable_neighbours(meta_network, n_steps, names(metab_input))
  n_inputs_after <- c(length(tf_input), length(metab_input))
}

# --- Compress redundant network paths -----------------------------------------

# Merge parent nodes that share identical sets of children, reducing the
# network size without losing information. This speeds up MOON scoring.

compressed <- compress_same_children(
  meta_network, sig_input = tf_input, metab_input = metab_input)

meta_network_compressed <- compressed$compressed_network
node_signatures         <- compressed$node_signatures
duplicated_parents      <- compressed$duplicated_signatures

meta_network_compressed <- meta_network_cleanup(meta_network_compressed)

# --- Run MOON scoring with TF-target coherence filtering ----------------------

# MOON (Meta-fOOtprint aNalysis) iteratively scores network nodes layer by
# layer using ULM (Univariate Linear Model), starting from the downstream
# metabolomic measurements and propagating upstream toward TF activity scores.
#
# After each scoring round, TF-target interactions incoherent with RNA
# expression data are removed from the network. The loop repeats until
# convergence (no more edges removed).

load("support/collectri_regulon_R.RData")

moon_network <- meta_network_compressed

n_edges_before <- 1
n_edges_after  <- 0
i <- 1
while (n_edges_before != n_edges_after & i < 10) {
  n_edges_before <- nrow(moon_network)
  moon_res <- moon(upstream_input = tf_input,
                   downstream_input = metab_input,
                   meta_network = moon_network,
                   n_layers = n_steps,
                   statistic = "ulm")

  moon_network <- filter_incohrent_TF_target(
    moon_res, regulons, moon_network, rna_input)
  n_edges_after <- nrow(moon_network)
  i <- i + 1
}

if (i < 10) {
  print(paste0("Converged after ", i - 1, " iterations"))
} else {
  print(paste0("Interrupted after ", i, " iterations. Convergence uncertain."))
}

# --- Decompress and threshold MOON results ------------------------------------

# Save full (compressed) MOON scores for reference
write_csv(moon_res,
          file = paste0("results/moon/", cell_line, "_ATT_moon_full.csv"))

# Expand compressed node groups back to individual nodes
moon_res <- decompress_moon_result(moon_res, compressed, moon_network)

# Inspect score distribution to inform cutoff choice
plot(density(moon_res$score))
abline(v = 1)
abline(v = -1)

# Reorder columns to match expected format: source, score, level
moon_res <- moon_res[, c(4, 2, 3)]
names(moon_res)[1] <- "source"

# --- Extract solution subnetwork ----------------------------------------------

# Extract solution subnetwork using a double threshold: nodes with scores above
# the primary threshold (1.5) are selected first, then the subnetwork is
# expanded to include neighboring nodes above the secondary threshold (0.5).

solution_network <- reduce_solution_network_double_thresh(
  decoupleRnival_res = moon_res,
  meta_network = meta_network,
  primary_thresh = 1.7,
  secondary_thresh = 1,
  upstream_input = tf_input,
  RNA_input = rna_input)

SIF <- solution_network$SIF
names(SIF)[3] <- "sign"
ATT <- solution_network$ATT

# --- Translate metabolite IDs to human-readable names -------------------------

data("HMDB_mapper_vec")

translated_res <- translate_res(SIF, ATT, HMDB_mapper_vec)

SIF <- translated_res[[1]]
ATT <- translated_res[[2]]

# --- Write final results ------------------------------------------------------

write_csv(SIF, file = paste0("results/", cell_line, "_dec_compressed_SIF.csv"))
write_csv(ATT, file = paste0("results/", cell_line, "_dec_compressed_ATT.csv"))
