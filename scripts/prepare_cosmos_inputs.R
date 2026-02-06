library(readr)
library(decoupleR)

# --- Prepare metabolomic inputs -----------------------------------------------

metab_vsn <- as.data.frame(read_csv("data/metabolomic/metabolomic_clean_vsn_common.csv"))
metab_matching <- as.data.frame(read_csv("support/metabolite_matching.csv"))

# Map metabolite names to HMDB IDs
metab_vsn <- merge(metab_matching, metab_vsn, by = "metabolite")
metab_vsn <- metab_vsn[, -1]

# Z-score scaling across cell lines
sds <- apply(metab_vsn[, -1], 1, function(x) sd(x, na.rm = TRUE))
means <- rowMeans(metab_vsn[, -1], na.rm = TRUE)

metab_scaled <- metab_vsn
metab_scaled[, -1] <- (metab_scaled[, -1] - means) / sds

row.names(metab_scaled) <- metab_scaled$hmdb
metab_scaled <- metab_scaled[, -1]

# Convert to per-cell-line named vectors (dropping NAs)
metab_per_cell_line <- apply(metab_scaled, 2, function(x) {
  x[which(!is.na(x))]
}, simplify = FALSE)

# --- Prepare RNA and TF activity inputs ---------------------------------------

rna_fpkm <- as.data.frame(read_csv("data/RNA/RNA_log2_FPKM_cleaned_common.csv"))
row.names(rna_fpkm) <- rna_fpkm$gene
rna_fpkm <- rna_fpkm[, -1]

# Remove genes with too many missing values (> 20 cell lines)
rna_fpkm <- rna_fpkm[rowSums(is.na(rna_fpkm)) <= 20, ]

# Z-score scaling across cell lines
sds <- apply(rna_fpkm, 1, function(x) sd(x, na.rm = TRUE))
means <- rowMeans(rna_fpkm, na.rm = TRUE)

rna_scaled <- (rna_fpkm - means) / sds

# Estimate TF activities from scaled RNA using DoRothEA regulons
# See https://github.com/saezlab/decoupler for details
doro_regulons <- get_dorothea(levels = c("A", "B"))

tf_activities <- apply(rna_scaled, 2, function(x) {
  x <- as.data.frame(x[which(!is.na(x))])
  tfs <- run_wmean(as.matrix(x), network = doro_regulons, times = 1000, minsize = 20)
  tfs <- as.data.frame(tfs)
  tfs <- tfs[which(tfs$statistic == "norm_wmean"), c(2, 4)]
  result <- tfs[, 2]
  names(result) <- tfs[, 1]
  return(result)
})

# Convert RNA to per-cell-line named vectors (dropping NAs)
rna_per_cell_line <- apply(rna_scaled, 2, function(x) {
  x[which(!is.na(x))]
}, simplify = FALSE)

# --- Assemble and save COSMOS inputs per cell line ----------------------------

cosmos_inputs <- sapply(names(tf_activities), function(cell_line) {
  list(
    RNA = rna_per_cell_line[[cell_line]],
    metabolomic = metab_per_cell_line[[cell_line]],
    TF_scores = tf_activities[[cell_line]]
  )
}, USE.NAMES = TRUE, simplify = FALSE)

save(cosmos_inputs, file = "data/cosmos/cosmos_inputs.RData")
