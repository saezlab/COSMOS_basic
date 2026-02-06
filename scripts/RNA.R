library(readr)

# --- Load RNA expression data ------------------------------------------------

rna_raw <- as.data.frame(
  read_delim("data/RNA/RNA_log2_FPKM.csv",
             delim = ";", escape_double = FALSE, trim_ws = TRUE))

# Drop annotation columns (keep only gene name + expression values)
rna_raw <- rna_raw[, -c(2:6)]

# Strip prefix from column headers (e.g. "prefix:cellline" -> "cellline")
names(rna_raw) <- gsub(".*[:]", "", names(rna_raw))

# Remove rows with any missing values
rna_raw <- rna_raw[complete.cases(rna_raw), ]

# Set gene names as row names and remove the gene name column
row.names(rna_raw) <- rna_raw$`Gene name d`
rna_raw <- rna_raw[, -1]

# Fix decimal separator (European locale uses comma instead of dot)
rna_raw <- as.data.frame(apply(rna_raw, 2, function(x) gsub(",", ".", x)))

# --- Filter lowly expressed genes ---------------------------------------------

hist(as.numeric(unlist(rna_raw)), breaks = 1000)

# Set values below threshold to NA (log2 FPKM < 1 considered not expressed)
rna_raw[rna_raw < 1] <- NA
hist(as.numeric(unlist(rna_raw)), breaks = 1000)

# Keep genes expressed in at least ~20% of cell lines (56 - 45 = 11)
n_cell_lines <- ncol(rna_raw)
max_missing <- n_cell_lines - 11
rna_filtered <- rna_raw[rowSums(is.na(rna_raw)) < max_missing, ]
hist(as.numeric(unlist(rna_filtered)), breaks = 1000)

# --- Align cell lines with metabolomics data ----------------------------------

metabolomic_vsn <- as.data.frame(read_csv("data/metabolomic/metabolomic_clean_vsn.csv"))

common_cell_lines <- intersect(names(metabolomic_vsn[, -1]), names(rna_filtered))

metabolomic_vsn <- metabolomic_vsn[, c("metabolite", common_cell_lines)]

rna_filtered <- rna_filtered[, common_cell_lines]
rna_filtered$gene <- row.names(rna_filtered)

# Reorder to put gene name column first
rna_filtered <- rna_filtered[, c("gene", common_cell_lines)]

# --- Write outputs ------------------------------------------------------------

write_csv(rna_filtered, file = "data/RNA/RNA_log2_FPKM_cleaned_common.csv")
write_csv(metabolomic_vsn, file = "data/metabolomic/metabolomic_clean_vsn_common.csv")
