library(readxl)
library(readr)

# --- Load RNA expression data from raw Excel file ----------------------------

# Source: NCI60 RNA-seq composite expression from CellMiner
# log2(FPKM + 1) values
# https://discover.nci.nih.gov/cellminer/datasets.do
rna_raw <- as.data.frame(
  read_excel("data/RNA/RNA__RNA_seq_composite_expression.xls",
             skip = 10, sheet = 1))

# Drop annotation columns (keep only gene name + expression values)
rna_raw <- rna_raw[, -c(2:6)]

# Set gene names as row names and remove the gene name column
row.names(rna_raw) <- rna_raw[[1]]
rna_raw <- rna_raw[, -1]

# Strip prefix from column headers (e.g. "BR:MCF7" -> "MCF7")
names(rna_raw) <- gsub(".*[:]", "", names(rna_raw))

# --- Filter lowly expressed genes ---------------------------------------------

# Set values below threshold to NA (log2 FPKM < 1 considered not expressed)
rna_raw[rna_raw < 1] <- NA

# Keep genes expressed in at least 11 cell lines
rna_filtered <- rna_raw[rowSums(is.na(rna_raw)) < 45, ]

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
