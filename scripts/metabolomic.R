library(readr)
library(reshape2)
library(pheatmap)
library(vsn)
library(metaboliteIDmapping)

# --- Load and filter raw metabolomics data ------------------------------------

metabs_raw <- as.data.frame(read_csv("data/metabolomic/WEB_DATA_METABOLON.TXT"))

# Remove ambiguous identifications
metabs_raw <- metabs_raw[!grepl("Isobar", metabs_raw$TITLE), ]
metabs_raw <- metabs_raw[!grepl("^X[-]", metabs_raw$TITLE), ]
metabs_raw <- metabs_raw[!grepl("^Possible", metabs_raw$TITLE), ]

# Pivot to metabolite x cell line matrix (average duplicates)
metabs_df <- dcast(metabs_raw, formula = TITLE ~ cellname, value.var = "VALUE",
                   fun.aggregate = mean)
row.names(metabs_df) <- metabs_df$TITLE
metabs_df <- metabs_df[, -1]

# --- QC and outlier removal ---------------------------------------------------

hist(as.numeric(unlist(log2(metabs_df))), breaks = 1000)
pheatmap(log2(metabs_df), show_colnames = FALSE, show_rownames = FALSE)

# Remove extreme values (|log2| > 5) before normalization
metabs_df[abs(log2(metabs_df)) > 5] <- NA

# --- VSN normalization --------------------------------------------------------

fit <- vsnMatrix(as.matrix(metabs_df))
meanSdPlot(fit)
metabs_vsn <- as.data.frame(vsn::predict(fit, as.matrix(metabs_df)))

pheatmap(metabs_vsn[complete.cases(metabs_vsn), ],
         show_colnames = TRUE, show_rownames = FALSE, cluster_rows = FALSE)

# --- Write VSN-normalized output ----------------------------------------------

metabs_output <- metabs_vsn
metabs_output$metabolite <- row.names(metabs_output)
metabs_output$metabolite <- gsub(",", "", metabs_output$metabolite)

# Reorder to put metabolite name column first
metabs_output <- metabs_output[, c("metabolite", setdiff(names(metabs_output), "metabolite"))]

write_csv(metabs_output, file = "data/metabolomic/metabolomic_clean_vsn.csv")

# --- Optional: z-score scaling for exploration --------------------------------

sds <- apply(metabs_vsn, 1, function(x) sd(x, na.rm = TRUE))
means <- rowMeans(metabs_vsn, na.rm = TRUE)

metabs_scaled <- (metabs_vsn - means) / sds
pheatmap(metabs_scaled[complete.cases(metabs_scaled), ],
         show_colnames = TRUE, show_rownames = FALSE, cluster_rows = FALSE)

hist(as.numeric(unlist(metabs_vsn)), breaks = 1000)
