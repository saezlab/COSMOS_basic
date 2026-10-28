# Optional regeneration of the bundled numerical reference, using base R only.
# This exports archived scores; it does not run TF activity inference.
# Run from any directory: Rscript /path/to/TF_activity/export_r_reference.R

args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(script) != 1L) stop("Run this file with Rscript.")
here <- dirname(normalizePath(script))
root <- normalizePath(file.path(here, ".."))
source_relative <- "data/cosmos/cosmos_inputs.RData"
source_file <- file.path(root, source_relative)

env <- new.env(parent = emptyenv())
load(source_file, envir = env)
if (!is.list(env$cosmos_inputs) || is.null(names(env$cosmos_inputs))) {
  stop("The source must contain a named list called cosmos_inputs.")
}

records <- lapply(names(env$cosmos_inputs), function(cell_line) {
  scores <- env$cosmos_inputs[[cell_line]]$TF_scores
  if (!is.numeric(scores) || is.null(names(scores)) ||
      anyNA(scores) || any(!is.finite(scores)) ||
      anyDuplicated(names(scores))) {
    stop("Invalid archived TF_scores vector for ", cell_line)
  }
  data.frame(source = names(scores), cell_line = cell_line,
             score = unname(scores), stringsAsFactors = FALSE)
})
reference <- do.call(rbind, records)
row.names(reference) <- NULL
if (anyDuplicated(reference[c("source", "cell_line")])) {
  stop("The reference contains duplicate TF / cell-line pairs.")
}

output_dir <- file.path(here, "data")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(output_dir, "r_reference_tf_scores.tsv.gz")
header <- "source\tcell_line\tscore"
# Seventeen significant digits preserve the stored double precision values.
lines <- paste(reference$source, reference$cell_line,
               sprintf("%.17g", reference$score), sep = "\t")
connection <- gzfile(output_file, open = "wt", compression = 9)
tryCatch(writeLines(c(header, lines), connection, useBytes = TRUE),
         finally = close(connection))

json_string <- function(x) encodeString(x, quote = '"')
description <- paste(
  "Text export of the archived TF_scores vectors from cosmos_inputs.RData.",
  "Used only to validate Python activity estimates; never used as an inference input.",
  "The R package versions that generated the original archive were not recorded."
)
provenance <- c(
  "{",
  paste0('  "description": ', json_string(description), ","),
  '  "generated_with": "Rscript TF_activity/export_r_reference.R",',
  paste0('  "source": ', json_string(source_relative), ","),
  paste0('  "source_md5": ', json_string(unname(tools::md5sum(source_file))), ","),
  '  "export": "TF_activity/data/r_reference_tf_scores.tsv.gz",',
  paste0('  "export_md5": ', json_string(unname(tools::md5sum(output_file))), ","),
  sprintf('  "n_records": %d,', nrow(reference)),
  sprintf('  "n_cell_lines": %d,', length(unique(reference$cell_line))),
  '  "score": "Archived decoupleR ULM t-statistic; no p-value filtering.",',
  '  "numeric_precision": "17 significant decimal digits preserve stored doubles."',
  "}"
)
writeLines(provenance, file.path(output_dir, "r_reference_provenance.json"),
           useBytes = TRUE)
message("Exported ", nrow(reference), " archived TF activity scores from ",
        length(unique(reference$cell_line)), " cell lines.")
