# Regenerate the bundled prior from the existing COSMOS teaching resource.
# Run from any directory: Rscript /path/to/GRN_comparison/export_collectri.R
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(script) != 1L) stop("Run this file with Rscript.")
here <- dirname(normalizePath(script))
input <- file.path(here, "..", "support", "collectri_regulon_R.RData")
env <- new.env()
load(input, envir = env)
stopifnot(all(c("source", "target", "mor") %in% names(env$regulons)))
prior <- unique(env$regulons[, c("source", "target", "mor")])
names(prior)[3] <- "score"
prior <- prior[order(prior$source, prior$target, prior$score), ]
dir.create(file.path(here, "data"), showWarnings = FALSE)
output <- file.path(here, "data", "collectri.tsv")
write.table(prior, output,
            sep = "\t", row.names = FALSE, quote = FALSE)
provenance <- c(
  '{',
  '  "description": "Text export of the local CollecTRI snapshot; no live database query.",',
  '  "generated_with": "Rscript GRN_comparison/export_collectri.R",',
  '  "source": "support/collectri_regulon_R.RData",',
  sprintf('  "source_md5": "%s",', unname(tools::md5sum(input))),
  sprintf('  "export_md5": "%s",', unname(tools::md5sum(output))),
  sprintf('  "rows": %d,', nrow(prior)),
  '  "source_release": "Not recorded in the existing RData; do not assume latest release.",',
  '  "score": "Original mor annotation; not a confidence rank."',
  '}'
)
writeLines(provenance, file.path(here, "data", "provenance.json"))
message("Exported ", nrow(prior), " prior rows; no network download required.")
