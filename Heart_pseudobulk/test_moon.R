# Small scientific checks for the RNA-only extraction/compression behavior.
# Run from the repository root: Rscript Heart_pseudobulk/test_moon.R
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
source(file.path(script_dir, "03_moon.R"))
stopifnot(requireNamespace("igraph", quietly = TRUE), requireNamespace("cosmosR", quietly = TRUE))

# A positive chain ending at observed TF must survive. A signed but disconnected
# branch, a feedback edge, and a sign-incoherent edge must not enter the display.
scores <- data.frame(source = c("A", "B", "TF", "X", "Y", "BAD", "ISO"),
                     score = c(3, 2, 4, 3, 3, -3, 3), level = c(2, 1, 0, 2, 1, 1, 0))
pkn <- data.frame(source = c("A", "B", "X", "Y", "BAD"), interaction = 1,
                  target = c("B", "TF", "Y", "X", "TF"))
solution <- extract_downstream_anchored(scores, pkn, cutoff = 1.5)
stopifnot(setequal(paste(solution$SIF$source, solution$SIF$target), c("A B", "B TF")),
          setequal(solution$ATT$nodes, c("A", "B", "TF", "ISO")),
          !solution$ATT$has_explanatory_edge[solution$ATT$nodes == "ISO"])
invisible(audit_solution(solution))

# Protected observed TF must never be renamed into an otherwise identical
# non-input parent group. Two unprotected duplicate parents may be compressed.
pkn <- data.frame(source = c("TF", "P", "Q", "R", "D"), interaction = 1,
                  target = c("C", "C", "D", "D", "E"))
comp <- safe_compression(pkn, c(TF = 3))
stopifnot("TF" %in% network_nodes(comp$pkn),
          comp$mapping$source_scoring[comp$mapping$source_original == "TF"] == "TF",
          comp$used,
          comp$mapping$source_scoring[comp$mapping$source_original == "Q"] ==
            comp$mapping$source_scoring[comp$mapping$source_original == "R"])
scoring_nodes <- unique(comp$mapping$source_scoring)
expanded <- expand_scores(data.frame(source = scoring_nodes, score = seq_along(scoring_nodes), level = 1),
                          comp$mapping)
stopifnot(setequal(expanded$source, network_nodes(pkn)), !anyDuplicated(expanded$source))

# With no duplicate groups the explicit expansion must still preserve all nodes.
no_comp <- safe_compression(pkn[1:2, ], c(TF = 3), enabled = FALSE)
stopifnot(setequal(expand_scores(data.frame(source = network_nodes(pkn[1:2, ]),
                                           score = 1, level = 0), no_comp$mapping)$source,
                   network_nodes(pkn[1:2, ])))
cat("Passed: downstream reachability, signed consecutive layers, isolated TF visibility, protected compression, complete score expansion.\n")
