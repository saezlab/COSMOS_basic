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

# Controlled selection keeps every finite shared TF, including a genuine zero,
# and never replaces a missing score or a missing PKN node by zero.
values <- c(TF1 = 0, TF2 = 1, TF3 = -2, TF4 = 0.7, TF5 = 3, MISSING = NA_real_)
controlled_pkn <- data.frame(
  source = c("A", "A", "B", "B", "C", "C", "D", "D", "U", "U", "V", "V"),
  interaction = c(1, -1, -1, 1, 1, 1, -1, -1, 1, -1, -1, 1),
  target = c("TF1", "TF2", "TF2", "TF3", "TF3", "TF4", "TF4", "TF5", "A", "B", "C", "D"))
shared <- names(values)[1:5]
selected <- select_tf_inputs(values, controlled_pkn, cutoff = 2, selected_tf_ids = shared)
stopifnot(identical(names(selected$input), shared), selected$input[["TF1"]] == 0,
          !"MISSING" %in% names(selected$input),
          !selected$coverage$above_activity_cutoff[selected$coverage$source == "TF1"])
stopifnot(inherits(try(select_tf_inputs(values, controlled_pkn, selected_tf_ids = "MISSING"),
                       silent = TRUE), "try-error"))
prepared <- prepare_scoring_network(controlled_pkn, selected$input, list(layers = 2, compress = TRUE))
a <- score_prepared_network(selected$input[prepared$downstream_ids], prepared, n_layers = 2)
b <- score_prepared_network(selected$input[prepared$downstream_ids], prepared, n_layers = 2)
stopifnot(identical(a, b), same_score_levels(a, b),
          a$score[a$source == "TF1" & a$level == 0] == 0)
paired <- compare_score_tables(a, b, "test", "controlled_shared_TFs")
stopifnot(all(paired$delta_score[paired$comparable] == 0),
          all(paired$delta_abs_score[paired$comparable] == 0),
          all(paired$rank_gain_scenicplus[paired$comparable] == 0),
          all(is.na(paired$delta_score[!paired$comparable])))

# Native runs can observe a node as a TF in only one method, score it at
# different layers, or omit it entirely. None receives a fictitious zero delta.
a <- data.frame(source = c("same", "changed", "observed", "only_a"),
                score = c(1, 2, 3, 4), level = c(1, 1, 0, 1))
b <- data.frame(source = c("same", "changed", "observed", "only_b"),
                score = c(-3, 3, 4, 5), level = c(1, 2, 1, 1))
paired <- compare_score_tables(a, b, "test", "native_full_TFs")
stopifnot(sum(paired$comparable) == 1,
          paired$delta_score[paired$source == "same"] == -4,
          paired$delta_abs_score[paired$source == "same"] == 2,
          paired$comparison_status[paired$source == "changed"] == "different_scoring_level",
          paired$comparison_status[paired$source == "observed"] == "observed_TF_in_one_run",
          is.na(paired$score_collectri[paired$source == "only_b"]),
          is.na(paired$score_scenicplus[paired$source == "only_a"]),
          all(is.na(paired$delta_score[!paired$comparable])))
cat("Passed: shared zero-valued TFs retained, identical inputs give identical controlled scores, level-0/changed-level/missing nodes excluded from upstream deltas.\n")
