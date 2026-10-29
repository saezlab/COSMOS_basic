#!/usr/bin/env Rscript
# RNA-only MOON: infer upstream activity hypotheses from downstream TF footprints.
# Run 01_tf_activity.py first. See MOON.md for the mapping and interpretation.
# No installation is performed by this script and no CARNIVAL solver is needed.

network_nodes <- function(pkn) unique(c(pkn$source, pkn$target))

ordered_pkn <- function(pkn) {
  pkn <- as.data.frame(pkn[, c("source", "interaction", "target")])
  pkn <- pkn[order(pkn$source, pkn$target, pkn$interaction), ]
  rownames(pkn) <- NULL
  pkn
}

validate_pkn <- function(pkn) {
  stopifnot(all(c("source", "interaction", "target") %in% names(pkn)),
            !anyNA(pkn), is.numeric(pkn$interaction),
            all(pkn$interaction %in% c(-1, 1)),
            all(nzchar(pkn$source)), all(nzchar(pkn$target)),
            !any(pkn$source == pkn$target),
            !anyDuplicated(paste(pkn$source, pkn$target, sep = "\t")))
  invisible(TRUE)
}

safe_compression <- function(pkn, downstream_input, enabled = TRUE) {
  # The package compressor uses column positions: source, target, interaction.
  # Put protected parents last; this avoids a current package ordering pitfall.
  protected <- intersect(names(downstream_input), network_nodes(pkn))
  original_nodes <- network_nodes(pkn)
  result <- list(pkn = pkn, mapping = data.frame(
    source_original = original_nodes, source_scoring = original_nodes),
    used = FALSE, reason = "disabled", package_result = NULL)
  if (!enabled) return(result)
  ordered <- pkn[order(pkn$source %in% protected, pkn$source,
                       pkn$target, pkn$interaction), ]
  protection <- setNames(rep(0, length(protected)), protected)
  compressed <- cosmosR::compress_same_children(
    ordered[, c("source", "target", "interaction")],
    sig_input = protection, metab_input = downstream_input)
  renamed <- compressed$duplicated_signatures
  scoring <- ordered_pkn(cosmosR::meta_network_cleanup(
    compressed$compressed_network[, c("source", "interaction", "target")]))
  safe <- !any(protected %in% names(renamed)) &&
    all(protected %in% network_nodes(scoring))
  if (!safe) {
    result$reason <- "protected TF identity check failed; used uncompressed PKN"
    return(result)
  }
  map <- setNames(original_nodes, original_nodes)
  map[names(renamed)] <- as.character(renamed)
  projected <- pkn
  projected$source <- unname(map[projected$source])
  projected$target <- unname(map[projected$target])
  projected <- ordered_pkn(cosmosR::meta_network_cleanup(projected))
  if (!identical(projected, scoring)) {
    result$reason <- "compressed edge projection check failed; used uncompressed PKN"
    return(result)
  }
  result$pkn <- scoring
  result$mapping$source_scoring <- unname(map[original_nodes])
  result$used <- length(renamed) > 0
  result$reason <- if (result$used) "safe same-children compression" else "no duplicate parent signatures"
  result$package_result <- compressed
  validate_pkn(result$pkn)
  result
}

expand_scores <- function(scores, mapping) {
  # Explicit mapping includes every original node, including sink TFs. This
  # also handles zero compressed groups without relying on positional indexing.
  expanded <- merge(mapping, scores, by.x = "source_scoring", by.y = "source")
  names(expanded)[names(expanded) == "source_original"] <- "source"
  expanded <- expanded[, c("source", "score", "level", "source_scoring")]
  stopifnot(!anyDuplicated(expanded$source))
  expanded[order(expanded$level, -abs(expanded$score), expanded$source), ]
}

extract_downstream_anchored <- function(scores, pkn, cutoff = 1.5, rna_input = NULL) {
  # The package reducers require upstream anchors. This RNA-only case has none.
  # Keep score/sign-consistent edges between successive scored layers, then
  # require a directed route to one of the observed level-0 TFs. Level descent
  # excludes feedback cycles from this small explanatory display, not scoring.
  selected <- scores[is.finite(scores$score) & abs(scores$score) > cutoff, ]
  score_map <- setNames(selected$score, selected$source)
  level_map <- setNames(selected$level, selected$source)
  anchors <- selected$source[selected$level == 0]
  sif <- pkn[pkn$source %in% selected$source & pkn$target %in% selected$source, ]
  if (nrow(sif)) {
    sif <- sif[level_map[sif$source] == level_map[sif$target] + 1 &
                 sign(score_map[sif$source] * score_map[sif$target]) == sif$interaction, ]
  }
  if (nrow(sif) && length(anchors)) {
    graph <- igraph::graph_from_data_frame(
      sif[, c("source", "target")], directed = TRUE)
    anchors_present <- intersect(anchors, igraph::V(graph)$name)
    reachable <- unique(unlist(lapply(anchors_present, function(node)
      igraph::as_ids(igraph::subcomponent(graph, node, mode = "in")))))
    sif <- sif[sif$source %in% reachable & sif$target %in% reachable, ]
  } else {
    sif <- sif[FALSE, ]
  }
  sif$consistency <- rep(TRUE, nrow(sif))
  sif <- sif[, c("source", "target", "interaction", "consistency")]
  # Retain selected measured TFs even if no surviving explanatory edge reaches
  # them. Their isolation is explicit, rather than silently losing observations.
  final_nodes <- unique(c(network_nodes(sif), anchors))
  att <- selected[selected$source %in% final_nodes, c("source", "score", "level")]
  names(att)[1] <- "nodes"
  att$RNA_input <- if (is.null(rna_input)) NA_real_ else unname(rna_input[att$nodes])
  att$role <- ifelse(att$level == 0, "observed TF footprint", "inferred upstream activity")
  att$has_explanatory_edge <- att$nodes %in% network_nodes(sif)
  list(SIF = sif, ATT = att)
}

audit_solution <- function(solution) {
  sif <- solution$SIF
  att <- solution$ATT
  score_map <- setNames(att$score, att$nodes)
  level_map <- setNames(att$level, att$nodes)
  stopifnot(!anyDuplicated(att$nodes), !anyNA(att[, c("nodes", "score", "level")]),
            all(sif$source %in% att$nodes), all(sif$target %in% att$nodes),
            all(sign(score_map[sif$source] * score_map[sif$target]) == sif$interaction),
            all(level_map[sif$source] == level_map[sif$target] + 1))
  anchors <- att$nodes[att$level == 0]
  path_ok <- setNames(rep(FALSE, nrow(att)), att$nodes)
  path_ok[anchors] <- TRUE
  if (nrow(sif)) {
    graph <- igraph::graph_from_data_frame(sif[, c("source", "target")], directed = TRUE)
    for (node in intersect(anchors, igraph::V(graph)$name)) {
      path_ok[igraph::as_ids(igraph::subcomponent(graph, node, mode = "in"))] <- TRUE
    }
    stopifnot(igraph::is_dag(graph))
    components <- igraph::components(graph, mode = "weak")
    stopifnot(all(vapply(split(names(components$membership), components$membership),
                        function(nodes) any(nodes %in% anchors), logical(1))))
  }
  stopifnot(all(path_ok))
  data.frame(nodes = att$nodes, reaches_observed_TF = unname(path_ok[att$nodes]),
             level = att$level, has_explanatory_edge = att$has_explanatory_edge)
}

read_score_matrix <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path, ". Run 01_tf_activity.py first.")
  mat <- as.matrix(utils::read.csv(path, row.names = 1, check.names = FALSE))
  storage.mode(mat) <- "double"
  stopifnot(!anyDuplicated(rownames(mat)), !anyDuplicated(colnames(mat)),
            !anyNA(rownames(mat)), !anyNA(colnames(mat)))
  mat
}

plot_diagnostics <- function(scores, cell_type, path, cutoff) {
  draw <- function() {
    par(mfrow = c(1, 2), mar = c(4.2, 5.2, 3, 1))
    inferred <- scores[scores$level > 0 & is.finite(scores$score), ]
    hist(inferred$score, breaks = 40, col = "grey75", border = "white",
         main = "All inferred MOON scores", xlab = "ULM score (activity alignment)")
    abline(v = c(-cutoff, cutoff), lty = 2, col = "grey40")
    best <- head(inferred[order(-abs(inferred$score)), ], 15)
    best <- best[order(best$score), ]
    par(mar = c(4.2, 7.5, 3, 1))
    barplot(best$score, names.arg = best$source, horiz = TRUE, las = 1,
            col = ifelse(best$score > 0, "#1f78b4", "#e66101"), border = NA,
            main = "Largest inferred scores", xlab = "ULM score")
    abline(v = 0, col = "grey40")
    mtext(cell_type, outer = TRUE, line = -1, font = 2)
  }
  grDevices::pdf(paste0(path, ".pdf"), width = 12, height = 5.8)
  draw()
  grDevices::dev.off()
  grDevices::png(paste0(path, ".png"), width = 1800, height = 870, res = 150)
  draw()
  grDevices::dev.off()
}

select_tf_inputs <- function(tf_all, pkn, cutoff = 2, selected_tf_ids = NULL) {
  # An explicit shared set bypasses the native activity cutoff, including real
  # zero scores. Missing/unmeasured scores are never filled with zeros.
  coverage <- data.frame(source = names(tf_all), score = unname(tf_all),
                         present_in_signaling_pkn = names(tf_all) %in% network_nodes(pkn))
  coverage$finite_activity <- is.finite(coverage$score)
  coverage$above_activity_cutoff <- coverage$finite_activity & abs(coverage$score) > cutoff
  coverage$selected_by_shared_set <- if (is.null(selected_tf_ids)) NA else coverage$source %in% selected_tf_ids
  chosen <- if (is.null(selected_tf_ids)) coverage$above_activity_cutoff else coverage$selected_by_shared_set
  coverage$retained <- coverage$present_in_signaling_pkn & coverage$finite_activity & chosen
  if (!is.null(selected_tf_ids) && !all(selected_tf_ids %in% coverage$source[coverage$retained])) {
    stop("Shared TF set contains missing/nonfinite/uncovered activity values.")
  }
  list(input = tf_all[coverage$retained], coverage = coverage)
}

prepare_scoring_network <- function(base_pkn, downstream_input, options) {
  # Network construction depends only on TF identities, never their scores.
  stopifnot(length(downstream_input) >= 3, !anyDuplicated(names(downstream_input)),
            all(is.finite(downstream_input)))
  pkn_uncompressed <- cosmosR:::keep_observable_neighbours(
    base_pkn, n_steps = options$layers, observed_nodes = names(downstream_input))
  pkn_uncompressed <- ordered_pkn(pkn_uncompressed)
  # moon() discards outgoing edges from observed TFs. Save this exact eligible
  # topology, rather than edges that the scoring algorithm will ignore.
  pkn_uncompressed <- pkn_uncompressed[!pkn_uncompressed$source %in% names(downstream_input), ]
  downstream_ids <- names(downstream_input)[names(downstream_input) %in% network_nodes(pkn_uncompressed)]
  if (length(downstream_ids) < 3) stop("Insufficient observable TF coverage after pruning.")
  validate_pkn(pkn_uncompressed)
  compressed <- safe_compression(pkn_uncompressed, downstream_input[downstream_ids],
                                 enabled = options$compress)
  stopifnot(all(downstream_ids %in% network_nodes(compressed$pkn)))
  list(pkn_uncompressed = pkn_uncompressed, compressed = compressed,
       downstream_ids = downstream_ids)
}

score_prepared_network <- function(downstream_input, prepared, n_layers) {
  stopifnot(identical(names(downstream_input), prepared$downstream_ids),
            all(is.finite(downstream_input)))
  # The default cutoff of zero would omit genuine zero-valued observations
  # from level 0. They must survive in the shared-input controlled experiment.
  scores <- cosmosR::moon(
    upstream_input = NULL, downstream_input = downstream_input,
    meta_network = prepared$compressed$pkn, n_layers = n_layers, statistic = "ulm",
    downstream_cutoff = -Inf)
  scores <- as.data.frame(scores[, c("source", "score", "level")])
  stopifnot(!anyDuplicated(scores$source), all(is.finite(scores$score)),
            setequal(scores$source[scores$level == 0], names(downstream_input)))
  scores
}

run_cell_type <- function(cell_type, tf_scores, rna_scaled, base_pkn, options, output_dir,
                          selected_tf_ids = NULL, prepared_network = NULL) {
  started <- proc.time()[[3]]
  slug <- gsub("(^_+|_+$)", "", gsub("[^A-Za-z0-9]+", "_", cell_type))
  out <- file.path(output_dir, slug)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  tf_all <- setNames(tf_scores[, cell_type], rownames(tf_scores))
  rna_input <- setNames(rna_scaled[, cell_type], rownames(rna_scaled))
  rna_input <- rna_input[is.finite(rna_input)]
  selected <- select_tf_inputs(tf_all, base_pkn, options$tf_cutoff, selected_tf_ids)
  coverage <- selected$coverage
  downstream_input <- selected$input
  if (length(downstream_input) < 3) stop("Fewer than 3 covered downstream TFs for ", cell_type)
  # There is no upstream layer, so never call keep_controllable_neighbours().
  prepared <- if (is.null(prepared_network)) {
    prepare_scoring_network(base_pkn, downstream_input, options)
  } else prepared_network
  pkn_uncompressed <- prepared$pkn_uncompressed
  compressed <- prepared$compressed
  pkn_scoring <- compressed$pkn
  stopifnot(all(prepared$downstream_ids %in% names(downstream_input)))
  downstream_input <- downstream_input[prepared$downstream_ids]
  coverage$retained_after_observability <- coverage$source %in% names(downstream_input)
  message("MOON ", options$regulon_label, " / ", cell_type, ": ", length(downstream_input), " TFs; ",
          nrow(pkn_scoring), " edges; ", options$layers, " layers.")
  moon_scores <- score_prepared_network(downstream_input, prepared, options$layers)
  original_scores <- expand_scores(moon_scores, compressed$mapping)
  stopifnot(setequal(original_scores$source[original_scores$level == 0], names(downstream_input)))
  solution <- extract_downstream_anchored(
    original_scores, pkn_uncompressed, cutoff = options$score_cutoff, rna_input = rna_input)
  audit <- audit_solution(solution)
  sensitivity <- do.call(rbind, lapply(sort(unique(c(1, 1.5, 2, options$score_cutoff))), function(cut) {
    s <- extract_downstream_anchored(original_scores, pkn_uncompressed, cut, rna_input)
    audit_solution(s)
    data.frame(score_cutoff = cut, n_nodes = nrow(s$ATT), n_edges = nrow(s$SIF),
               n_observed_TFs = sum(s$ATT$level == 0))
  }))
  top_regulators <- original_scores[original_scores$level > 0, ]
  top_regulators <- top_regulators[order(-abs(top_regulators$score)), ]
  top_regulators$RNA_zscore <- unname(rna_input[top_regulators$source])
  top_regulators$in_display_network <- top_regulators$source %in% solution$ATT$nodes
  outputs <- list(tf_input_coverage = coverage, pkn_uncompressed = pkn_uncompressed,
                  pkn_scoring = pkn_scoring, compression_mapping = compressed$mapping,
                  moon_scores_scoring_ids = moon_scores, moon_scores = original_scores,
                  top_regulators = top_regulators, SIF = solution$SIF, ATT = solution$ATT,
                  network_audit = audit, threshold_sensitivity = sensitivity)
  for (name in names(outputs)) {
    utils::write.csv(outputs[[name]], file.path(out, paste0(name, ".csv")), row.names = FALSE, na = "")
  }
  saveRDS(compressed, file.path(out, "compression.rds"), compress = "gzip")
  plot_diagnostics(original_scores, cell_type, file.path(out, "moon_diagnostics"), options$score_cutoff)
  # Immutable raw PKN provenance is retained once at the parent run directory;
  # this cell's exact final scoring and original-ID PKNs are saved above.
  description <- utils::packageDescription("cosmosR")
  ledger <- list(
    cell_type = cell_type, timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    cosmosR_version = as.character(utils::packageVersion("cosmosR")),
    cosmosR_commit = if (is.null(description$RemoteSha)) NA_character_ else description$RemoteSha,
    decoupleR_version = as.character(utils::packageVersion("decoupleR")),
    input_files = list(directory = normalizePath(options$input_dir),
                       md5 = as.list(tools::md5sum(file.path(options$input_dir,
                                                          c("tf_scores.csv", "rna_scaled.csv"))))),
    mapping = list(upstream_input = NULL, downstream_input = paste(options$regulon_label, "ULM TF activity t scores"),
                   identifiers = "gene symbols and original signaling-complex IDs",
                   RNA = "annotation only; not propagated as protein activity",
                   baseline = "gene-wise Z scores across six pooled cell types; not a treated-control contrast"),
    pkn = list(source = "cosmosR::meta_network", signaling_edges_before_pruning = nrow(base_pkn),
               uncompressed_edges = nrow(pkn_uncompressed), scoring_edges = nrow(pkn_scoring),
               expression_pruning = FALSE,
               expression_pruning_reason = "HVG/TF-selected gene subset cannot establish absence of unmeasured intermediates",
               no_upstream_controllability_prune = TRUE),
    settings = list(statistic = "ulm", n_steps = options$layers, n_layers = options$layers,
                    TF_activity_cutoff = if (is.null(selected_tf_ids)) options$tf_cutoff else NULL,
                    TF_selection = if (is.null(selected_tf_ids)) "native finite, PKN-covered TFs above activity cutoff" else "fixed shared finite, PKN-covered TF identities; all scores including zero",
                    downstream_cutoff_for_level_zero = "-Inf (retain genuine observed zero scores)",
                    display_score_cutoff = options$score_cutoff,
                    thresholds_are_heuristics_not_significance = TRUE),
    compression = list(used = compressed$used, reason = compressed$reason,
                       protected_input_ids_unchanged = TRUE,
                       n_renamed_original_nodes = sum(compressed$mapping$source_original != compressed$mapping$source_scoring)),
    coherence = list(applied = FALSE, iterations = 1L, convergence = "not applicable; fixed signaling PKN",
                     reason = "No separate TF-target transcript layer is added; RNA Z scores are annotations and the TF footprints already use RNA"),
    coverage = list(n_available_TFs = sum(is.finite(tf_all)), n_downstream_TFs = length(downstream_input),
                    input_ids = names(downstream_input),
                    lost_input_ids = coverage$source[(if (is.null(selected_tf_ids)) coverage$above_activity_cutoff else coverage$selected_by_shared_set) &
                                                    !coverage$retained_after_observability]),
    score_levels = as.list(table(original_scores$level)),
    extraction = list(rule = "abs(score)>cutoff; signed edges from level k to k-1; directed route to a level-0 TF",
                      sign_consistency_passed = TRUE, downstream_reachability_passed = TRUE,
                      n_nodes = nrow(solution$ATT), n_edges = nrow(solution$SIF),
                      n_isolated_observed_TFs = sum(!solution$ATT$has_explanatory_edge)),
    elapsed_seconds = unname(proc.time()[[3]] - started),
    outputs = c(paste0(names(outputs), ".csv"), "compression.rds", "moon_diagnostics.pdf", "moon_diagnostics.png", "run_info.json", "sessionInfo.txt"),
    interpretation = "Prior-knowledge-compatible upstream hypotheses; no measured upstream constraints or causal validation")
  jsonlite::write_json(ledger, file.path(out, "run_info.json"), pretty = TRUE, auto_unbox = TRUE, null = "null")
  writeLines(capture.output(utils::sessionInfo()), file.path(out, "sessionInfo.txt"))
  message("Saved ", nrow(solution$SIF), " audited edges and ", nrow(solution$ATT), " nodes to ", out)
  invisible(list(ledger = ledger, scores = original_scores,
                 pkn_scoring = pkn_scoring, pkn_uncompressed = pkn_uncompressed,
                 compression_mapping = compressed$mapping, downstream_input = downstream_input,
                 coverage = coverage))
}

compare_score_tables <- function(collectri, scenicplus, cell_type, design) {
  # Retain union rows and explicit missing values. A node observed as input in
  # one run and inferred upstream in the other is not a valid upstream pair.
  a <- collectri[, c("source", "score", "level")]
  if ("source_scoring" %in% names(collectri)) a$source_scoring_collectri <- collectri$source_scoring
  b <- scenicplus[, c("source", "score", "level")]
  if ("source_scoring" %in% names(scenicplus)) b$source_scoring_scenicplus <- scenicplus$source_scoring
  names(a)[2:3] <- c("score_collectri", "level_collectri")
  names(b)[2:3] <- c("score_scenicplus", "level_scenicplus")
  joined <- merge(a, b, by = "source", all = TRUE, sort = TRUE)
  finite_a <- is.finite(joined$score_collectri) & is.finite(joined$level_collectri)
  finite_b <- is.finite(joined$score_scenicplus) & is.finite(joined$level_scenicplus)
  both <- finite_a & finite_b
  upstream_both <- both & joined$level_collectri > 0 & joined$level_scenicplus > 0
  joined$comparable <- upstream_both & joined$level_collectri == joined$level_scenicplus
  joined$comparison_status <- "comparable"
  joined$comparison_status[!finite_a] <- "not_scored_collectri"
  joined$comparison_status[!finite_b] <- "not_scored_scenicplus"
  joined$comparison_status[both & joined$level_collectri == 0 & joined$level_scenicplus == 0] <- "observed_TF_in_both"
  joined$comparison_status[both & xor(joined$level_collectri == 0, joined$level_scenicplus == 0)] <- "observed_TF_in_one_run"
  joined$comparison_status[upstream_both & joined$level_collectri != joined$level_scenicplus] <- "different_scoring_level"
  for (name in c("delta_score", "delta_abs_score", "rank_abs_collectri", "rank_abs_scenicplus", "rank_gain_scenicplus")) {
    joined[[name]] <- NA_real_
  }
  keep <- which(joined$comparable)
  joined$delta_score[keep] <- joined$score_scenicplus[keep] - joined$score_collectri[keep]
  joined$delta_abs_score[keep] <- abs(joined$score_scenicplus[keep]) - abs(joined$score_collectri[keep])
  # Rank identical eligible nodes at the identical scoring level; avoid mixing
  # ULM tests whose available target/background sets differ between layers.
  for (level in sort(unique(joined$level_collectri[keep]))) {
    idx <- which(joined$comparable & joined$level_collectri == level)
    joined$rank_abs_collectri[idx] <- rank(-abs(joined$score_collectri[idx]), ties.method = "min")
    joined$rank_abs_scenicplus[idx] <- rank(-abs(joined$score_scenicplus[idx]), ties.method = "min")
  }
  joined$rank_gain_scenicplus[keep] <- joined$rank_abs_collectri[keep] - joined$rank_abs_scenicplus[keep]
  joined$cell_type <- cell_type
  joined$design <- design
  joined[, c("cell_type", "design", "source", setdiff(names(joined), c("cell_type", "design", "source")))]
}

same_score_levels <- function(a, b) {
  a <- a[order(a$source), c("source", "level")]
  b <- b[order(b$source), c("source", "level")]
  rownames(a) <- rownames(b) <- NULL
  identical(a, b)
}

parse_options <- function(args) {
  options <- list(cell_type = "Ventricular Cardiomyocyte", all = FALSE, layers = 4L,
                  tf_cutoff = 2, score_cutoff = 1.5, compress = TRUE, mode = "compare",
                  input_dir = NULL, scenicplus_input_dir = NULL, output_dir = NULL)
  mode_explicit <- FALSE
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--help") {
      cat("Usage: Rscript Heart_pseudobulk/03_moon.R [options]\n",
          "  --mode compare|collectri|scenicplus  Default: compare both TF activity sources\n",
          "  --cell-type NAME       Default: Ventricular Cardiomyocyte\n",
          "  --all                  Run all cell types\n",
          "  --layers N             Default: 4 upstream scoring layers\n",
          "  --tf-cutoff NUMBER     Native runs keep abs(TF activity)>2; shared-input control uses all shared scores\n",
          "  --score-cutoff NUMBER  Display abs(MOON score)>1.5 by default\n",
          "  --no-compression       Score the stable uncompressed PKN\n",
          "  --input-dir PATH       Single-mode input, or CollecTRI input in compare mode\n",
          "  --scenicplus-input-dir PATH  Default: results/tf_activity_scenicplus\n",
          "  --output-dir PATH      Single-mode output, or comparison output root\n",
          "Default native outputs: results/moon (CollecTRI), results/moon_scenicplus.\n",
          "Default comparison outputs: results/moon_comparison. No upstream input is supplied.\n", sep = "")
      return(NULL)
    } else if (arg == "--all") {
      options$all <- TRUE
    } else if (arg == "--no-compression") {
      options$compress <- FALSE
    } else {
      known <- c("--cell-type", "--layers", "--tf-cutoff", "--score-cutoff", "--input-dir",
                 "--scenicplus-input-dir", "--output-dir", "--mode")
      if (!arg %in% known || i == length(args)) stop("Unknown/incomplete argument: ", arg, ". Use --help.")
      i <- i + 1L
      name <- gsub("-", "_", substring(arg, 3))
      options[[name]] <- args[[i]]
      if (arg == "--mode") mode_explicit <- TRUE
    }
    i <- i + 1L
  }
  # Backwards compatibility: a custom single input directory still runs once.
  if (!is.null(options$input_dir) && !mode_explicit) options$mode <- "collectri"
  options$layers <- as.numeric(options$layers)
  options$tf_cutoff <- as.numeric(options$tf_cutoff)
  options$score_cutoff <- as.numeric(options$score_cutoff)
  stopifnot(options$mode %in% c("compare", "collectri", "scenicplus"),
            is.finite(options$layers), options$layers >= 1, options$layers <= 10,
            options$layers == as.integer(options$layers),
            is.finite(options$tf_cutoff), options$tf_cutoff >= 0,
            is.finite(options$score_cutoff), options$score_cutoff >= 0)
  options
}

write_pkn_provenance <- function(output_dir, raw_pkn, signaling, base_pkn) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(raw_pkn, file.path(output_dir, "pkn_raw.rds"), compress = "gzip")
  utils::write.csv(data.frame(stage = c("package raw", "symbol/complex signaling", "cleaned signaling"),
                              n_edges = c(nrow(raw_pkn), nrow(signaling), nrow(base_pkn)),
                              n_nodes = c(length(network_nodes(raw_pkn)), length(network_nodes(signaling)),
                                          length(network_nodes(base_pkn)))),
                   file.path(output_dir, "pkn_preparation.csv"), row.names = FALSE)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  options <- parse_options(args)
  if (is.null(options)) return(invisible(NULL))
  for (package in c("cosmosR", "decoupleR", "igraph", "jsonlite")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Missing R package '", package, "'. See Heart_pseudobulk/MOON.md for setup.")
    }
  }
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(file_arg)) stop("Run this script with Rscript (see --help).")
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
  result_dir <- file.path(script_dir, "results")
  collectri_dir <- if (is.null(options$input_dir)) file.path(result_dir, "tf_activity") else options$input_dir
  scenicplus_dir <- if (is.null(options$scenicplus_input_dir)) file.path(result_dir, "tf_activity_scenicplus") else options$scenicplus_input_dir
  if (options$mode == "scenicplus" && !is.null(options$input_dir)) scenicplus_dir <- options$input_dir
  input_dirs <- list(collectri = collectri_dir, scenicplus = scenicplus_dir)
  labels <- c(collectri = "CollecTRI", scenicplus = "SCENIC+")
  methods <- if (options$mode == "compare") names(input_dirs) else options$mode
  inputs <- lapply(input_dirs[methods], function(path) {
    tf_scores <- read_score_matrix(file.path(path, "tf_scores.csv"))
    rna_scaled <- read_score_matrix(file.path(path, "rna_scaled.csv"))
    stopifnot(all(colnames(tf_scores) %in% colnames(rna_scaled)))
    list(tf_scores = tf_scores, rna_scaled = rna_scaled)
  })
  cells_all <- colnames(inputs[[1]]$tf_scores)
  if (length(inputs) == 2L) {
    stopifnot(setequal(cells_all, colnames(inputs[[2]]$tf_scores)))
    stopifnot(isTRUE(all.equal(inputs[[1]]$rna_scaled,
                              inputs[[2]]$rna_scaled[rownames(inputs[[1]]$rna_scaled), cells_all, drop = FALSE],
                              tolerance = 1e-12, check.attributes = TRUE)))
  }
  cell_types <- if (options$all) cells_all else options$cell_type
  if (!all(cell_types %in% cells_all)) stop("Unknown cell type. Available: ", paste(cells_all, collapse = "; "))
  pkn_env <- new.env(parent = emptyenv())
  utils::data("meta_network", package = "cosmosR", envir = pkn_env)
  raw_pkn <- as.data.frame(pkn_env$meta_network[, c("source", "interaction", "target")])
  stopifnot(!anyNA(raw_pkn), all(raw_pkn$interaction %in% c(-1, 1)))
  # RNA is an annotation: the HVG/TF-selected subset cannot establish absence
  # of unmeasured signaling intermediates. Keep the protein/complex layer.
  metabolic_id <- function(x) grepl("^(Metab__|Gene[0-9]+__)", x)
  signaling <- raw_pkn[!metabolic_id(raw_pkn$source) & !metabolic_id(raw_pkn$target), ]
  base_pkn <- ordered_pkn(cosmosR::meta_network_cleanup(signaling))
  validate_pkn(base_pkn)
  native <- list()
  for (method in methods) {
    output_dir <- if (length(methods) == 1L && !is.null(options$output_dir)) options$output_dir else
      file.path(result_dir, if (method == "collectri") "moon" else "moon_scenicplus")
    method_options <- options
    method_options$input_dir <- input_dirs[[method]]
    method_options$regulon_label <- labels[[method]]
    write_pkn_provenance(output_dir, raw_pkn, signaling, base_pkn)
    native[[method]] <- setNames(lapply(cell_types, run_cell_type,
      tf_scores = inputs[[method]]$tf_scores, rna_scaled = inputs[[method]]$rna_scaled,
      base_pkn = base_pkn, options = method_options, output_dir = output_dir), cell_types)
  }
  if (length(methods) == 1L) return(invisible(native))

  comparison_dir <- if (is.null(options$output_dir)) file.path(result_dir, "moon_comparison") else options$output_dir
  dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)
  # Define one global identity set using every supplied profile, so selecting a
  # different cell type does not alter the controlled graph or eligible TFs.
  shared_ids <- intersect(rownames(inputs$collectri$tf_scores), rownames(inputs$scenicplus$tf_scores))
  finite_all <- function(mat) apply(mat[shared_ids, cells_all, drop = FALSE], 1, function(x) all(is.finite(x)))
  shared_ids <- sort(shared_ids[finite_all(inputs$collectri$tf_scores) & finite_all(inputs$scenicplus$tf_scores) &
                                  shared_ids %in% network_nodes(base_pkn)])
  if (length(shared_ids) < 3) stop("Fewer than three globally shared finite, PKN-covered TF activities.")
  # Zero values below are only placeholders for topology construction. They are
  # replaced by each actual method/profile's measured TF activities for scoring.
  prepared <- prepare_scoring_network(base_pkn, setNames(rep(0, length(shared_ids)), shared_ids), options)
  controlled <- list()
  for (method in methods) {
    out <- file.path(comparison_dir, "controlled", method)
    method_options <- options
    method_options$input_dir <- input_dirs[[method]]
    method_options$regulon_label <- labels[[method]]
    write_pkn_provenance(out, raw_pkn, signaling, base_pkn)
    controlled[[method]] <- setNames(lapply(cell_types, run_cell_type,
      tf_scores = inputs[[method]]$tf_scores, rna_scaled = inputs[[method]]$rna_scaled,
      base_pkn = base_pkn, options = method_options, output_dir = out,
      selected_tf_ids = shared_ids, prepared_network = prepared), cell_types)
  }
  controlled_pairs <- native_pairs <- checks <- coverage <- list()
  for (cell in cell_types) {
    a <- controlled$collectri[[cell]]
    b <- controlled$scenicplus[[cell]]
    matched <- list(identical_input_ids = identical(names(a$downstream_input), names(b$downstream_input)),
                    identical_scoring_pkn = identical(a$pkn_scoring, b$pkn_scoring),
                    identical_uncompressed_pkn = identical(a$pkn_uncompressed, b$pkn_uncompressed),
                    identical_compression_mapping = identical(a$compression_mapping, b$compression_mapping),
                    identical_scored_node_levels = same_score_levels(a$scores, b$scores))
    if (!all(unlist(matched))) stop("Controlled comparison topology/level mismatch for ", cell,
                                    "; no score differences will be interpreted.")
    checks[[cell]] <- data.frame(cell_type = cell, as.data.frame(matched),
                                 n_shared_input_TFs = length(a$downstream_input),
                                 n_upstream_nodes = sum(a$scores$level > 0),
                                 n_scoring_edges = nrow(a$pkn_scoring),
                                 scoring_pkn_md5 = unname(tools::md5sum(file.path(comparison_dir, "controlled", "collectri",
                                   gsub("(^_+|_+$)", "", gsub("[^A-Za-z0-9]+", "_", cell)), "pkn_scoring.csv"))))
    controlled_pairs[[cell]] <- compare_score_tables(a$scores, b$scores, cell, "controlled_shared_TFs")
    native_pairs[[cell]] <- compare_score_tables(native$collectri[[cell]]$scores,
                                                  native$scenicplus[[cell]]$scores, cell, "native_full_TFs")
    for (method in methods) {
      for (design in c("native_full_TFs", "controlled_shared_TFs")) {
        run <- if (design == "native_full_TFs") native[[method]][[cell]] else controlled[[method]][[cell]]
        coverage[[paste(cell, method, design)]] <- data.frame(
          cell_type = cell, method = method, design = design,
          n_available_TFs = sum(is.finite(inputs[[method]]$tf_scores[, cell])),
          n_PKN_covered_TFs = sum(run$coverage$finite_activity & run$coverage$present_in_signaling_pkn),
          n_native_above_cutoff_PKN_TFs = sum(run$coverage$finite_activity & run$coverage$present_in_signaling_pkn & run$coverage$above_activity_cutoff),
          n_input_TFs = length(run$downstream_input), n_zero_input_TFs = sum(run$downstream_input == 0),
          n_inferred_upstream_nodes = sum(run$scores$level > 0))
      }
    }
  }
  comparisons <- list(controlled_shared_tf_comparison = do.call(rbind, controlled_pairs),
                       native_full_tf_comparison = do.call(rbind, native_pairs),
                       controlled_topology_checks = do.call(rbind, checks),
                       input_coverage = do.call(rbind, coverage),
                       shared_tf_ids = data.frame(source = shared_ids,
                         retained_after_observability = shared_ids %in% prepared$downstream_ids))
  for (name in names(comparisons)) utils::write.csv(comparisons[[name]],
     file.path(comparison_dir, paste0(name, ".csv")), row.names = FALSE, na = "")
  manifest <- list(
    primary_design = "controlled_shared_TFs",
    question = "Which upstream PKN regulators align more strongly with the same TF identities when their RNA footprints are estimated using SCENIC+ versus CollecTRI?",
    mapping = list(upstream_input = NULL, downstream_input = "signed TF footprint ULM scores", RNA = "annotation only"),
    shared_input_rule = "TFs finite in both methods in all supplied cell types and covered by signaling PKN; identical globally fixed identity set; no activity cutoff; real zero scores retained",
    n_shared_TFs_before_observability = length(shared_ids), n_shared_TFs_scored = length(prepared$downstream_ids),
    shared_input_ids = prepared$downstream_ids,
    controls = list(same_RNA_preprocessing = TRUE, same_TF_input_identities = TRUE,
                    same_final_scoring_graph = TRUE, same_layers = TRUE, same_compression_mapping = TRUE,
                    same_scored_node_levels_verified = TRUE, no_missing_score_imputation = TRUE,
                    n_layers = options$layers, compression_requested = options$compress,
                    same_graph_across_all_requested_cell_types = TRUE),
    signed_delta = "score_scenicplus - score_collectri; positive means more positive, not necessarily stronger alignment",
    magnitude_delta = "abs(score_scenicplus) - abs(score_collectri); positive means stronger alignment to SCENIC+-derived TF footprints",
    rank_scope = "Within cell type and scoring level, absolute-score ranks among common comparable upstream nodes; smaller rank means larger absolute score",
    rank_gain_scenicplus = "rank_abs_collectri - rank_abs_scenicplus; positive means higher ranking with SCENIC+",
    comparison_eligibility = "finite score in both runs, identical level, level>0; level-0 input TFs excluded; incomparable values remain missing",
    native_supplement = "All finite covered TFs with abs(activity)>tf_cutoff for each method. Coverage and final graph may differ; changed levels and observed-vs-inferred nodes are incomparable",
    interpretation = "Descriptive prior-knowledge alignment and sensitivity to TF regulon choice; not probabilities, statistical significance, independent validation, or proof of causal upstream activity. The two ULM footprint methods can differ in target coverage, sign calibration and t-score scale even with fixed TF identities.",
    files = paste0(names(comparisons), ".csv"),
    cell_types = cell_types)
  jsonlite::write_json(manifest, file.path(comparison_dir, "run_info.json"), pretty = TRUE, auto_unbox = TRUE, null = "null")
  message("Saved controlled and native MOON comparisons to ", comparison_dir)
  message("Plot comparisons with: python Heart_pseudobulk/04_compare_moon.py")
  invisible(list(native = native, controlled = controlled, comparisons = comparisons, manifest = manifest))
}

if (sys.nframe() == 0L) main()
