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

run_cell_type <- function(cell_type, tf_scores, rna_scaled, base_pkn, options, output_dir) {
  started <- proc.time()[[3]]
  slug <- gsub("(^_+|_+$)", "", gsub("[^A-Za-z0-9]+", "_", cell_type))
  out <- file.path(output_dir, slug)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  tf_all <- setNames(tf_scores[, cell_type], rownames(tf_scores))
  rna_input <- setNames(rna_scaled[, cell_type], rownames(rna_scaled))
  rna_input <- rna_input[is.finite(rna_input)]
  coverage <- data.frame(source = names(tf_all), score = unname(tf_all),
                         present_in_signaling_pkn = names(tf_all) %in% network_nodes(base_pkn))
  coverage$above_activity_cutoff <- is.finite(coverage$score) & abs(coverage$score) > options$tf_cutoff
  coverage$retained <- coverage$present_in_signaling_pkn & coverage$above_activity_cutoff
  downstream_input <- tf_all[coverage$retained]
  if (length(downstream_input) < 3) stop("Fewer than 3 covered downstream TFs for ", cell_type)
  stopifnot(!anyDuplicated(names(downstream_input)), all(is.finite(downstream_input)))

  # There is no upstream seed layer, so do not call keep_controllable_neighbours.
  pkn_uncompressed <- cosmosR:::keep_observable_neighbours(
    base_pkn, n_steps = options$layers, observed_nodes = names(downstream_input))
  pkn_uncompressed <- ordered_pkn(pkn_uncompressed)
  downstream_input <- downstream_input[names(downstream_input) %in% network_nodes(pkn_uncompressed)]
  # Edges outgoing from level-0 TFs are ignored inside moon(). Remove them here
  # explicitly so the saved scoring PKN is exactly the topology eligible to score.
  pkn_uncompressed <- pkn_uncompressed[!pkn_uncompressed$source %in% names(downstream_input), ]
  downstream_input <- downstream_input[names(downstream_input) %in% network_nodes(pkn_uncompressed)]
  coverage$retained_after_observability <- coverage$source %in% names(downstream_input)
  if (length(downstream_input) < 3) stop("Insufficient observable TF coverage after pruning.")
  validate_pkn(pkn_uncompressed)
  compressed <- safe_compression(pkn_uncompressed, downstream_input, enabled = options$compress)
  pkn_scoring <- compressed$pkn
  stopifnot(all(names(downstream_input) %in% network_nodes(pkn_scoring)))
  message("MOON ", cell_type, ": ", length(downstream_input), " TFs; ",
          nrow(pkn_scoring), " edges; ", options$layers, " layers.")
  moon_scores <- cosmosR::moon(
    upstream_input = NULL, downstream_input = downstream_input,
    meta_network = pkn_scoring, n_layers = options$layers, statistic = "ulm")
  moon_scores <- as.data.frame(moon_scores[, c("source", "score", "level")])
  stopifnot(!anyDuplicated(moon_scores$source), all(is.finite(moon_scores$score)),
            setequal(moon_scores$source[moon_scores$level == 0], names(downstream_input)))
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
    mapping = list(upstream_input = NULL, downstream_input = "CollecTRI ULM TF activity t scores",
                   identifiers = "gene symbols and original signaling-complex IDs",
                   RNA = "annotation only; not propagated as protein activity",
                   baseline = "gene-wise Z scores across six pooled cell types; not a treated-control contrast"),
    pkn = list(source = "cosmosR::meta_network", signaling_edges_before_pruning = nrow(base_pkn),
               uncompressed_edges = nrow(pkn_uncompressed), scoring_edges = nrow(pkn_scoring),
               expression_pruning = FALSE,
               expression_pruning_reason = "HVG/TF-selected gene subset cannot establish absence of unmeasured intermediates",
               no_upstream_controllability_prune = TRUE),
    settings = list(statistic = "ulm", n_steps = options$layers, n_layers = options$layers,
                    TF_activity_cutoff = options$tf_cutoff, display_score_cutoff = options$score_cutoff,
                    thresholds_are_heuristics_not_significance = TRUE),
    compression = list(used = compressed$used, reason = compressed$reason,
                       protected_input_ids_unchanged = TRUE,
                       n_renamed_original_nodes = sum(compressed$mapping$source_original != compressed$mapping$source_scoring)),
    coherence = list(applied = FALSE, iterations = 1L, convergence = "not applicable; fixed signaling PKN",
                     reason = "No separate TF-target transcript layer is added; RNA Z scores are annotations and the TF footprints already use RNA"),
    coverage = list(n_available_TFs = sum(is.finite(tf_all)), n_downstream_TFs = length(downstream_input),
                    input_ids = names(downstream_input),
                    lost_input_ids = coverage$source[coverage$above_activity_cutoff & !coverage$retained_after_observability]),
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
  invisible(ledger)
}

parse_options <- function(args) {
  options <- list(cell_type = "Ventricular Cardiomyocyte", all = FALSE, layers = 4L,
                  tf_cutoff = 2, score_cutoff = 1.5, compress = TRUE,
                  input_dir = NULL, output_dir = NULL)
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--help") {
      cat("Usage: Rscript Heart_pseudobulk/03_moon.R [options]\n",
          "  --cell-type NAME       Default: Ventricular Cardiomyocyte\n",
          "  --all                  Run all cell types\n",
          "  --layers N             Default: 4 upstream scoring layers\n",
          "  --tf-cutoff NUMBER     Keep downstream abs(TF activity)>2 by default\n",
          "  --score-cutoff NUMBER  Display abs(MOON score)>1.5 by default\n",
          "  --no-compression       Score the stable uncompressed PKN\n",
          "  --input-dir PATH       Default: Heart_pseudobulk/results/tf_activity\n",
          "  --output-dir PATH      Default: Heart_pseudobulk/results/moon\n", sep = "")
      return(NULL)
    } else if (arg == "--all") {
      options$all <- TRUE
    } else if (arg == "--no-compression") {
      options$compress <- FALSE
    } else {
      known <- c("--cell-type", "--layers", "--tf-cutoff", "--score-cutoff", "--input-dir", "--output-dir")
      if (!arg %in% known || i == length(args)) stop("Unknown/incomplete argument: ", arg, ". Use --help.")
      i <- i + 1L
      name <- gsub("-", "_", substring(arg, 3))
      options[[name]] <- args[[i]]
    }
    i <- i + 1L
  }
  options$layers <- as.numeric(options$layers)
  options$tf_cutoff <- as.numeric(options$tf_cutoff)
  options$score_cutoff <- as.numeric(options$score_cutoff)
  stopifnot(is.finite(options$layers), options$layers >= 1, options$layers <= 10,
            options$layers == as.integer(options$layers),
            is.finite(options$tf_cutoff), options$tf_cutoff >= 0,
            is.finite(options$score_cutoff), options$score_cutoff >= 0)
  options
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
  input_dir <- if (is.null(options$input_dir)) file.path(script_dir, "results", "tf_activity") else options$input_dir
  output_dir <- if (is.null(options$output_dir)) file.path(script_dir, "results", "moon") else options$output_dir
  options$input_dir <- input_dir
  tf_scores <- read_score_matrix(file.path(input_dir, "tf_scores.csv"))
  rna_scaled <- read_score_matrix(file.path(input_dir, "rna_scaled.csv"))
  stopifnot(all(colnames(tf_scores) %in% colnames(rna_scaled)))
  cell_types <- if (options$all) colnames(tf_scores) else options$cell_type
  if (!all(cell_types %in% colnames(tf_scores))) {
    stop("Unknown cell type. Available: ", paste(colnames(tf_scores), collapse = "; "))
  }
  pkn_env <- new.env(parent = emptyenv())
  utils::data("meta_network", package = "cosmosR", envir = pkn_env)
  raw_pkn <- as.data.frame(pkn_env$meta_network[, c("source", "interaction", "target")])
  stopifnot(!anyNA(raw_pkn), all(raw_pkn$interaction %in% c(-1, 1)))
  # Retain the protein/complex signaling layer. Metabolites and reaction-specific
  # enzyme copies are outside this RNA-only teaching question.
  metabolic_id <- function(x) grepl("^(Metab__|Gene[0-9]+__)", x)
  signaling <- raw_pkn[!metabolic_id(raw_pkn$source) & !metabolic_id(raw_pkn$target), ]
  base_pkn <- ordered_pkn(cosmosR::meta_network_cleanup(signaling))
  validate_pkn(base_pkn)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(raw_pkn, file.path(output_dir, "pkn_raw.rds"), compress = "gzip")
  utils::write.csv(data.frame(stage = c("package raw", "symbol/complex signaling", "cleaned signaling"),
                              n_edges = c(nrow(raw_pkn), nrow(signaling), nrow(base_pkn)),
                              n_nodes = c(length(network_nodes(raw_pkn)), length(network_nodes(signaling)),
                                          length(network_nodes(base_pkn)))),
                   file.path(output_dir, "pkn_preparation.csv"), row.names = FALSE)
  ledgers <- lapply(cell_types, run_cell_type, tf_scores = tf_scores, rna_scaled = rna_scaled,
                   base_pkn = base_pkn, options = options, output_dir = output_dir)
  invisible(ledgers)
}

if (sys.nframe() == 0L) main()
