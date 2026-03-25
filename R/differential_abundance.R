#' Run differential abundance analysis via Seurat FindAllMarkers.
#'
#' @param obj Seurat object with identities set to the grouping variable.
#' @param method One of: "bimod", "ROC", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2".
#' @return Data frame from FindAllMarkers with columns: p_val, avg_log2FC, pct.1, pct.2,
#'   p_val_adj, cluster, gene (and myAUC for ROC).
#' @export
diff_abundance_seurat <- function(obj, method) {
  valid_methods <- c("bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2")
  # FindAllMarkers uses lowercase "roc" for test.use but we accept "ROC" from the user
  test_use <- if (tolower(method) == "roc") "roc" else method

  if (!test_use %in% valid_methods) {
    stop(paste("Unknown differential abundance method:", method,
               "— supported:", paste(valid_methods, collapse = ", ")))
  }

  # Check for optional Bioconductor packages
  if (test_use == "MAST" && !requireNamespace("MAST", quietly = TRUE)) {
    stop("MAST method requires the MAST Bioconductor package")
  }
  if (test_use == "DESeq2" && !requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 method requires the DESeq2 Bioconductor package")
  }

  markers <- Seurat::FindAllMarkers(obj, test.use = test_use, only.pos = FALSE, verbose = FALSE)
  markers
}

#' Convert FindAllMarkers output to the platform's results_by_group format.
#'
#' @param markers Data frame from FindAllMarkers.
#' @param method The method name (used to determine score computation).
#' @param n_top Number of top features per group.
#' @return List with `groups` (character vector) and `results_by_group` (named list).
#' @export
format_markers_to_results <- function(markers, method, n_top = 50L) {
  groups <- sort(unique(as.character(markers$cluster)))
  results_by_group <- list()

  for (g in groups) {
    gdf <- markers[as.character(markers$cluster) == g, , drop = FALSE]

    # Compute scores
    if (tolower(method) == "roc" && "myAUC" %in% colnames(gdf)) {
      scores <- gdf$myAUC
    } else {
      # score = -log10(p_val) * sign(avg_log2FC), handle p_val == 0
      pvals <- pmax(gdf$p_val, .Machine$double.xmin)
      scores <- -log10(pvals) * sign(gdf$avg_log2FC)
    }

    # Sort by absolute score descending
    ord <- order(abs(scores), decreasing = TRUE)
    gdf <- gdf[ord, , drop = FALSE]
    scores <- scores[ord]

    # Take top N
    n <- min(n_top, nrow(gdf))
    gdf <- gdf[seq_len(n), , drop = FALSE]
    scores <- scores[seq_len(n)]

    results_by_group[[g]] <- list(
      names          = as.character(gdf$gene),
      scores         = as.numeric(scores),
      logfoldchanges = as.numeric(gdf$avg_log2FC),
      pvals          = as.numeric(gdf$p_val),
      pvals_adj      = as.numeric(gdf$p_val_adj)
    )
  }

  list(groups = groups, results_by_group = results_by_group)
}
