#' Load a feature chunk from COO parquet via DuckDB.
#'
#' @param parquet_path Path to COO sparse parquet (row, col, value).
#' @param feature_indices Integer vector of 0-based feature indices to load.
#' @param n_samples Number of samples (columns in the dense matrix).
#' @param sample_index_map Named integer vector mapping col values to 1-based positions.
#' @return Dense matrix (n_samples x chunk_size).
#' @export
load_feature_chunk <- function(parquet_path, feature_indices, n_samples, sample_index_map) {
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  idx_str <- paste(feature_indices, collapse = ",")
  query <- sprintf(
    "SELECT col, row, COALESCE(value, 0.0) AS value FROM read_parquet('%s') WHERE row IN (%s)",
    parquet_path, idx_str
  )
  df <- DBI::dbGetQuery(con, query)

  chunk_size <- length(feature_indices)
  feat_local <- setNames(seq_along(feature_indices), as.character(feature_indices))

  mat <- matrix(0.0, nrow = n_samples, ncol = chunk_size)
  for (k in seq_len(nrow(df))) {
    si <- sample_index_map[as.character(df$col[k])]
    fi <- feat_local[as.character(df$row[k])]
    if (!is.na(si) && !is.na(fi)) {
      mat[si, fi] <- df$value[k]
    }
  }

  mat
}


#' Run chunked differential abundance using Seurat per chunk.
#'
#' @param parquet_path COO parquet path.
#' @param row_names Feature names.
#' @param column_names Sample names.
#' @param group_labels Named character vector: sample -> group.
#' @param method Seurat test method.
#' @param reference Reference group.
#' @param chunk_size Features per chunk.
#' @return List with groups, all_results (data.frame with group, name, score, logfoldchange, pval).
#' @export
diff_abundance_chunked <- function(parquet_path, row_names, column_names,
                                    group_labels, method, reference, chunk_size) {
  valid_methods <- c("bimod", "roc", "negbinom", "poisson", "LR", "MAST", "DESeq2")
  test_use <- if (tolower(method) == "roc") "roc" else method
  if (!test_use %in% valid_methods) {
    stop(paste("Unknown method:", method, "— supported:", paste(valid_methods, collapse = ", ")))
  }

  if (test_use == "MAST" && !requireNamespace("MAST", quietly = TRUE))
    stop("MAST method requires the MAST Bioconductor package")
  if (test_use == "DESeq2" && !requireNamespace("DESeq2", quietly = TRUE))
    stop("DESeq2 method requires the DESeq2 Bioconductor package")

  n_features <- length(row_names)
  n_samples <- length(column_names)
  groups <- sort(unique(group_labels))

  # Build sample index map
  sample_index_map <- setNames(seq_along(column_names), as.character(column_names))
  # Also map integer indices
  for (i in seq_along(column_names)) {
    sample_index_map[as.character(i - 1L)] <- i
  }

  n_chunks <- ceiling(n_features / chunk_size)
  all_results <- list()

  for (ci in seq_len(n_chunks)) {
    start <- (ci - 1L) * chunk_size
    end <- min(start + chunk_size, n_features)
    feat_idx <- seq(start, end - 1L)
    feat_names <- row_names[feat_idx + 1L]

    message(sprintf("Chunk %d/%d: features [%d:%d]", ci, n_chunks, start, end - 1L))

    chunk_mat <- load_feature_chunk(parquet_path, feat_idx, n_samples, sample_index_map)
    colnames(chunk_mat) <- feat_names
    rownames(chunk_mat) <- column_names

    # Build Seurat object for this chunk
    t_mat <- t(chunk_mat)  # Seurat wants features x cells
    obj <- Seurat::CreateSeuratObject(counts = t_mat)
    obj <- Seurat::SetAssayData(obj, layer = "data", new.data = Seurat::GetAssayData(obj, layer = "counts"))
    obj$group <- group_labels[colnames(obj)]
    valid_cells <- !is.na(obj$group)
    if (sum(valid_cells) < ncol(obj)) {
      obj <- subset(obj, cells = colnames(obj)[valid_cells])
    }
    Seurat::Idents(obj) <- "group"

    # Run FindAllMarkers on this chunk
    markers <- tryCatch(
      Seurat::FindAllMarkers(obj, test.use = test_use, only.pos = FALSE, verbose = FALSE),
      error = function(e) {
        message(sprintf("FindAllMarkers failed on chunk %d: %s", ci, conditionMessage(e)))
        NULL
      }
    )

    if (!is.null(markers) && nrow(markers) > 0) {
      # Compute scores
      if (tolower(method) == "roc" && "myAUC" %in% colnames(markers)) {
        scores <- markers$myAUC
      } else {
        pvals_safe <- pmax(markers$p_val, .Machine$double.xmin)
        scores <- -log10(pvals_safe) * sign(markers$avg_log2FC)
      }

      chunk_df <- data.frame(
        group = as.character(markers$cluster),
        name = as.character(markers$gene),
        score = as.numeric(scores),
        logfoldchange = as.numeric(markers$avg_log2FC),
        pval = as.numeric(markers$p_val),
        stringsAsFactors = FALSE
      )
      all_results[[ci]] <- chunk_df
    }

    rm(chunk_mat, t_mat, obj)
    gc()
  }

  combined <- do.call(rbind, all_results)
  list(groups = groups, all_results = combined)
}


#' Apply BH FDR and format into results_by_group.
#'
#' @param groups Sorted group names.
#' @param all_results data.frame with group, name, score, logfoldchange, pval.
#' @param method Test method name.
#' @return List with groups and results_by_group.
#' @export
adjust_and_format_results <- function(groups, all_results, method) {
  all_results$pval_adj <- p.adjust(all_results$pval, method = "BH")

  results_by_group <- list()
  for (g in groups) {
    gdf <- all_results[all_results$group == g, , drop = FALSE]
    ord <- order(abs(gdf$score), decreasing = TRUE)
    gdf <- gdf[ord, , drop = FALSE]

    results_by_group[[g]] <- list(
      names          = as.character(gdf$name),
      scores         = as.numeric(gdf$score),
      logfoldchanges = as.numeric(gdf$logfoldchange),
      pvals          = as.numeric(gdf$pval),
      pvals_adj      = as.numeric(gdf$pval_adj)
    )
  }

  list(groups = groups, results_by_group = results_by_group)
}
