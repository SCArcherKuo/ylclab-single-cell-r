#' Run graph-based clustering via Seurat FindClusters.
#'
#' @param obj Seurat object with a neighbor graph already computed (via FindNeighbors).
#' @param method One of: "Louvain-refined", "SLM".
#' @param resolution Clustering resolution (higher = more clusters).
#' @return Seurat object with cluster assignments in \code{seurat_clusters}.
#' @export
cluster_seurat <- function(obj, method, resolution = 1.0) {
  algo <- switch(method,
    "Louvain-refined" = 2L,
    "SLM"             = 3L,
    stop(paste("Unknown clustering method:", method,
               "— supported: Louvain-refined, SLM"))
  )
  obj <- Seurat::FindClusters(obj, algorithm = algo, resolution = resolution, verbose = FALSE)
  obj
}
