#' Run ICA dimension reduction via Seurat.
#'
#' @param obj Seurat object with variable features and scaled data
#' @param nics Number of independent components to compute
#' @return Seurat object with ICA reduction
#' @export
reduce_ica_seurat <- function(obj, nics = 50L) {
  obj <- Seurat::RunICA(obj, nics = as.integer(nics), verbose = FALSE)
  obj
}

#' Run SPCA (Supervised PCA) dimension reduction via Seurat.
#'
#' Requires SNN graph: runs FindNeighbors first, then RunSPCA.
#'
#' @param obj Seurat object with variable features and scaled data, PCA already computed
#' @param npcs Number of supervised PCs to compute
#' @param n_neighbors Number of neighbors for SNN graph
#' @return Seurat object with SPCA reduction
#' @export
reduce_spca_seurat <- function(obj, npcs = 50L, n_neighbors = 20L) {
  npcs <- as.integer(npcs)
  n_neighbors <- as.integer(n_neighbors)
  # Build SNN graph (required for SPCA supervision)
  obj <- Seurat::FindNeighbors(obj, dims = 1:npcs, k.param = n_neighbors, verbose = FALSE)
  obj <- Seurat::RunSPCA(obj, npcs = npcs, graph = "RNA_snn", verbose = FALSE)
  obj
}
