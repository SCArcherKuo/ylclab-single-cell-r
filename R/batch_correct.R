#' Compute PCA coordinates and variance ratios from a Seurat reduction.
#'
#' @param obj Seurat object with a PCA reduction.
#' @param reduction Name of the reduction slot (default "pca").
#' @return List with `points` (list of [PC1, PC2] per cell) and `variance_ratio` (numeric vector).
#' @export
compute_pca_coords <- function(obj, reduction = "pca") {
  emb <- Seurat::Embeddings(obj, reduction = reduction)[, 1:2]
  points <- lapply(seq_len(nrow(emb)), function(i) as.numeric(emb[i, ]))
  stdev <- Seurat::Stdev(obj, reduction = reduction)
  variance_ratio <- (stdev^2) / sum(stdev^2)
  list(points = points, variance_ratio = as.numeric(variance_ratio))
}

#' Apply Seurat IntegrateLayers batch correction.
#'
#' @param obj Seurat object with layers already split by batch.
#' @param method One of: "CCA", "RPCA", "JointPCA", "FastMNN".
#' @param batch_column Metadata column identifying batches.
#' @return Integrated Seurat object with a joined RNA assay.
#' @export
batch_correct_seurat <- function(obj, method, batch_column) {
  method_map <- list(
    CCA      = Seurat::CCAIntegration,
    RPCA     = Seurat::RPCAIntegration,
    JointPCA = Seurat::JointPCAIntegration,
    FastMNN  = Seurat::FastMNNIntegration
  )
  integration_fn <- method_map[[method]]
  if (is.null(integration_fn)) {
    stop(paste("Unknown batch correction method:", method,
               "— supported: CCA, RPCA, JointPCA, FastMNN"))
  }
  obj <- Seurat::IntegrateLayers(
    object  = obj,
    method  = integration_fn,
    orig.reduction = "pca",
    new.reduction  = "integrated.dr",
    verbose = FALSE
  )
  obj[["RNA"]] <- Seurat::JoinLayers(obj[["RNA"]])
  obj
}
