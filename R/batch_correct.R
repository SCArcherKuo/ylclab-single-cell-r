#' Compute PCA coordinates and variance ratios from a Seurat reduction.
#'
#' @param obj Seurat object with a PCA reduction.
#' @param reduction Name of the reduction slot (default "pca").
#' @return List with `points` (list of [PC1, PC2] per cell) and `variance_ratio`
#'   (numeric vector; fraction of variance explained among the PCs computed, not
#'   total data variance — values sum to 1.0 across all returned PCs).
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
#' @param batch_column Metadata column identifying batches. The object's RNA
#'   layers must already be split by this column before calling this function
#'   (Seurat v5 infers batch structure from split layer names, not this argument).
#' @return Integrated Seurat object with a joined RNA assay.
#' @export
batch_correct_seurat <- function(obj, method, batch_column) {
  # Resolve function name lazily so unavailable methods (e.g. FastMNNIntegration
  # in Seurat core) don't crash jobs that use a different method.
  fn_name <- switch(method,
    CCA      = "CCAIntegration",
    RPCA     = "RPCAIntegration",
    JointPCA = "JointPCAIntegration",
    FastMNN  = "FastMNNIntegration",
    stop(paste("Unknown batch correction method:", method,
               "— supported: CCA, RPCA, JointPCA, FastMNN"))
  )
  integration_fn <- utils::getFromNamespace(fn_name, "Seurat")
  npcs_used <- ncol(Seurat::Embeddings(obj, "pca"))
  obj <- Seurat::IntegrateLayers(
    object  = obj,
    method  = integration_fn,
    orig.reduction = "pca",
    new.reduction  = "integrated.dr",
    dims    = seq_len(npcs_used),
    verbose = FALSE
  )
  obj[["RNA"]] <- Seurat::JoinLayers(obj[["RNA"]])
  obj
}
