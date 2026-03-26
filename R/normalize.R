#' Run SCTransform normalization via Seurat.
#'
#' @param obj Seurat object (counts in RNA assay)
#' @param variable_features_n Number of variable features to select
#' @return Seurat object with SCTransform normalization applied
#' @export
sctransform_seurat <- function(obj, variable_features_n = 3000L) {
  obj <- Seurat::SCTransform(
    obj,
    vst.flavor    = "v2",
    variable.features.n = as.integer(variable_features_n),
    verbose       = FALSE
  )
  obj
}

#' Run EigenMS normalization via ProteoMM.
#'
#' Requires treatment group labels for ANOVA model.
#' Returns a normalized dense matrix (samples x features).
#'
#' @param mat Dense matrix (features x samples) --- standard R orientation
#' @param group_labels Character vector of treatment group labels per sample
#' @return Dense matrix (features x samples) after EigenMS normalization
#' @export
eigenms_normalize <- function(mat, group_labels) {
  if (!requireNamespace("ProteoMM", quietly = TRUE)) {
    stop("ProteoMM package is required for EigenMS normalization. ",
         "Install via BiocManager::install('ProteoMM')")
  }
  # ProteoMM::eig_norm1 expects log-transformed data in a data frame with
  # a "treatment" column and protein columns. We use eig_norm2 which is
  # the two-step wrapper: eig_norm1 (bias estimation) then eig_norm2 (correction).

  # Input: mat is features x samples. EigenMS wants samples x proteins (features).
  df <- as.data.frame(t(mat))
  df$treatment <- group_labels

  # eig_norm1: estimate bias eigenvectors
  grp <- factor(group_labels)
  # Create the model matrix for ANOVA
  m_ints <- as.matrix(df[, seq_len(ncol(df) - 1)])  # exclude treatment column

  # Use log2 transform for EigenMS (expects log-scale data)
  m_ints_log <- log2(m_ints + 1)

  ep <- ProteoMM::eig_norm1(m = m_ints_log, treatment = grp, prot.info = data.frame(
    prot_name = colnames(m_ints),
    stringsAsFactors = FALSE
  ))

  # eig_norm2: apply bias correction
  en <- ProteoMM::eig_norm2(rv = ep)

  # Extract normalized matrix and convert back from log2
  normalized_log <- en$norm_m
  normalized <- (2^normalized_log) - 1
  normalized[normalized < 0] <- 0

  # Return features x samples matrix
  t(as.matrix(normalized))
}
