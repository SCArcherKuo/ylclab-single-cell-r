library(testthat)
library(ylclabscm)

# Helper: build a tiny in-memory Seurat object (10 features × 20 cells, 2 batches)
make_test_seurat <- function() {
  set.seed(42)
  mat <- matrix(rpois(200, lambda = 5), nrow = 10, ncol = 20)
  rownames(mat) <- paste0("feat", seq_len(nrow(mat)))
  colnames(mat) <- paste0("cell", seq_len(ncol(mat)))
  meta <- data.frame(
    batch = rep(c("B1", "B2"), each = 10),
    bio   = rep(c("T1", "T2"), times = 10),
    row.names = colnames(mat)
  )
  obj <- Seurat::CreateSeuratObject(counts = mat, meta.data = meta)
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, verbose = FALSE)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, npcs = 5, verbose = FALSE)
  obj
}

test_that("compute_pca_coords returns a list with points and variance_ratio", {
  skip_if_not_installed("Seurat")
  obj <- make_test_seurat()
  result <- compute_pca_coords(obj, reduction = "pca")
  expect_type(result, "list")
  expect_named(result, c("points", "variance_ratio"))
  expect_equal(length(result$points), ncol(obj))
  expect_true(length(result$variance_ratio) > 0)
  expect_equal(sum(result$variance_ratio), 1.0, tolerance = 1e-6)
})

test_that("batch_correct_seurat returns a Seurat object for CCA method", {
  skip_if_not_installed("Seurat")
  obj <- make_test_seurat()
  # Split layers by batch to satisfy IntegrateLayers requirement
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj$batch)
  result <- batch_correct_seurat(obj, method = "CCA", batch_column = "batch")
  expect_s4_class(result, "Seurat")
})

test_that("batch_correct_seurat errors on unknown method", {
  skip_if_not_installed("Seurat")
  obj <- make_test_seurat()
  expect_error(
    batch_correct_seurat(obj, method = "unknown", batch_column = "batch"),
    "Unknown batch correction method"
  )
})
