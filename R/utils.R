#' Write metadata sentinel JSON to a file.
#'
#' @param data Named list of metadata fields.
#' @param path Absolute path for the output JSON file.
#' @export
write_sentinel <- function(data, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  jsonlite::write_json(data, path, auto_unbox = TRUE, pretty = FALSE)
}


#' Load COO sparse parquet with name TSVs that may be multi-column or headerless.
#'
#' Handles two TSV formats:
#' - Multi-column with 'index': e.g., 'sample\tindex' or 'mode\tmass\tindex'
#' - Single-column headerless: just names, one per line (row order = index)
#'
#' @param parquet_path Path to COO parquet (row, col, value)
#' @param row_path Path to row (feature) names TSV
#' @param column_path Path to column (sample) names TSV
#' @return List with sparse Matrix 'mat', character vectors 'row_names', 'column_names'
#' @export
load_coo_parquet <- function(parquet_path, row_path, column_path) {
  sparse_df <- arrow::read_parquet(parquet_path)

  read_name_tsv <- function(path, name_col) {
    # Peek at first line to detect format
    first_line <- readLines(path, n = 1)
    has_tab <- grepl("\t", first_line)

    if (has_tab) {
      # Multi-column format with header
      df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      if ("index" %in% colnames(df)) {
        df <- df[order(df$index), ]
        name_cols <- setdiff(colnames(df), "index")
        if (length(name_cols) == 1) {
          return(df[[name_cols]])
        } else {
          return(apply(df[, name_cols, drop = FALSE], 1, paste, collapse = "_"))
        }
      }
    }

    # Single-column headerless format (or fallback)
    lines <- readLines(path)
    return(lines)
  }

  row_names    <- read_name_tsv(row_path, "feature")
  column_names <- read_name_tsv(column_path, "sample")

  mat <- Matrix::sparseMatrix(
    i    = sparse_df$row + 1L,
    j    = sparse_df$col + 1L,
    x    = sparse_df$value,
    dims = c(length(row_names), length(column_names))
  )
  rownames(mat) <- row_names
  colnames(mat) <- column_names

  list(mat = mat, row_names = row_names, column_names = column_names)
}
