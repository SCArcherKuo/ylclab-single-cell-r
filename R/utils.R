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

  # Read a name TSV and return list(names, index_map)
  # index_map: original_index -> 1-based consecutive position
  read_name_tsv <- function(path) {
    first_line <- readLines(path, n = 1)
    has_tab <- grepl("\t", first_line)

    if (has_tab) {
      # Multi-column format with header (from preparation)
      df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      if ("index" %in% colnames(df)) {
        df <- df[order(df$index), ]
        name_cols <- setdiff(colnames(df), "index")
        if (length(name_cols) == 1) {
          names <- df[[name_cols]]
        } else {
          names <- apply(df[, name_cols, drop = FALSE], 1, paste, collapse = "_")
        }
        # Build map: original index -> 1-based position
        index_map <- setNames(seq_along(df$index), as.character(df$index))
        return(list(names = names, index_map = index_map))
      }
    }

    # Single-column headerless format (from analysis step outputs)
    # Row order = 0-based index
    lines <- readLines(path)
    index_map <- setNames(seq_along(lines), as.character(seq(0L, length(lines) - 1L)))
    return(list(names = lines, index_map = index_map))
  }

  row_info <- read_name_tsv(row_path)
  col_info <- read_name_tsv(column_path)

  row_names    <- row_info$names
  column_names <- col_info$names

  # Remap sparse indices through the name TSV mapping
  # Original indices may be non-consecutive (e.g., after filtration)
  mapped_i <- row_info$index_map[as.character(sparse_df$row)]
  mapped_j <- col_info$index_map[as.character(sparse_df$col)]

  # Drop entries that don't map (shouldn't happen, but defensive)
  valid <- !is.na(mapped_i) & !is.na(mapped_j)
  if (sum(!valid) > 0) {
    warning(sprintf("load_coo_parquet: dropped %d entries with unmapped indices", sum(!valid)))
  }

  mat <- Matrix::sparseMatrix(
    i    = as.integer(mapped_i[valid]),
    j    = as.integer(mapped_j[valid]),
    x    = sparse_df$value[valid],
    dims = c(length(row_names), length(column_names))
  )
  rownames(mat) <- row_names
  colnames(mat) <- column_names

  list(mat = mat, row_names = row_names, column_names = column_names)
}
