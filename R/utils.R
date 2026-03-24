#' Write metadata sentinel JSON to a file.
#'
#' @param data Named list of metadata fields.
#' @param path Absolute path for the output JSON file.
#' @export
write_sentinel <- function(data, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  jsonlite::write_json(data, path, auto_unbox = TRUE, pretty = FALSE)
}
