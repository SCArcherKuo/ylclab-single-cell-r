library(testthat)
library(ylclabscm)

test_that("write_sentinel writes valid JSON to the specified path", {
  tmp <- tempfile(fileext = ".json")
  on.exit(unlink(tmp))

  data <- list(status = "completed", serial_number = "SN001", dataset_uid = 42L)
  write_sentinel(data, tmp)

  parsed <- jsonlite::fromJSON(tmp)
  expect_equal(parsed$status, "completed")
  expect_equal(parsed$serial_number, "SN001")
  expect_equal(parsed$dataset_uid, 42L)
})

test_that("write_sentinel creates parent directories if missing", {
  tmp <- file.path(tempdir(), "nested", "dir", "out.json")
  on.exit(unlink(dirname(dirname(tmp)), recursive = TRUE))

  write_sentinel(list(status = "completed"), tmp)
  expect_true(file.exists(tmp))
})
