test_that("is_number() accepts a single finite number", {
  expect_all_true(vapply(list(1, 2L, -0.5), is_number, logical(1L)))
  expect_all_false(vapply(list(NA_real_, Inf, c(1, 2), "1", NULL, TRUE), is_number, logical(1L)))
})

test_that("is_whole() accepts whole numbers from a lower bound", {
  expect_all_true(vapply(list(0, c(1, 2), 3L, integer(0)), is_whole, logical(1L)))
  expect_all_false(vapply(list(2.5, -1, NA, Inf, "1", NULL), is_whole, logical(1L)))
  expect_all_false(vapply(list(0, -1), is_whole, logical(1L), lower = 1))
})

test_that("ipol() checks the grid of the grid methods", {
  expect_snapshot(error = TRUE, {
    ipol(sum, method = "fh")
    ipol(sum, grid = list(1), method = "multilinear")
    ipol(sum, grid = c(0, 2, 1), method = "general")
  })
})
