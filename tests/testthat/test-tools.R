test_that("evalongrid() evaluates scalar and vector valued functions", {
  dims <- c(8, 7, 6)
  f <- function(x) exp(-sum(x^2))
  g <- function(x) c(sum(x), prod(x), exp(-sum(x^2)))
  expect_equal(dim(evalongrid(f, dims)), dims)
  expect_equal(dim(evalongrid(g, dims)), c(3, dims))
})

test_that("evalongrid() evaluates on a supplied grid", {
  grid <- list(c(0, 1), c(0, 2, 4))
  expect_equal(evalongrid(sum, grid = grid), outer(grid[[1]], grid[[2]], `+`), ignore_attr = TRUE)
})
