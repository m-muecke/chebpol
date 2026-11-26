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

test_that("evalongrid() accepts functions returning integer or logical values", {
  expect_equal(evalongrid(\(x) 1L, 4), array(1, 4))
  expect_equal(evalongrid(\(x) x > 0, 3), array(c(1, 0, 0), 3))
})

test_that("evalongrid() errors when the function value is not usable", {
  calls <- 0
  growing <- function(x) {
    calls <<- calls + 1
    if (calls > 2) c(1, 2) else 1
  }
  expect_snapshot(error = TRUE, {
    evalongrid(\(x) "a", 3)
    evalongrid(growing, 4)
  })
})

test_that("evalongrid() errors on an empty grid dimension", {
  expect_snapshot(error = TRUE, {
    evalongrid(sum, grid = list(numeric(0), 1:3))
    evalongrid(sum, dims = c(2, 0))
    evalongrid(sum, grid = list())
    evalongrid(sum, dims = integer(0))
  })
})

test_that("evalongrid() errors on a grid with too many points", {
  axis <- seq_len(65536)
  expect_snapshot(error = TRUE, {
    evalongrid(sum, grid = list(axis, axis, axis, axis))
    evalongrid(\(x) numeric(8), grid = rep(list(axis[1:1024]), 5))
  })
})

test_that("evalongrid() accepts an integer grid", {
  f <- function(x) sum(x^2)
  expect_equal(evalongrid(f, grid = list(0:4, 0:3)), evalongrid(f, grid = list(seq(0, 4), seq(0, 3))))
})

test_that("evalongrid() calls the function once at each grid point", {
  calls <- 0
  f <- function(x) {
    calls <<- calls + 1
    calls
  }
  expect_equal(evalongrid(f, grid = list(1:3, 1:2)), matrix(1:6, 3L))
  expect_equal(calls, 6)
})

test_that("evalongrid() does not overwrite an argument that the function kept", {
  seen <- list()
  f <- function(x) {
    seen[[length(seen) + 1L]] <<- x
    sum(x)
  }
  evalongrid(f, grid = list(1:3, 1:2))
  expect_equal(seen, list(c(1, 1), c(2, 1), c(3, 1), c(1, 2), c(2, 2), c(3, 2)))
})
