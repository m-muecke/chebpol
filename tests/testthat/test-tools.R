test_that("interpolants accept a point, a matrix, or a vector of several points", {
  grid <- list(c(0, 1, 2), c(0, 1))
  ml <- ipol(evalongrid(sum, grid = grid), grid = grid, method = "multilinear")
  m <- rbind(c(0.5, 1.5, 2), c(0.2, 0.4, 1))
  expect_equal(ml(c(0.5, 0.2)), 0.7)
  expect_equal(ml(m), colSums(m))
  expect_equal(ml(c(1L, 1L)), 2)
  expect_snapshot(out <- ml(c(m)))
  expect_equal(out, colSums(m))
  ml1 <- ipol(c(0, 2, 4), grid = list(c(0, 1, 2)), method = "multilinear")
  expect_equal(ml1(c(0.5, 1.5)), c(1, 3))
})

test_that("interpolants error when the input length does not fit the number of variables", {
  grid <- list(c(0, 1), c(0, 1), c(0, 1))
  ml <- ipol(evalongrid(sum, grid = grid), grid = grid, method = "multilinear")
  expect_snapshot(ml(c(0.1, 0.2)), error = TRUE)
})

test_that("ipol() sets the default blend of the interpolant", {
  grid <- list(c(0, 1, 2))
  ml <- ipol(c(0, 1, 0), grid = grid, method = "multilinear", blend = "cubic")
  expect_equal(eval(formals(ml)$blend), c("cubic", "linear", "sigmoid", "parodic", "square"))
  expect_equal(ml(0.25), ml(0.25, blend = "cubic"))
  expect_equal(attr(ml, "arity"), 1)
})

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
  })
})

test_that("evalongrid() accepts an integer grid", {
  f <- function(x) sum(x^2)
  expect_equal(evalongrid(f, grid = list(0:4, 0:3)), evalongrid(f, grid = list(seq(0, 4), seq(0, 3))))
})

test_that("evalongrid() does not overwrite an argument that the function kept", {
  seen <- list()
  f <- function(x) {
    seen[[length(seen) + 1L]] <<- x
    sum(x)
  }
  evalongrid(f, grid = list(1:3, 1:2))
  expect_equal(unique(seen), list(c(1, 1), c(2, 1), c(3, 1), c(1, 2), c(2, 2), c(3, 2)))
})
