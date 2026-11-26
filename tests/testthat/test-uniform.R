test_that("uniform grid interpolation reproduces values at the grid points", {
  f <- function(x) 1 / (1 + 25 * x^2)
  grid <- seq(-1, 1, length.out = 15)
  uc <- ipol(f(grid), method = "uniform")
  expect_equal(uc(grid), f(grid))
})

test_that("the uniform grid mapping is exact at integer and half-integer phases", {
  expect_identical(uniform_to_chebyshev(c(-4, -2, 0, 2, 4), 2), c(0, 1, 0, -1, 0))
})

test_that("uniform grid interpolation is the same from values and from a function", {
  withr::local_seed(42)
  f <- function(x) 1 / (1 + 25 * x^2)
  uc_val <- ipol(f(seq(-1, 1, length.out = 15)), method = "uniform")
  uc_fun <- ipol(f, dims = 15, intervals = c(-1, 1), method = "uniform")
  a <- runif(100, -1, 1)
  expect_equal(uc_fun(a), uc_val(a))
  expect_equal(uc_fun(a), f(a), tolerance = 0.05)
})

test_that("uniform grid interpolation maps each dimension of a matrix of points", {
  withr::local_seed(42)
  f <- function(x) exp(x[1] - x[2] / 4)
  intervals <- list(c(0, 1), c(3, -2))
  uc <- ipol(f, dims = c(12, 9), intervals = intervals, method = "uniform")
  m <- rbind(runif(30), runif(30, -2, 3))
  expect_equal(uc(m), apply(m, 2L, uc))
  expect_equal(uc(m), apply(m, 2L, f), tolerance = 0.05)
  expect_equal(uc(m[, 1L, drop = FALSE]), uc(m[, 1L]))
  ud <- ipol(array(runif(6 * 4), c(6, 4)), method = "uniform")
  m <- matrix(runif(20, -1, 1), 2)
  expect_equal(ud(m), apply(m, 2L, ud))
})

test_that("uniform grid interpolation from values accepts a single interval as a vector", {
  f <- function(x) exp(-x^2)
  values <- f(seq(0, 2, length.out = 15))
  uc <- ipol(values, intervals = c(0, 2), method = "uniform")
  expect_equal(uc(0.7), f(0.7), tolerance = 1e-4)
  expect_equal(uc(0.7), ipol(values, intervals = list(c(0, 2)), method = "uniform")(0.7))
  expect_equal(attr(uc, "domain"), attr(ipol(values, intervals = list(c(0, 2)), method = "uniform"), "domain"))
})

test_that("uniform grid interpolation works on a reversed interval", {
  uc <- ipol(exp, dims = 8, intervals = c(1, 0), method = "uniform")
  expect_equal(uc(0.2), ipol(exp, dims = 8, intervals = c(0, 1), method = "uniform")(0.2))
})

test_that("ipol() errors on malformed uniform intervals", {
  f <- function(x) sum(x)
  expect_snapshot(error = TRUE, {
    ipol(f, dims = c(5, 5), intervals = list(c(0, 1)), method = "uniform")
    ipol(f, dims = c(5, 5), intervals = list(c(0, 1), c(0, 2), c(0, 3)), method = "uniform")
    ipol(array(0, c(3, 4)), intervals = list(c(0, 1)), method = "uniform")
    ipol(array(0, c(3, 4)), intervals = "a", method = "uniform")
  })
})
