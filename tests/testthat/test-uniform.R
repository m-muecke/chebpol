test_that("uniform grid interpolation reproduces values at the grid points", {
  f <- function(x) 1 / (1 + 25 * x^2)
  grid <- seq(-1, 1, length.out = 15)
  uc <- ipol(f(grid), method = "uniform")
  expect_equal(uc(grid), f(grid))
})

test_that("the uniform grid mapping is exact at integer and half-integer phases", {
  expect_identical(ugm(c(-4, -2, 0, 2, 4), 2), c(0, 1, 0, -1, 0))
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
