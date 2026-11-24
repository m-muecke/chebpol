test_that("simplex linear interpolation reproduces values at the knots", {
  withr::local_seed(42)
  f <- function(x) exp(-sum(x^2))
  knots <- matrix(runif(200, -1, 1), 2)
  sl <- ipol(f, knots = knots, method = "simplexlinear")
  expect_equal(sl(knots), apply(knots, 2, f))
})

test_that("simplex linear interpolation is exact for a linear function with linear blending", {
  withr::local_seed(42)
  f <- function(x) 1 + 2 * x[1] - 3 * x[2]
  knots <- cbind(matrix(runif(200, -1, 1), 2), c(-1, -1), c(-1, 1), c(1, -1), c(1, 1))
  sl <- ipol(f, knots = knots, method = "simplexlinear")
  s <- matrix(runif(100, -0.9, 0.9), 2)
  expect_equal(sl(s, blend = "linear"), apply(s, 2, f))
})

test_that("simplex linear interpolation approximates a smooth function", {
  withr::local_seed(42)
  f <- function(x) exp(-sum(x^2))
  knots <- cbind(matrix(runif(200, -1, 1), 2), c(-1, -1), c(-1, 1), c(1, -1), c(1, 1))
  sl <- ipol(f, knots = knots, method = "simplexlinear")
  s <- matrix(runif(100, -0.9, 0.9), 2)
  expect_equal(sl(s), apply(s, 2, f), tolerance = 0.1)
})

test_that("simplex linear interpolation only extrapolates when asked", {
  withr::local_seed(42)
  f <- function(x) 1 + 2 * x[1] - 3 * x[2]
  knots <- cbind(matrix(runif(200, -1, 1), 2), c(-1, -1), c(-1, 1), c(1, -1), c(1, 1))
  sl <- ipol(f, knots = knots, method = "simplexlinear")
  outside <- c(2, 2)
  expect_equal(sl(outside), NA_real_)
  expect_equal(sl(outside, epol = TRUE), f(outside))
})

test_that("simplex linear interpolation builds knots from a grid", {
  f <- function(x) exp(-sum(x^2))
  grid <- list(seq(-1, 1, length.out = 9), seq(-1, 1, length.out = 7))
  pts <- grid_points(grid)
  sl <- ipol(f, grid = grid, method = "simplexlinear")
  expect_equal(sl(pts), apply(pts, 2, f))
})

test_that("ipol() errors without knots for simplex linear interpolation", {
  f <- function(x) exp(-sum(x^2))
  expect_snapshot(ipol(f, method = "simplexlinear"), error = TRUE)
})

test_that("simplex linear interpolants error on input of the wrong dimension", {
  withr::local_seed(42)
  f <- function(x) exp(-sum(x^2))
  sl <- ipol(f, knots = matrix(runif(200, -1, 1), 2L), method = "simplexlinear")
  expect_snapshot(error = TRUE, {
    sl(matrix(0, 1, 5))
    sl(matrix(0, 3, 2))
  })
})
