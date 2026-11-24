test_that("Floater-Hormann interpolation approximates a smooth function", {
  withr::local_seed(42)
  f <- function(x) sum(log(1 + x^2))
  grid <- list(seq(-1, 1, length.out = 50), seq(-1, 1, length.out = 30)^3L, chebknots(20)[[1]])
  fh <- ipol(f, grid = grid, method = "fh", k = c(3, 4, 4))
  x <- matrix(runif(300), 3)
  expect_equal(fh(x), apply(x, 2, f), tolerance = 1e-3)
  expect_equal(fh(x, threads = 2), fh(x))
})

test_that("Floater-Hormann interpolation is continuous at the grid points", {
  f <- function(x) sum(log(1 + x^2))
  grid <- list(seq(-1, 1, length.out = 50), seq(-1, 1, length.out = 30)^3L, chebknots(20)[[1]])
  fh <- ipol(f, grid = grid, method = "fh", k = c(3, 4, 4))
  p <- vapply(grid, \(x) x[4], numeric(1L))
  expect_equal(fh(p), f(p))
  expect_equal(fh(p + 1e-9), f(p), tolerance = 1e-7)
})

test_that("Floater-Hormann interpolation forwards extra arguments to the function", {
  f <- function(x, a) a * sum(x^2)
  grid <- list(seq(0, 1, length.out = 6), seq(0, 1, length.out = 6))
  fh <- ipol(f, grid = grid, k = 2, method = "fh", a = 3)
  expect_equal(fh(c(0.2, 0.4)), f(c(0.2, 0.4), a = 3))
})

test_that("Floater-Hormann interpolation works at the smallest and largest degree", {
  f <- function(x) exp(sum(x))
  grid <- list(seq(0, 1, length.out = 6), seq(0, 1, length.out = 4))
  pts <- grid_points(grid)
  for (k in list(0, c(5, 3))) {
    fh <- ipol(f, grid = grid, method = "fh", k = k)
    expect_equal(fh(pts), apply(pts, 2L, f))
    expect_equal(anyNA(fh(c(0.3, 0.6))), FALSE)
  }
})

test_that("ipol() rejects a Floater-Hormann degree outside the valid range", {
  f <- function(x) exp(sum(x))
  grid <- list(seq(0, 1, length.out = 6), seq(0, 1, length.out = 4))
  expect_snapshot(error = TRUE, {
    ipol(f, grid = grid, method = "fh", k = -1)
    ipol(f, grid = grid, method = "fh", k = 4)
    ipol(f, grid = grid, method = "fh", k = c(6, 3))
  })
})
