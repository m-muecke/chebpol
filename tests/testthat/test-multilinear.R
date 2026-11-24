test_that("multilinear interpolation reproduces values at the grid points", {
  f <- function(x) exp(-sum(x^2))
  su <- seq(0, 1, length.out = 11)
  grid <- list(su, su^2, su^3)
  fv <- evalongrid(f, grid = grid)
  lip <- ipol(fv, grid = grid, method = "multilinear")
  expect_equal(lip(grid_points(grid)), c(fv))
})

test_that("multilinear interpolation approximates a smooth function", {
  withr::local_seed(42)
  f <- function(x) exp(-sum(x^2))
  su <- seq(0, 1, length.out = 11)
  grid <- list(su, su^2, su^3)
  lip <- ipol(f, grid = grid, method = "multilinear")
  s <- matrix(runif(300), 3)
  expect_equal(lip(s), apply(s, 2, f), tolerance = 0.05)
})

test_that("multilinear blends reproduce values at the grid points", {
  f <- function(x) exp(-sum(x^2))
  grid <- list(seq(-1, 1, length.out = 9), seq(-1, 1, length.out = 7))
  pts <- grid_points(grid)
  lip <- ipol(f, grid = grid, method = "multilinear")
  for (blend in c("linear", "cubic", "sigmoid", "parodic", "square")) {
    expect_equal(lip(pts, blend = blend), apply(pts, 2, f))
  }
})

test_that("ipol() rejects a blend that multilinear interpolation does not support", {
  f <- function(x) exp(-sum(x^2))
  grid <- list(seq(-1, 1, length.out = 9), seq(-1, 1, length.out = 7))
  expect_snapshot(ipol(f, grid = grid, method = "multilinear", blend = "mean"), error = TRUE)
})
