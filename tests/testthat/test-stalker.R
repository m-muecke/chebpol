test_that("oldstalker interpolant errors on input of the wrong dimension", {
  f <- function(x) exp(-sum(x^2))
  grid <- list(seq(-1, 1, length.out = 9), seq(-1, 1, length.out = 9))
  st <- ipol(f, grid = grid, method = "oldstalker")
  expect_snapshot(st(matrix(0, 3, 4)), error = TRUE)
})

test_that("stalker splines reproduce values at the grid points", {
  f <- function(x) exp(-sum(x^2))
  grid <- list(seq(-1, 1, length.out = 9), seq(-1, 1, length.out = 7))
  pts <- grid_points(grid)
  for (method in c("stalker", "hstalker", "oldstalker")) {
    st <- ipol(f, grid = grid, method = method)
    expect_equal(st(pts), apply(pts, 2, f))
  }
})

test_that("stalker splines reproduce values on a non-uniform grid", {
  skip_if_not(havegsl())
  f <- function(x) exp(-sum(x^2))
  grid <- list(seq(-1, 1, length.out = 9), seq(-1, 1, length.out = 7)^3)
  pts <- grid_points(grid)
  for (method in c("stalker", "hstalker", "oldstalker")) {
    st <- ipol(f, grid = grid, method = method)
    expect_equal(st(pts), apply(pts, 2, f))
  }
})

test_that("stalker splines approximate a smooth function", {
  withr::local_seed(42)
  f <- function(x) exp(-sum(x^2))
  grid <- list(seq(-1, 1, length.out = 9), seq(-1, 1, length.out = 7))
  s <- matrix(runif(200, -1, 1), 2)
  for (method in c("stalker", "hstalker", "oldstalker")) {
    st <- ipol(f, grid = grid, method = method)
    expect_equal(st(s), apply(s, 2, f), tolerance = 0.05)
  }
})

test_that("stalker splines preserve monotonicity", {
  x <- seq(0, 1, length.out = 8)
  y <- c(0, 0.01, 0.02, 0.5, 0.9, 0.95, 0.96, 1)
  s <- seq(0, 1, length.out = 500)
  for (method in c("stalker", "hstalker")) {
    st <- ipol(y, grid = list(x), method = method)
    expect_gte(min(diff(st(s))), -1e-12)
  }
})

test_that("stalker splines give the same result for a vector, a matrix and several threads", {
  f <- function(x) exp(-sum(x^2))
  grid <- list(seq(-1, 1, length.out = 9), seq(-1, 1, length.out = 7))
  st <- ipol(f, grid = grid, method = "stalker")
  p <- c(0.3, 0.2)
  expect_equal(st(matrix(p, 2)), st(p))
  expect_equal(st(p, threads = 2), st(p))
})

test_that("blend in ipol() sets the default blending of the interpolant", {
  withr::local_seed(42)
  f <- function(x) exp(-sum(x^2))
  grid <- list(seq(-1, 1, length.out = 9), seq(-1, 1, length.out = 7))
  s <- matrix(runif(200, -1, 1), 2)
  st <- ipol(f, grid = grid, method = "stalker")
  st_linear <- ipol(f, grid = grid, method = "stalker", blend = "linear")
  expect_equal(eval(formals(st_linear)$blend)[1], "linear")
  expect_equal(st_linear(s), st(s, blend = "linear"))
})

test_that("interpolating blends reproduce values at the grid points", {
  f <- function(x) exp(-sum(x^2))
  grid <- list(seq(-1, 1, length.out = 9), seq(-1, 1, length.out = 7))
  pts <- grid_points(grid)
  st <- ipol(f, grid = grid, method = "stalker")
  for (blend in c("cubic", "linear", "sigmoid", "parodic", "square")) {
    expect_equal(st(pts, blend = blend), apply(pts, 2, f))
  }
})

test_that("ipol() errors on invalid stalker input", {
  f <- function(x) exp(-sum(x^2))
  grid <- list(seq(-1, 1, length.out = 9))
  expect_snapshot(error = TRUE, {
    ipol(f, method = "stalker")
    ipol(f, grid = list(c(1, 0, 2)), method = "stalker")
    ipol(f, grid = grid, method = "stalker", blend = "nope")
  })
})

test_that("ipol() errors on a stalker grid dimension with a single point", {
  grid <- list(c(0, 1, 2), 5)
  expect_snapshot(error = TRUE, {
    ipol(sum, grid = grid, method = "stalker")
    ipol(sum, grid = grid, method = "hstalker")
    ipol(sum, grid = grid, method = "oldstalker")
  })
})

test_that("stalker interpolants error on input of the wrong dimension", {
  f <- function(x) exp(-sum(x^2))
  grid <- list(seq(-1, 1, length.out = 9), seq(-1, 1, length.out = 7))
  st <- ipol(f, grid = grid, method = "stalker")
  hst <- ipol(f, grid = grid, method = "hstalker")
  expect_snapshot(error = TRUE, {
    st(matrix(0, 1, 5))
    st(matrix(0, 3, 2))
    hst(matrix(0, 1, 5))
    hst(matrix(0, 3, 2))
  })
})
