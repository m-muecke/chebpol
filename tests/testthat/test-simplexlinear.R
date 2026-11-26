test_that("simplex linear interpolation reproduces values at the knots", {
  skip_if_not_installed("geometry")
  withr::local_seed(42)
  f <- function(x) exp(-sum(x^2))
  knots <- matrix(runif(200, -1, 1), 2)
  sl <- ipol(f, knots = knots, method = "simplexlinear")
  expect_equal(sl(knots), apply(knots, 2, f))
})

test_that("simplex linear interpolation is exact for a linear function with linear blending", {
  skip_if_not_installed("geometry")
  withr::local_seed(42)
  f <- function(x) 1 + 2 * x[1] - 3 * x[2]
  knots <- cbind(matrix(runif(200, -1, 1), 2), c(-1, -1), c(-1, 1), c(1, -1), c(1, 1))
  sl <- ipol(f, knots = knots, method = "simplexlinear")
  s <- matrix(runif(100, -0.9, 0.9), 2)
  expect_equal(sl(s, blend = "linear"), apply(s, 2, f))
})

test_that("simplex linear interpolation approximates a smooth function", {
  skip_if_not_installed("geometry")
  withr::local_seed(42)
  f <- function(x) exp(-sum(x^2))
  knots <- cbind(matrix(runif(200, -1, 1), 2), c(-1, -1), c(-1, 1), c(1, -1), c(1, 1))
  sl <- ipol(f, knots = knots, method = "simplexlinear")
  s <- matrix(runif(100, -0.9, 0.9), 2)
  expect_equal(sl(s), apply(s, 2, f), tolerance = 0.1)
})

test_that("simplex linear interpolation only extrapolates when asked", {
  skip_if_not_installed("geometry")
  withr::local_seed(42)
  f <- function(x) 1 + 2 * x[1] - 3 * x[2]
  knots <- cbind(matrix(runif(200, -1, 1), 2), c(-1, -1), c(-1, 1), c(1, -1), c(1, 1))
  sl <- ipol(f, knots = knots, method = "simplexlinear")
  outside <- c(2, 2)
  expect_equal(sl(outside), NA_real_)
  expect_equal(sl(outside, epol = TRUE), f(outside))
})

test_that("simplex linear interpolants error when epol is not TRUE or FALSE", {
  skip_if_not_installed("geometry")
  sl <- ipol(sum, knots = cbind(c(0, 0), c(1, 0), c(0, 1)), method = "simplexlinear")
  expect_snapshot(error = TRUE, {
    sl(c(2, 2), epol = NA)
    sl(c(2, 2), epol = "yes")
  })
})

test_that("simplex linear interpolation finds the same simplices with and without the grid index", {
  skip_if_not_installed("geometry")
  withr::local_seed(42)
  for (d in 1:4) {
    knots <- matrix(runif(d * 60), d)
    sl <- ipol(\(x) sum(sin(3 * x)), knots = knots, method = "simplexlinear")
    # without the index, as in interpolants saved by earlier versions, all simplices are scanned
    scan <- sl
    environment(scan) <- list2env(as.list(environment(sl)), parent = parent.env(environment(sl)))
    environment(scan)$adata <- environment(scan)$adata[1:3]
    x <- matrix(runif(d * 500, -0.1, 1.1), d)
    expect_identical(sl(x), scan(x))
    expect_identical(sl(x, epol = TRUE), scan(x, epol = TRUE))
  }
})

test_that("simplex linear interpolation builds knots from a grid", {
  skip_if_not_installed("geometry")
  f <- function(x) exp(-sum(x^2))
  grid <- list(seq(-1, 1, length.out = 9), seq(-1, 1, length.out = 7))
  pts <- grid_points(grid)
  sl <- ipol(f, grid = grid, method = "simplexlinear")
  expect_equal(sl(pts), apply(pts, 2, f))
})

test_that("simplex linear interpolation accepts one-dimensional knots as a vector", {
  skip_if_not_installed("geometry")
  sl <- ipol(sqrt, knots = c(0, 1, 4, 9), method = "simplexlinear")
  expect_equal(sl(c(0.5, 2.5)), c(0.5, 1.5))
})

test_that("ipol() errors without knots for simplex linear interpolation", {
  f <- function(x) exp(-sum(x^2))
  expect_snapshot(ipol(f, method = "simplexlinear"), error = TRUE)
})

test_that("simplex linear interpolants error on input of the wrong dimension", {
  skip_if_not_installed("geometry")
  withr::local_seed(42)
  f <- function(x) exp(-sum(x^2))
  sl <- ipol(f, knots = matrix(runif(200, -1, 1), 2L), method = "simplexlinear")
  expect_snapshot(error = TRUE, {
    sl(matrix(0, 1, 5))
    sl(matrix(0, 3, 2))
  })
})

test_that("simplex linear interpolation forwards extra arguments to the function", {
  skip_if_not_installed("geometry")
  withr::local_seed(42)
  f <- function(x, a) a * sum(x^2)
  knots <- matrix(runif(200, -1, 1), 2L)
  sl <- ipol(f, knots = knots, method = "simplexlinear", a = 3)
  expect_equal(sl(knots), apply(knots, 2L, f, a = 3))
})

test_that("simplex linear blends reproduce values at the knots", {
  skip_if_not_installed("geometry")
  withr::local_seed(42)
  f <- function(x) exp(-sum(x^2))
  knots <- cbind(matrix(runif(200, -1, 1), 2L), c(-1, -1), c(-1, 1), c(1, -1), c(1, 1))
  sl <- ipol(f, knots = knots, method = "simplexlinear")
  for (blend in c("cubic", "linear", "sigmoid", "parodic", "square")) {
    expect_equal(sl(knots, blend = blend), apply(knots, 2L, f))
  }
})

test_that("square blending returns the value at the vertex with the largest barycentric coordinate", {
  skip_if_not_installed("geometry")
  knots <- cbind(c(0, 0), c(1, 0), c(0, 1))
  sl <- ipol(c(10, 20, 40), knots = knots, method = "simplexlinear")
  x <- cbind(c(0.1, 0.1), c(0.3, 0.25), c(0.25, 0.45), c(0.7, 0.2))
  expect_equal(sl(x, blend = "square"), c(10, 10, 40, 20))
})

test_that("simplex linear blends are finite on the edges of a grid triangulation", {
  skip_if_not_installed("geometry")
  f <- function(x) exp(-sum(x^2))
  grid <- list(seq(-1, 1, length.out = 5), seq(-1, 1, length.out = 5))
  sl <- ipol(f, grid = grid, method = "simplexlinear")
  x <- grid_points(list(seq(-1, 1, length.out = 9), seq(-1, 1, length.out = 9)))
  for (blend in c("cubic", "linear", "sigmoid", "parodic", "square", "mean")) {
    expect_false(anyNA(sl(x, blend = blend)))
  }
})

test_that("simplex linear blends work in three dimensions", {
  skip_if_not_installed("geometry")
  withr::local_seed(42)
  f <- function(x) sum(x^2)
  corners <- grid_points(list(c(-1, 1), c(-1, 1), c(-1, 1)))
  knots <- cbind(matrix(runif(300, -1, 1), 3L), corners)
  sl <- ipol(f, knots = knots, method = "simplexlinear")
  x <- matrix(runif(900, -0.9, 0.9), 3L)
  for (blend in c("cubic", "linear", "sigmoid", "parodic", "square")) {
    expect_equal(sl(knots, blend = blend), apply(knots, 2L, f))
  }
  for (blend in c("cubic", "linear", "sigmoid", "parodic", "square", "mean")) {
    expect_false(anyNA(sl(x, blend = blend)))
  }
})

test_that("simplex linear interpolation returns NA for missing input", {
  skip_if_not_installed("geometry")
  withr::local_seed(42)
  f <- function(x) exp(-sum(x^2))
  knots <- cbind(matrix(runif(200, -1, 1), 2L), c(-1, -1), c(-1, 1), c(1, -1), c(1, 1))
  sl <- ipol(f, knots = knots, method = "simplexlinear")
  x <- cbind(c(NaN, 0.1), c(0.1, NA), c(0.1, 0.2))
  for (blend in c("cubic", "linear", "sigmoid", "parodic", "square", "mean")) {
    expect_equal(is.na(sl(x, blend = blend)), c(TRUE, TRUE, FALSE))
  }
})

test_that("ipol() errors when values and knots differ in number for simplex linear interpolation", {
  skip_if_not_installed("geometry")
  knots <- matrix(seq(0, 1, length.out = 20), 2L)
  expect_snapshot(error = TRUE, {
    ipol(1:3, knots = knots, method = "simplexlinear")
  })
})

test_that("ipol() errors for simplex linear interpolation without geometry", {
  skip_if(requireNamespace("geometry", quietly = TRUE))
  expect_snapshot(error = TRUE, {
    ipol(sum, knots = matrix(runif(20), 2L), method = "simplexlinear")
  })
})
