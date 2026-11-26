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

test_that("oldstalker interpolation with a non-linear blend is linear outside the grid", {
  grid <- seq(0, 1, length.out = 6)
  x <- c(-0.5, -0.01, 1.01, 1.5)
  expected <- c(sin(0.2) * x[1:2] / 0.2, sin(1) + (sin(1) - sin(0.8)) * (x[3:4] - 1) / 0.2)
  for (blend in c("cubic", "sigmoid", "parodic", "square")) {
    st <- ipol(sin, grid = grid, method = "oldstalker", blend = blend)
    expect_equal(st(x), expected)
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

test_that("ipol() errors when the values do not match the stalker grid", {
  grid <- list(seq(-1, 1, length.out = 9))
  expect_snapshot(error = TRUE, {
    ipol(1:5, grid = grid, method = "stalker")
    ipol(1:5, grid = grid, method = "hstalker")
    ipol(1:5, grid = grid, method = "oldstalker")
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

test_that("the hyperbolic stalker does not depend on the scale of the values or the grid", {
  f <- function(x) x^3 - x
  grid <- seq(0, 1, length.out = 6)
  x <- seq(0, 1, length.out = 101)
  expected <- ipol(f, grid = list(grid), method = "hstalker")(x)
  for (scale in c(1e-150, 1e-7, 1e5, 1e150)) {
    values <- ipol(\(x) scale * f(x), grid = list(grid), method = "hstalker")
    expect_equal(values(x) / scale, expected)
    stretched <- ipol(\(x) f(x / scale), grid = list(grid * scale), method = "hstalker")
    expect_equal(stretched(x * scale), expected)
  }
})

test_that("the stalker spline does not depend on the scale of the values or the grid", {
  f <- function(x) exp(2 * x) * sin(5 * x)
  x <- seq(0, 1, length.out = 101)
  for (grid in list(seq(0, 1, length.out = 6), c(0, 0.1, 0.35, 0.5, 0.8, 1))) {
    expected <- ipol(f, grid = list(grid), method = "stalker")(x)
    for (scale in c(1e-150, 1e-7, 1e5, 1e150)) {
      values <- ipol(\(x) scale * f(x), grid = list(grid), method = "stalker")
      expect_equal(values(x) / scale, expected)
    }
    for (scale in c(1e-7, 1e5)) {
      stretched <- ipol(\(x) f(x / scale), grid = list(grid * scale), method = "stalker")
      expect_equal(stretched(x * scale), expected)
    }
  }
})

test_that("the old stalker spline with a solved degree does not depend on the scale of the grid", {
  skip_if_not(havegsl())
  f <- function(x) exp(2 * x) * sin(5 * x)
  grid <- c(0, 0.1, 0.35, 0.5, 0.8, 1)
  x <- seq(0, 1, length.out = 101)
  expected <- ipol(f, grid = list(grid), method = "oldstalker", k = NA)(x)
  for (scale in c(1e-9, 1e9)) {
    stretched <- ipol(\(x) f(x / scale), grid = list(grid * scale), method = "oldstalker", k = NA)
    expect_equal(stretched(x * scale), expected)
  }
})

test_that("the old hyperbolic stalker spline does not depend on the scale of the values or the grid", {
  f <- function(x) exp(2 * x) * sin(5 * x)
  grid <- c(0, 0.1, 0.35, 0.5, 0.8, 1)
  x <- seq(0, 1, length.out = 101)
  expected <- ipol(f, grid = list(grid), method = "oldstalker", k = 0)(x)
  for (scale in c(1e-12, 1e-6, 1e6, 1e12)) {
    values <- ipol(\(x) scale * f(x), grid = list(grid), method = "oldstalker", k = 0)
    expect_equal(values(x) / scale, expected)
    stretched <- ipol(\(x) f(x / scale), grid = list(grid * scale), method = "oldstalker", k = 0)
    expect_equal(stretched(x * scale), expected)
  }
})

test_that("the hyperbolic stalker has no pole next to a nearly flat neighbor", {
  hst <- ipol(c(1e-6, 0, 1, 2), grid = list(c(-1, 0, 1, 2)), method = "hstalker")
  x <- c(0.999997, 0.999998, 0.999999)
  expect_equal(hst(x, blend = "linear"), x, tolerance = 1e-5)
})

test_that("the hyperbolic stalker is finite at the grid points when neighboring values tie", {
  values <- cbind(c(1, 0, 0, 1), c(2, 2, 0, 3), c(0, 0, 0.1 + 0.2, 0.3))
  grid <- list(1:4, 1:3)
  hst <- ipol(values, grid = grid, method = "hstalker")
  x <- grid_points(grid)
  for (blend in c("cubic", "linear", "sigmoid", "parodic", "square")) {
    expect_equal(hst(x, blend = blend), c(values))
  }
  expect_false(anyNA(hst(x, blend = "mean")))
  hst1 <- ipol(c(1, 0, 0, 1), grid = list(1:4), method = "hstalker")
  expect_equal(hst1(1:4, blend = "mean"), c(0.5, 0, 0, 0.5))
})
