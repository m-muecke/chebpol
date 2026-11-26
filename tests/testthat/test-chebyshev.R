test_that("chebknots() returns Chebyshev knots", {
  expect_equal(chebknots(3), list(c(sqrt(3) / 2, 0, -sqrt(3) / 2)))
  expect_identical(chebknots(3)[[1L]][2L], 0)
  expect_equal(lengths(chebknots(c(x = 3, y = 4, z = 5))), c(x = 3L, y = 4L, z = 5L))
  expect_equal(range(chebknots(4, c(2, 3))[[1]]), 2.5 + c(-0.5, 0.5) * cospi(1 / 8))
})

test_that("chebknots() errors on malformed dims and intervals", {
  expect_snapshot(error = TRUE, {
    chebknots(2.5)
    chebknots(c(2, 2), list(c(0, 1)))
    chebknots(2, list(c(0, 1), c(5, 6)))
    chebknots(3, c(0, 1, 2))
  })
})

test_that("chebcoef() keeps dimensions and scales the constant term", {
  withr::local_seed(42)
  a <- array(rnorm(24), c(2, 3, 4))
  expect_equal(dim(chebcoef(a)), dim(a))
  expect_equal(chebcoef(a)[1], mean(a))
  expect_equal(chebcoef(a, dct = TRUE)[1], 8 * sum(a))
})

test_that("chebcoef() errors when dct is not TRUE or FALSE", {
  expect_snapshot(error = TRUE, {
    chebcoef(1:3, dct = NA)
    chebcoef(1:3, dct = c(TRUE, FALSE))
  })
})

test_that("Chebyshev interpolation reproduces values at the knots", {
  f <- function(x) exp(-sum(x^2))
  dims <- c(8, 7, 6)
  ch <- ipol(f, dims = dims, method = "chebyshev")
  pts <- grid_points(chebknots(dims))
  expect_equal(ch(pts), apply(pts, 2, f))
})

test_that("Chebyshev interpolation approximates a smooth function", {
  withr::local_seed(42)
  f <- function(x) exp(-sum(x^2))
  ch <- ipol(f, dims = c(8, 7, 6), method = "chebyshev")
  s <- matrix(runif(300, -1, 1), 3)
  expect_equal(ch(s), apply(s, 2, f), tolerance = 0.01)
})

test_that("Chebyshev interpolation works on custom intervals", {
  f <- function(x) exp(-sum(x^2))
  ch <- ipol(f, dims = c(8, 7, 6), intervals = list(c(1, 2), c(1, 4), c(1, 3)))
  s <- cbind(c(1.4, 2.3, 1.9), c(1.01, 1.01, 1.01))
  expect_equal(ch(s), apply(s, 2, f), tolerance = 0.01)
})

test_that("Chebyshev interpolation handles long one-dimensional grids", {
  skip_if_not(has_fftw())
  f <- function(x) ifelse(x == 0, 0, sin(1 / x))
  ch <- ipol(f, dims = 50000, method = "chebyshev")
  expect_equal(ch(c(0.03, 0.031)), f(c(0.03, 0.031)), tolerance = 0.01)
})

test_that("general grid interpolation reproduces values at the grid points", {
  f <- function(x) exp(-sum(x^2))
  su <- seq(0, 1, length.out = 11)
  grid <- list(su, su^2, su^3)
  fv <- evalongrid(f, grid = grid)
  ch <- ipol(fv, grid = grid, method = "general")
  expect_equal(ch(grid_points(grid)), c(fv))
})

test_that("general grid interpolation approximates a smooth function", {
  withr::local_seed(42)
  s <- seq(0, 1, length.out = 10)^3
  ch <- ipol(exp(s), grid = list(s), method = "general")
  r <- runif(50)
  expect_equal(ch(r), exp(r), tolerance = 0.05)
})

test_that("general grid interpolation maps each dimension of a matrix of points", {
  withr::local_seed(42)
  f <- function(x) exp(-sum(x^2)) + sin(3 * x[1]) * x[2]
  grid <- list(seq(0, 1, length.out = 12)^2, seq(-2, 3, length.out = 9)^3 / 9)
  ch <- ipol(f, grid = grid, method = "general")
  m <- rbind(runif(30), runif(30, -1, 3))
  expect_equal(ch(m), apply(m, 2L, ch))
  expect_equal(ch(m[, 1L, drop = FALSE]), ch(m[, 1L]))
})

test_that("general grid interpolation is constant along a dimension with a single grid point", {
  f <- function(x) x[1]^2 + x[2]
  grid <- list(c(0, 0.3, 1), 5)
  ch <- ipol(f, grid = grid, method = "general")
  expect_equal(ch(c(0.3, 5)), f(c(0.3, 5)))
  expect_equal(ch(cbind(c(0.3, 5), c(1, 7))), c(f(c(0.3, 5)), f(c(1, 5))))
})

test_that("chebcoef() and chebeval() accept integer input", {
  values <- array(seq_len(20), c(5L, 4L))
  coef <- chebcoef(values + 0)
  expect_equal(chebcoef(values), coef)
  expect_equal(chebeval(c(0L, 1L), coef), chebeval(c(0, 1), coef))
  expect_equal(chebeval(c(0.1, 0.2), values), chebeval(c(0.1, 0.2), values + 0))
})

test_that("chebeval() maps points into the given intervals", {
  f <- function(x) exp(-sum(x^2))
  intervals <- list(c(0, 2), c(1, 3))
  coef <- chebcoef(evalongrid(f, c(8, 8), intervals))
  ch <- ipol(f, dims = c(8, 8), intervals = intervals, method = "chebyshev")
  x <- cbind(c(0.5, 1.5), c(1.2, 2.2), c(1.9, 1.1))
  expect_equal(chebeval(x, coef, intervals), ch(x))
  expect_equal(chebeval(x[, 1L], coef, intervals), ch(x[, 1L]))
})

test_that("ipol() errors on malformed Chebyshev intervals", {
  val <- array(0, c(3, 4))
  expect_snapshot(error = TRUE, {
    ipol(val, intervals = "a", method = "chebyshev")
    ipol(val, intervals = list(c(0, 1)), method = "chebyshev")
    ipol(val, intervals = list(c(0, 1), c(0, 1, 2)), method = "chebyshev")
    ipol(sum, dims = 3, intervals = list(c(0, 1), c(0, 2)), method = "chebyshev")
    ipol(sin, dims = 2.5, method = "chebyshev")
  })
})

test_that("chebeval() errors on malformed intervals", {
  coef <- chebcoef(matrix(1:9, 3))
  expect_snapshot(error = TRUE, {
    chebeval(c(0.1, 0.2), coef, intervals = list(c(0, 1)))
    chebeval(c(0.1, 0.2), coef, intervals = list(c(0, 1), c(0, 1, 2)))
  })
})

test_that("chebeval() accepts a single interval as a vector", {
  f <- function(x) exp(-x^2)
  coef <- chebcoef(evalongrid(f, 12, c(0, 2)))
  expect_equal(chebeval(0.7, coef, c(0, 2)), f(0.7), tolerance = 1e-6)
  expect_equal(chebeval(0.7, coef, list(c(0, 2))), chebeval(0.7, coef, c(0, 2)))
})

test_that("chebeval() treats coefficients without a dim attribute as one-dimensional", {
  coef <- chebcoef(sin(chebknots(10)[[1L]]))
  expect_equal(chebeval(0.3, as.numeric(coef)), chebeval(0.3, coef))
  expect_equal(chebeval(0.3, c(1, 0.5, 0.25)), 1 + 0.5 * 0.3 + 0.25 * (2 * 0.3^2 - 1))
})

test_that("chebeval() errors on empty coefficients", {
  expect_snapshot(error = TRUE, {
    chebeval(0, numeric())
    chebeval(c(0, 0), matrix(numeric(), nrow = 0, ncol = 2))
  })
})
