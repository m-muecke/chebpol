test_that("chebknots() returns Chebyshev knots", {
  expect_equal(chebknots(3), list(c(sqrt(3) / 2, 0, -sqrt(3) / 2)))
  expect_equal(lengths(chebknots(c(x = 3, y = 4, z = 5))), c(x = 3L, y = 4L, z = 5L))
  expect_equal(range(chebknots(4, c(2, 3))[[1]]), 2.5 + c(-0.5, 0.5) * cos(pi / 8))
})

test_that("chebcoef() keeps dimensions and scales the constant term", {
  withr::local_seed(42)
  a <- array(rnorm(24), c(2, 3, 4))
  expect_equal(dim(chebcoef(a)), dim(a))
  expect_equal(chebcoef(a)[1], mean(a))
  expect_equal(chebcoef(a, dct = TRUE)[1], 8 * sum(a))
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
  skip_if_not(havefftw())
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
