test_that("uniform grid interpolation reproduces values at the grid points", {
  f <- function(x) 1 / (1 + 25 * x^2)
  grid <- seq(-1, 1, length.out = 15)
  uc <- ipol(f(grid), method = "uniform")
  expect_equal(uc(grid), f(grid))
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
