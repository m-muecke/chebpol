test_that("polyharmonic splines reproduce values at the knots", {
  withr::local_seed(42)
  f <- function(x) 10 / (10 + sum(sqrt(x)))
  knots <- matrix(runif(6000, 0, 20), 6)
  phs <- ipol(f, knots = knots, k = 3, method = "polyharmonic")
  expect_equal(phs(knots), apply(knots, 2, f))
})

test_that("polyharmonic splines approximate a smooth function", {
  withr::local_seed(42)
  f <- function(x) 10 / (10 + sum(sqrt(x)))
  knots <- matrix(runif(6000, 0, 20), 6)
  a <- matrix(runif(600, 0, 20), 6)
  for (k in c(3, 5, -1)) {
    phs <- ipol(f, knots = knots, k = k, method = "polyharmonic", normalize = TRUE)
    expect_equal(phs(a), apply(a, 2, f), tolerance = 0.1)
  }
})

test_that("polyharmonic splines forward extra arguments to the function", {
  withr::local_seed(42)
  f <- function(x, a) a * sum(x^2)
  knots <- matrix(runif(200, -1, 1), 2L)
  phs <- ipol(f, knots = knots, method = "polyharmonic", a = 3)
  expect_equal(phs(knots), apply(knots, 2L, f, a = 3))
})
