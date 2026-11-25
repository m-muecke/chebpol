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

test_that("polyharmonic splines fit a single knot", {
  for (k in c(3, 2, -1)) {
    phs <- ipol(5, knots = matrix(0.5), k = k, method = "polyharmonic", nowarn = TRUE)
    expect_equal(phs(c(0.5, 0.2)), c(5, 5))
  }
})

test_that("polyharmonic splines forward extra arguments to the function", {
  withr::local_seed(42)
  f <- function(x, a) a * sum(x^2)
  knots <- matrix(runif(200, -1, 1), 2L)
  phs <- ipol(f, knots = knots, method = "polyharmonic", a = 3)
  expect_equal(phs(knots), apply(knots, 2L, f, a = 3))
})

test_that("polyharmonic splines accept a vector or matrix for normalize", {
  withr::local_seed(42)
  f <- function(x) 10 / (10 + sum(sqrt(x)))
  knots <- matrix(runif(400, 0, 20), 2L)
  for (normalize in list(c(1 / 20, 0), cbind(c(1 / 20, 1 / 10), c(0, 5)))) {
    phs <- ipol(f, knots = knots, method = "polyharmonic", normalize = normalize)
    expect_equal(phs(knots), apply(knots, 2L, f))
  }
})

test_that("normalize = TRUE equals the matrix it stands for", {
  withr::local_seed(42)
  f <- function(x) 10 / (10 + sum(sqrt(x)))
  knots <- matrix(runif(400, 0, 20), 2L)
  x <- matrix(runif(100, 1, 19), 2L)
  scale <- 1 / apply(knots, 1L, \(k) diff(range(k)))
  shift <- apply(knots, 1L, min)
  by_flag <- ipol(f, knots = knots, method = "polyharmonic", normalize = TRUE)
  by_matrix <- ipol(f, knots = knots, method = "polyharmonic", normalize = cbind(scale, shift))
  expect_equal(by_matrix(x), by_flag(x))
})

test_that("ipol() errors on a malformed normalize", {
  f <- function(x) sum(x)
  knots <- matrix(seq(0, 1, length.out = 20), 2L)
  expect_snapshot(error = TRUE, {
    ipol(f, knots = knots, method = "polyharmonic", normalize = c(1, 2, 3))
    ipol(f, knots = knots, method = "polyharmonic", normalize = matrix(1, 3L, 2L))
  })
})

test_that("ipol() errors when values and knots differ in number", {
  knots <- matrix(seq(0, 1, length.out = 20), 2L)
  expect_snapshot(error = TRUE, {
    ipol(1:3, knots = knots, method = "polyharmonic")
  })
})

test_that("the domain of a polyharmonic spline is the range of the knots", {
  withr::local_seed(42)
  knots <- matrix(runif(60, 5, 9), 2L)
  for (normalize in c(TRUE, FALSE)) {
    phs <- ipol(sum, knots = knots, method = "polyharmonic", normalize = normalize)
    expect_equal(as.matrix(attr(phs, "domain")), apply(knots, 1L, range), ignore_attr = TRUE)
  }
})

test_that("polyharmonic interpolation warns for the ipol() call when it falls back to least squares", {
  knots <- cbind(c(0, 0), c(1, 0), c(0, 1), c(0, 0))
  expect_snapshot(ph <- ipol(c(1, 2, 3, 4), knots = knots, method = "polyharmonic"))
  expect_no_warning(ipol(c(1, 2, 3, 4), knots = knots, method = "polyharmonic", nowarn = TRUE))
})

test_that("polyharmonic interpolation handles knots on a line or in a plane", {
  withr::local_seed(1)
  f <- function(x) x[1]^2
  knots <- rbind(runif(20, 0, 5), rep(2, 20))
  for (normalize in list(NA, FALSE)) {
    expect_no_warning(ph <- ipol(f, knots = knots, normalize = normalize, method = "polyharmonic"))
    expect_equal(ph(knots), apply(knots, 2L, f))
  }
  t <- runif(20)
  line <- rbind(t, 2 * t + 1)
  expect_no_warning(ph <- ipol(f, knots = line, method = "polyharmonic"))
  expect_equal(ph(line), apply(line, 2L, f))
  plane <- rbind(runif(30), runif(30))
  plane <- rbind(plane, plane[1L, ] - plane[2L, ])
  expect_no_warning(ph <- ipol(f, knots = plane, method = "polyharmonic"))
  expect_equal(ph(plane), apply(plane, 2L, f))
  single <- ipol(25, knots = matrix(c(5, 2), 2), method = "polyharmonic")
  expect_equal(single(c(1, 1)), 25)
})
