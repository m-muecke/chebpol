test_that("interpolants accept integer points", {
  f <- function(x) sum(x^2)
  grid <- list(seq(0, 4), seq(0, 3))
  knots <- rbind(c(0, 4, 0, 4, 2, 1), c(0, 0, 3, 3, 1, 2))
  interpolants <- list(
    ipol(f, dims = c(5, 4), method = "chebyshev"),
    ipol(f, dims = c(5, 4), method = "uniform"),
    ipol(f, grid = grid, method = "general"),
    ipol(f, grid = grid, method = "multilinear"),
    ipol(f, grid = grid, method = "fh"),
    ipol(f, grid = grid, method = "stalker"),
    ipol(f, grid = grid, method = "hstalker"),
    ipol(f, grid = grid, method = "oldstalker"),
    ipol(f, knots = knots, method = "polyharmonic")
  )
  if (requireNamespace("geometry", quietly = TRUE)) {
    interpolants <- c(interpolants, ipol(f, knots = knots, method = "simplexlinear"))
  }
  for (ip in interpolants) {
    expect_equal(ip(c(1L, 1L)), ip(c(1, 1)))
    expect_equal(ip(matrix(1L, 2L, 3L)), ip(matrix(1, 2L, 3L)))
  }
})

test_that("ipol() accepts integer values", {
  values <- array(seq_len(20), c(5L, 4L))
  grid <- list(seq(0, 4), seq(0, 3))
  x <- c(1.5, 1.5)
  for (method in c("chebyshev", "uniform")) {
    expect_equal(ipol(values, method = method)(x / 4), ipol(values + 0, method = method)(x / 4))
  }
  for (method in c("general", "multilinear", "fh", "stalker", "hstalker", "oldstalker")) {
    expect_equal(ipol(values, grid = grid, method = method)(x), ipol(values + 0, grid = grid, method = method)(x))
  }
  knots <- rbind(c(0, 4, 0, 4, 2, 1), c(0, 0, 3, 3, 1, 2))
  values <- c(1L, 5L, 2L, 7L, 3L, 4L)
  for (method in c("polyharmonic", if (requireNamespace("geometry", quietly = TRUE)) "simplexlinear")) {
    expect_equal(ipol(values, knots = knots, method = method)(x), ipol(values + 0, knots = knots, method = method)(x))
  }
})

test_that("ipol() accepts logical values", {
  values <- c(TRUE, FALSE, TRUE, TRUE, FALSE)
  grid <- list(seq(0, 4))
  for (method in c("general", "multilinear", "fh", "stalker", "hstalker", "oldstalker")) {
    expect_equal(ipol(values, grid = grid, method = method)(1.5), ipol(values + 0, grid = grid, method = method)(1.5))
  }
})

test_that("ipol() accepts integer grids and knots", {
  f <- function(x) sum(x^2)
  grid <- list(0:4, 0:3)
  x <- c(1.5, 1.5)
  for (method in c(
    "general",
    "multilinear",
    "fh",
    "stalker",
    "hstalker",
    "oldstalker",
    if (requireNamespace("geometry", quietly = TRUE)) "simplexlinear"
  )) {
    expect_equal(
      ipol(f, grid = grid, method = method)(x),
      ipol(f, grid = lapply(grid, as.numeric), method = method)(x)
    )
  }
  knots <- rbind(c(0L, 4L, 0L, 4L, 2L, 1L), c(0L, 0L, 3L, 3L, 1L, 2L))
  for (method in c("polyharmonic", if (requireNamespace("geometry", quietly = TRUE)) "simplexlinear")) {
    expect_equal(ipol(f, knots = knots, method = method)(x), ipol(f, knots = knots + 0, method = method)(x))
  }
})

test_that("ipol() accepts a one-dimensional grid as a bare vector", {
  grid <- seq(0, 1, length.out = 10)
  x <- c(0.13, 0.5, 0.87)
  for (method in c("general", "multilinear", "fh", "stalker", "hstalker", "oldstalker", "polyharmonic")) {
    expect_equal(ipol(exp, grid = grid, method = method)(x), ipol(exp, grid = list(grid), method = method)(x))
  }
  suppressWarnings({
    uc <- ipol(exp(grid), grid = grid, method = "uniform")
    uc_list <- ipol(exp(grid), grid = list(grid), method = "uniform")
  })
  expect_equal(uc(x), uc_list(x))
})

test_that("ipol() errors on a grid without dimensions", {
  expect_snapshot(error = TRUE, {
    ipol(1, grid = list(), method = "multilinear")
    ipol(sum, grid = list(), method = "fh")
  })
})

test_that("ipol() errors on a grid with missing or infinite values", {
  expect_snapshot(error = TRUE, {
    ipol(1:3, grid = c(0, NA, 1), method = "multilinear")
    ipol(sum, grid = list(0:2, c(0, NaN, 1)), method = "fh")
    ipol(1:2, grid = c("a", "b"), method = "fh")
    ipol(1:3, grid = c(0, 1, Inf), method = "multilinear")
  })
})

test_that("ipol() errors when the function does not return a single number", {
  expect_snapshot(error = TRUE, {
    ipol(\(x) c(x, 1), dims = 3, method = "chebyshev")
    ipol(\(x) c(x, 1), dims = c(3, 4), method = "uniform")
    ipol(\(x) numeric(0), dims = 3, method = "chebyshev")
  })
})

test_that("ipol() errors when val is neither values nor a function", {
  grid <- list(1:3, 1:2)
  expect_snapshot(error = TRUE, {
    ipol("sum", grid = grid, method = "fh")
    ipol(as.list(1:6), grid = grid, method = "multilinear")
    ipol(data.frame(a = 1:6), dims = 6, method = "chebyshev")
  })
})
