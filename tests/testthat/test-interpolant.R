test_that("every method returns a chebpol_interpolant that records its method", {
  f <- function(x) sum(x^2)
  grid <- list(seq(0, 1, length.out = 5), seq(0, 2, length.out = 4))
  knots <- grid_points(grid)
  ips <- list(
    chebyshev = ipol(f, dims = c(5, 4), method = "chebyshev"),
    uniform = ipol(f, dims = c(5, 4), method = "uniform"),
    general = ipol(f, grid = grid, method = "general"),
    fh = ipol(f, grid = grid, method = "fh"),
    multilinear = ipol(f, grid = grid, method = "multilinear"),
    stalker = ipol(f, grid = grid, method = "stalker"),
    hstalker = ipol(f, grid = grid, method = "hstalker"),
    oldstalker = ipol(f, grid = grid, method = "oldstalker"),
    polyharmonic = ipol(f, knots = knots, method = "polyharmonic")
  )
  for (ip in ips) {
    expect_s3_class(ip, "chebpol_interpolant")
  }
  expect_equal(vapply(ips, \(ip) attr(ip, "method"), ""), names(ips), ignore_attr = TRUE)
})

test_that("an interpolant created with a blend keeps its class and attributes", {
  grid <- list(seq(0, 1, length.out = 5), seq(0, 2, length.out = 4))
  ml <- ipol(function(x) sum(x), grid = grid, method = "multilinear", blend = "cubic")
  expect_s3_class(ml, "chebpol_interpolant")
  expect_equal(attr(ml, "method"), "multilinear")
  expect_equal(attr(ml, "arity"), 2)
  expect_equal(nrow(attr(ml, "domain")), 2)
})

test_that("domain columns are named after the variables, or X1, X2, ... without names", {
  f <- function(x) sum(x)
  grid <- list(seq(0, 1, length.out = 4), seq(0, 2, length.out = 3))
  knots <- grid_points(grid)
  domain_names <- function(ip) names(attr(ip, "domain"))
  expect_equal(domain_names(ipol(f, grid = grid, method = "fh")), c("X1", "X2"))
  expect_equal(domain_names(ipol(f, grid = list(a = grid[[1]], grid[[2]]), method = "multilinear")), c("a", "X2"))
  expect_equal(domain_names(ipol(f, grid = grid[[1]], method = "stalker")), "X1")
  expect_equal(domain_names(ipol(f, dims = c(4, 3), method = "chebyshev")), c("X1", "X2"))
  expect_equal(domain_names(ipol(f, dims = c(a = 4, b = 3), method = "chebyshev")), c("a", "b"))
  expect_equal(domain_names(ipol(1:12, dims = c(a = 4, b = 3), method = "chebyshev")), c("a", "b"))
  iv <- list(a = c(0, 1), b = c(0, 2))
  expect_equal(domain_names(ipol(f, dims = c(4, 3), intervals = iv, method = "chebyshev")), c("a", "b"))
  expect_equal(domain_names(ipol(f, dims = c(a = 4, b = 3), method = "uniform")), c("a", "b"))
  expect_equal(domain_names(ipol(f, dims = c(a = 4, b = 3), intervals = unname(iv), method = "uniform")), c("a", "b"))
  expect_equal(domain_names(ipol(f, knots = unname(knots), method = "polyharmonic")), c("X1", "X2"))
  expect_equal(domain_names(ipol(f, grid = grid, method = "polyharmonic")), c("X1", "X2"))
  expect_equal(domain_names(ipol(f, grid = list(a = grid[[1]], b = grid[[2]]), method = "polyharmonic")), c("a", "b"))
  expect_equal(domain_names(ipol(\(x) x, knots = matrix(c(0, 0.5, 1), 1), method = "polyharmonic")), "X1")
  rownames(knots) <- c("a", "b")
  expect_equal(domain_names(ipol(f, knots = knots, method = "polyharmonic")), c("a", "b"))
})

test_that("domain keeps the names of the variables as given", {
  grid <- list(`my var` = c(0, 1, 2), `my var` = c(0, 3))
  ip <- ipol(function(x) sum(x), grid = grid, method = "fh")
  expect_named(attr(ip, "domain"), c("my var", "my var"))
})

test_that("domain has the rows lower and upper, also for a reversed interval", {
  ip <- ipol(function(x) x, dims = 4, intervals = c(1, 0), method = "chebyshev")
  expect_equal(attr(ip, "domain"), data.frame(X1 = c(0, 1), row.names = c("lower", "upper")))
})

test_that("a simplex linear interpolant records the range of its knots", {
  skip_if_not_installed("geometry")
  sl <- ipol(function(x) sum(x), grid = list(c(0, 1, 2), c(0, 3)), method = "simplexlinear")
  expect_equal(attr(sl, "domain"), data.frame(X1 = c(0, 2), X2 = c(0, 3), row.names = c("lower", "upper")))
})

test_that("an interpolant does not keep the extra arguments of the function", {
  f <- function(x, data) sum(x) + data[1]
  data <- rep(1, 1e6)
  grid <- list(seq(0, 1, length.out = 5), seq(0, 2, length.out = 4))
  ips <- list(
    ipol(f, dims = c(5, 4), method = "chebyshev", data = data),
    ipol(f, dims = c(5, 4), method = "uniform", data = data),
    ipol(f, grid = grid, method = "general", data = data),
    ipol(f, grid = list(c(0, 1, 2), 5), method = "general", data = data),
    ipol(f, grid = grid, method = "fh", data = data),
    ipol(f, grid = grid, method = "multilinear", data = data),
    ipol(f, grid = grid, method = "stalker", data = data),
    ipol(f, grid = grid, method = "hstalker", data = data),
    ipol(f, grid = grid, method = "oldstalker", data = data),
    ipol(f, grid = grid, method = "polyharmonic", data = data)
  )
  for (ip in ips) {
    expect_lt(length(serialize(ip, NULL)), length(serialize(data, NULL)) / 10)
  }
})

test_that("a simplex linear interpolant does not keep the extra arguments of the function", {
  skip_if_not_installed("geometry")
  f <- function(x, data) sum(x) + data[1]
  data <- rep(1, 1e6)
  sl <- ipol(f, grid = list(c(0, 1, 2), c(0, 1)), method = "simplexlinear", data = data)
  expect_lt(length(serialize(sl, NULL)), length(serialize(data, NULL)) / 10)
})

test_that("printing an interpolant summarizes it", {
  grid <- list(seq(0, 1, length.out = 5), seq(0, 2, length.out = 4))
  expect_snapshot({
    ipol(function(x) sum(x), grid = grid, method = "multilinear", blend = "cubic")
    ipol(function(x) sum(x), dims = 6, method = "chebyshev")
    ipol(function(x) sum(x), grid = grid, method = "oldstalker", k = 1.5)
  })
})

test_that("printing a simplex linear interpolant shows the range of the knots", {
  skip_if_not_installed("geometry")
  knots <- grid_points(list(c(0, 1, 2), c(0, 1)))
  expect_snapshot(ipol(function(x) sum(x), knots = knots, method = "simplexlinear"))
})
