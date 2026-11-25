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
  expect_equal(ips$multilinear(c(0.5, 1)), 1.25, tolerance = 0.1)
})

test_that("an interpolant created with a blend keeps its class and attributes", {
  grid <- list(seq(0, 1, length.out = 5), seq(0, 2, length.out = 4))
  ml <- ipol(function(x) sum(x), grid = grid, method = "multilinear", blend = "cubic")
  expect_s3_class(ml, "chebpol_interpolant")
  expect_equal(attr(ml, "method"), "multilinear")
  expect_equal(attr(ml, "arity"), 2)
  expect_equal(nrow(attr(ml, "domain")), 2)
})

test_that("printing an interpolant summarizes it", {
  grid <- list(seq(0, 1, length.out = 5), seq(0, 2, length.out = 4))
  expect_snapshot({
    ipol(function(x) sum(x), grid = grid, method = "multilinear", blend = "cubic")
    ipol(function(x) sum(x), dims = 6, method = "chebyshev")
    ipol(function(x) sum(x), grid = grid, method = "oldstalker", k = 1.5)
  })
})

test_that("printing a simplex linear interpolant omits the domain", {
  skip_if_not_installed("geometry")
  knots <- grid_points(list(c(0, 1, 2), c(0, 1)))
  expect_snapshot(ipol(function(x) sum(x), knots = knots, method = "simplexlinear"))
})
