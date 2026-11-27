test_that("interpolants accept a point, a matrix, or a vector of several points", {
  grid <- list(c(0, 1, 2), c(0, 1))
  ml <- ipol(evalongrid(sum, grid = grid), grid = grid, method = "multilinear")
  m <- rbind(c(0.5, 1.5, 2), c(0.2, 0.4, 1))
  expect_equal(ml(c(0.5, 0.2)), 0.7)
  expect_equal(ml(m), colSums(m))
  expect_equal(ml(c(1L, 1L)), 2)
  expect_snapshot(out <- ml(c(m)))
  expect_equal(out, colSums(m))
  ml1 <- ipol(c(0, 2, 4), grid = list(c(0, 1, 2)), method = "multilinear")
  expect_equal(ml1(c(0.5, 1.5)), c(1, 3))
})

test_that("interpolants error when the input length does not fit the number of variables", {
  grid <- list(c(0, 1), c(0, 1), c(0, 1))
  ml <- ipol(evalongrid(sum, grid = grid), grid = grid, method = "multilinear")
  gen <- ipol(evalongrid(sum, grid = grid), grid = grid, method = "general")
  uni <- ipol(sum, dims = c(2, 2, 2), method = "uniform")
  expect_snapshot(error = TRUE, {
    ml(c(0.1, 0.2))
    gen(c(0.1, 0.2))
    gen(matrix(0, 2, 4))
    uni(matrix(0, 2, 4))
  })
})

test_that("ipol() sets the default blend of the interpolant", {
  grid <- list(c(0, 1, 2))
  ml <- ipol(c(0, 1, 0), grid = grid, method = "multilinear", blend = "cubic")
  expect_equal(eval(formals(ml)$blend), c("cubic", "linear", "sigmoid", "parodic", "square"))
  for (method in c("multilinear", "stalker", "hstalker", "oldstalker")) {
    ip <- ipol(c(0, 1, 0), grid = grid, method = method, blend = "sigmoid")
    expect_equal(eval(formals(ip)$blend)[1], "sigmoid")
    expect_equal(ip(0.25), ip(0.25, blend = "sigmoid"))
  }
})

test_that("ipol() sets the default blend of a simplex linear interpolant", {
  skip_if_not_installed("geometry")
  knots <- grid_points(list(c(0, 1, 2), c(0, 1)))
  sl <- ipol(c(0, 1, 0, 1, 2, 1), knots = knots, method = "simplexlinear", blend = "sigmoid")
  expect_equal(eval(formals(sl)$blend)[1], "sigmoid")
  expect_equal(sl(c(0.3, 0.6)), sl(c(0.3, 0.6), blend = "sigmoid"))
})

test_that("crate() keeps only the named objects in the environment of a function", {
  a <- 1
  unused <- 2
  f <- crate(function(x) x + a, a)
  expect_equal(f(1), 2)
  expect_equal(ls(environment(f)), "a")
})
