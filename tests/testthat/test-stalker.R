test_that("oldstalker interpolant errors on input of the wrong dimension", {
  f <- function(x) exp(-sum(x^2))
  grid <- list(seq(-1, 1, length.out = 9), seq(-1, 1, length.out = 9))
  st <- ipol(f, grid = grid, method = "oldstalker")
  expect_snapshot(st(matrix(0, 3, 4)), error = TRUE)
})
