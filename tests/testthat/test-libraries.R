test_that("has_fftw() and has_gsl() return a single logical", {
  expect_type(has_fftw(), "logical")
  expect_length(has_fftw(), 1L)
  expect_type(has_gsl(), "logical")
  expect_length(has_gsl(), 1L)
})
