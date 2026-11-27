test_that("havefftw() and havealglib() are deprecated", {
  expect_snapshot({
    . <- havefftw()
    . <- havealglib()
  })
})

test_that("havefftw() returns what has_fftw() does and havealglib() returns FALSE", {
  suppressWarnings({
    expect_identical(havefftw(), has_fftw())
    expect_identical(havealglib(), FALSE)
  })
})
