test_that("havefftw() and havealglib() are deprecated", {
  expect_snapshot({
    . <- havefftw()
    . <- havealglib()
  })
})

test_that("havefftw() and havealglib() return what has_fftw() and has_alglib() do", {
  suppressWarnings({
    expect_identical(havefftw(), has_fftw())
    expect_identical(havealglib(), has_alglib())
  })
})
