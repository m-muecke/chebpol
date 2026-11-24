test_that(".onLoad() keeps a chebpol.threads option that is already set", {
  withr::local_options(chebpol.threads = 4L)
  withr::local_envvar(CHEBPOL_THREADS = "3")
  .onLoad()
  expect_equal(getOption("chebpol.threads"), 4L)
})

test_that(".onLoad() takes chebpol.threads from CHEBPOL_THREADS, or uses 1", {
  withr::local_options(chebpol.threads = NULL)
  withr::local_envvar(CHEBPOL_THREADS = "3")
  .onLoad()
  expect_equal(getOption("chebpol.threads"), 3L)

  options(chebpol.threads = NULL)
  withr::local_envvar(CHEBPOL_THREADS = NA)
  .onLoad()
  expect_equal(getOption("chebpol.threads"), 1L)
})
