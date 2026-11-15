test_that("extractSSsCoordsEEJs returns fiveSSs and threeSSs", {
  eej <- system.file("extdata", "test_eejs.txt", package = "huGSVcalleR")
  skip_if(eej == "", "test EEJ not available")
  res <- extractSSsCoordsEEJs(
    x = basename(eej),
    eejDir = NULL,
    thr = 5,
    workDir = dirname(eej)
  )
  expect_type(res, "list")
  expect_named(res, c("fiveSSs", "threeSSs"))
  expect_s4_class(res$fiveSSs, "GRanges")
  expect_s4_class(res$threeSSs, "GRanges")
  expect_equal(length(res$fiveSSs), 2)
  expect_equal(length(res$threeSSs), 2)
})
