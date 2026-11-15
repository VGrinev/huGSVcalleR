test_that("extractSSsCoordsGTF returns fiveSSs and threeSSs", {
    gtf <- system.file("extdata", "test_annotation.gtf", package = "huGSVcalleR")
    skip_if(gtf == "", "test GTF not available")
    res <- extractSSsCoordsGTF(
        x = basename(gtf),
        gtfDir = NULL,
        workDir = dirname(gtf)
    )
    expect_type(res, "list")
    expect_named(res, c("fiveSSs", "threeSSs"))
    expect_s4_class(res$fiveSSs, "GRanges")
    expect_s4_class(res$threeSSs, "GRanges")
})
