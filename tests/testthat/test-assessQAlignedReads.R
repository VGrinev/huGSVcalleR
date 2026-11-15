test_that("assessQAlignedReads returns list with six slots", {
    bam <- system.file("extdata", "test.bam", package = "huGSVcalleR")
    bai <- system.file("extdata", "test.bam.bai", package = "huGSVcalleR")
    skip_if(bam == "" || bai == "", "test BAM not available")
    td <- tempdir()
    file.copy(bam, file.path(td, "test.bam"), overwrite = TRUE)
    file.copy(bai, file.path(td, "test.bam.bai"), overwrite = TRUE)
    qal <- assessQAlignedReads(
        bamDir  = NULL,
        bamFile = "test.bam",
        workDir = td
    )
    expect_type(qal, "list")
    expect_named(qal, c("SN", "FN", "TNRs", "MAPQs", "TLENs", "FLAGs"))
    expect_equal(qal$SN, "test")
    expect_equal(qal$FN, "test.bam")
    expect_equal(qal$TNRs, 4)
    expect_true(length(qal$MAPQs) > 0)
    expect_true("isPaired" %in% names(qal$FLAGs) ||
        "isUnmappedQuery" %in% names(qal$FLAGs))
})
