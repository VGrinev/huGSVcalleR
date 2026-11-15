test_that("sortBamFile sorts and indexes BAM", {
    bam <- system.file("extdata", "test.bam", package = "huGSVcalleR")
    bai <- system.file("extdata", "test.bam.bai", package = "huGSVcalleR")
    skip_if(bam == "" || bai == "", "test BAM not available")
    td <- tempdir()
    file.copy(bam, file.path(td, "test.bam"), overwrite = TRUE)
    file.copy(bai, file.path(td, "test.bam.bai"), overwrite = TRUE)
    sortBamFile(
        bamDir  = NULL,
        bamFile = "test.bam",
        byQname = FALSE,
        threads = 1,
        workDir = td
    )
    expect_true(file.exists(file.path(td, "test.sorted.bam")))
    expect_true(file.exists(file.path(td, "test.sorted.bam.bai")))
})
