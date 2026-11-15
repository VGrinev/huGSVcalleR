test_that("reportQAlResults writes two PDF files", {
    bam <- system.file("extdata", "test.bam", package = "huGSVcalleR")
    bai <- system.file("extdata", "test.bam.bai", package = "huGSVcalleR")
    skip_if(bam == "" || bai == "", "test BAM not available")
    td <- tempdir()
    file.copy(bam, file.path(td, "test.bam"), overwrite = TRUE)
    file.copy(bai, file.path(td, "test.bam.bai"), overwrite = TRUE)
    qal <- assessQAlignedReads(
        bamDir = NULL, bamFile = "test.bam",
        workDir = td
    )
    reportQAlResults(x = qal, output = NULL, workDir = td)
    expect_true(file.exists(
        file.path(td, "Mapping quality metrics of test.pdf")
    ))
    expect_true(file.exists(
        file.path(td, "Mapping quality plots of test.pdf")
    ))
})
