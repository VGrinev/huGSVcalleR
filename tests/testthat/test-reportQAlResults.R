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
    metrics <- file.path(td, "Mapping quality metrics of test.pdf")
    plots <- file.path(td, "Mapping quality plots of test.pdf")
    expect_gt(file.info(metrics)$size, 500)
    expect_gt(file.info(plots)$size, 500)
})
