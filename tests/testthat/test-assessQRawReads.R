test_that("assessQRawReads returns list with all 15 slots", {
    fq <- system.file("extdata", "test_R1.fastq", package = "huGSVcalleR")
    skip_if(fq == "", "test FASTQ not available")
    res <- assessQRawReads(
        fastqDir     = NULL,
        fastq        = basename(fq),
        adapters     = NULL,
        contaminants = NULL,
        workDir      = dirname(fq)
    )
    expect_type(res, "list")
    expect_named(res, c(
        "SN", "FN", "TNRs", "TNBs", "RLs",
        "PERQs", "PCQRs", "BCs", "PCCBs",
        "GCs", "DSs", "FRs", "ORs",
        "ASs", "FSs"
    ))
    expect_equal(res$TNRs, 10)
    expect_equal(res$TNBs, 500)
    expect_true(is.na(res$ASs))
    expect_true(is.na(res$FSs))
})

test_that("assessQRawReads finds adapters", {
    fq <- system.file("extdata", "test_R1.fastq", package = "huGSVcalleR")
    ad <- system.file("extdata", "adapters.txt", package = "huGSVcalleR")
    skip_if(fq == "" || ad == "", "test data not available")
    res <- assessQRawReads(
        fastqDir     = NULL,
        fastq        = basename(fq),
        adapters     = "adapters.txt",
        contaminants = NULL,
        workDir      = dirname(fq)
    )
    expect_s3_class(res$ASs, "data.frame")
    expect_true(nrow(res$ASs) > 0)
})

test_that("assessQRawReads stops on invalid extension", {
    expect_error(
        assessQRawReads(
            fastqDir = NULL, fastq = "file.txt",
            workDir = tempdir()
        ),
        "Invalid file format"
    )
})
