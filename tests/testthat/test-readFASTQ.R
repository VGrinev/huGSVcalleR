test_that("readFASTQ returns list with reads and QScores", {
    fq <- system.file("extdata", "test_R1.fastq", package = "huGSVcalleR")
    skip_if(fq == "", "test FASTQ not available")
    res <- readFASTQ(
        fastqDir = NULL,
        fastq = basename(fq),
        workDir = dirname(fq)
    )
    expect_type(res, "list")
    expect_named(res, c("reads", "QScores"))
    expect_s4_class(res$reads, "DNAStringSet")
    expect_s4_class(res$QScores, "PhredQuality")
    expect_equal(length(res$reads), 10)
    expect_equal(width(res$reads), rep(50, 10))
})

test_that("readFASTQ respects n argument", {
    fq <- system.file("extdata", "test_R1.fastq", package = "huGSVcalleR")
    skip_if(fq == "", "test FASTQ not available")
    res <- readFASTQ(
        fastqDir = NULL,
        fastq = basename(fq),
        n = 3,
        workDir = dirname(fq)
    )
    expect_equal(length(res$reads), 3)
})

test_that("readFASTQ stops on invalid file extension", {
    expect_error(
        readFASTQ(
            fastqDir = NULL, fastq = "file.txt",
            workDir = tempdir()
        ),
        "Invalid file format"
    )
})
