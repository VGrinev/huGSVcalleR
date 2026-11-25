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

test_that("per-cycle base composition keeps every 5-cycle window", {
    fq <- system.file("extdata", "test_R1.fastq", package = "huGSVcalleR")
    skip_if(fq == "", "test FASTQ not available")
    res <- assessQRawReads(
        fastqDir = NULL,
        fastq = basename(fq),
        adapters = NULL,
        contaminants = NULL,
        workDir = dirname(fq)
    )
    expect_equal(nrow(res$PCCBs), 10)
    expect_false(anyNA(res$PCCBs[, -1]))
    expect_equal(rowSums(res$PCCBs[, -1]), rep(1, 10), tolerance = 1e-8)
})

test_that("per-cycle summaries keep unequal cycle lengths", {
    quality <- data.frame(
        Cycle = c(1, 2, 3),
        Quality = c("I", "I", "I"),
        Score = c(10, 20, 30),
        Count = c(4, 4, 2),
        lane = "x"
    )
    base_call <- data.frame(
        Cycle = c(1, 2, 3),
        Base = c("A", "C", "G"),
        Count = c(4, 4, 2),
        lane = "x"
    )
    qa <- list(perCycle = list(quality = quality, baseCall = base_call))
    expect_silent(scores <- huGSVcalleR:::.expandPerCycleQuality(qa))
    expect_equal(lengths(scores), c("1" = 4, "2" = 4, "3" = 2))
    binned <- huGSVcalleR:::.binPerCycleQuality(scores)
    expected <- quantile(
        c(rep(10, 4), rep(20, 4), rep(30, 2)),
        probs = c(0.1, 0.25, 0.5, 0.75, 0.9)
    )
    expect_equal(as.numeric(binned[[1]][-1]), as.numeric(expected))
    expect_silent(bases <- huGSVcalleR:::.expandPerCycleBases(qa))
    expect_equal(lengths(bases), c("1" = 4, "2" = 4, "3" = 2))
    counted <- huGSVcalleR:::.binPerCycleBases(bases)
    expect_equal(as.numeric(counted$A[1]), 4)
    expect_equal(as.numeric(counted$C[1]), 4)
    expect_equal(as.numeric(counted$G[1]), 2)
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
