test_that("cleanRawReads keeps reads intact by default", {
    fq <- system.file("extdata", "test_R1.fastq", package = "huGSVcalleR")
    skip_if(fq == "", "test FASTQ not available")
    td <- tempdir()
    file.copy(fq, file.path(td, "test_R1.fastq"), overwrite = TRUE)
    cleanRawReads(
        fastqDir  = NULL,
        fastq1    = "test_R1.fastq",
        postfix   = "plain",
        overwrite = TRUE,
        workDir   = td
    )
    out <- ShortRead::readFastq(file.path(td, "test_R1.plain.fastq.gz"))
    raw <- ShortRead::readFastq(fq)
    expect_equal(length(out), length(raw))
    expect_equal(as.integer(Biostrings::width(out)), as.integer(Biostrings::width(raw)))
})

test_that("cleanRawReads performs fixed trimming", {
    fq <- system.file("extdata", "test_R1.fastq", package = "huGSVcalleR")
    skip_if(fq == "", "test FASTQ not available")
    td <- tempdir()
    file.copy(fq, file.path(td, "test_R1.fastq"), overwrite = TRUE)
    cleanRawReads(
        fastqDir  = NULL,
        fastq1    = "test_R1.fastq",
        tr_start  = 5,
        tr_end    = 5,
        postfix   = "fixed",
        overwrite = TRUE,
        workDir   = td
    )
    expect_true(file.exists(file.path(td, "test_R1.fixed.fastq.gz")))
})

test_that("cleanRawReads performs dynamic trimming", {
    fq <- system.file("extdata", "test_R1.fastq", package = "huGSVcalleR")
    ph <- system.file("extdata", "phred_scores.txt", package = "huGSVcalleR")
    skip_if(fq == "" || ph == "", "test data not available")
    td <- tempdir()
    file.copy(fq, file.path(td, "test_R1.fastq"), overwrite = TRUE)
    file.copy(ph, file.path(td, "phred_scores.txt"), overwrite = TRUE)
    cleanRawReads(
        fastqDir  = NULL,
        fastq1    = "test_R1.fastq",
        tr_score  = 20,
        phred     = "phred_scores.txt",
        k         = 3,
        postfix   = "dyn",
        overwrite = TRUE,
        workDir   = td
    )
    expect_true(file.exists(file.path(td, "test_R1.dyn.fastq.gz")))
})

test_that("cleanRawReads trims adapters", {
    fq <- system.file("extdata", "test_R1.fastq", package = "huGSVcalleR")
    ad <- system.file("extdata", "adapters.txt", package = "huGSVcalleR")
    skip_if(fq == "" || ad == "", "test data not available")
    td <- tempdir()
    file.copy(fq, file.path(td, "test_R1.fastq"), overwrite = TRUE)
    file.copy(ad, file.path(td, "adapters.txt"), overwrite = TRUE)
    cleanRawReads(
        fastqDir  = NULL,
        fastq1    = "test_R1.fastq",
        adapters  = "adapters.txt",
        error     = 0.2,
        postfix   = "adapt",
        overwrite = TRUE,
        workDir   = td
    )
    expect_true(file.exists(file.path(td, "test_R1.adapt.fastq.gz")))
})

test_that("cleanRawReads handles paired-end", {
    fq1 <- system.file("extdata", "test_R1.fastq", package = "huGSVcalleR")
    fq2 <- system.file("extdata", "test_R2.fastq", package = "huGSVcalleR")
    skip_if(fq1 == "" || fq2 == "", "test FASTQ not available")
    td <- tempdir()
    file.copy(fq1, file.path(td, "test_R1.fastq"), overwrite = TRUE)
    file.copy(fq2, file.path(td, "test_R2.fastq"), overwrite = TRUE)
    cleanRawReads(
        fastqDir  = NULL,
        fastq1    = "test_R1.fastq",
        fastq2    = "test_R2.fastq",
        postfix   = "paired",
        overwrite = TRUE,
        workDir   = td
    )
    expect_true(file.exists(file.path(td, "test_R1.paired.fastq.gz")))
    expect_true(file.exists(file.path(td, "test_R2.paired.fastq.gz")))
})

test_that("cleanRawReads stops if output exists and overwrite = FALSE", {
    fq <- system.file("extdata", "test_R1.fastq", package = "huGSVcalleR")
    skip_if(fq == "", "test FASTQ not available")
    td <- tempdir()
    file.copy(fq, file.path(td, "test_R1.fastq"), overwrite = TRUE)
    writeLines("placeholder", file.path(td, "test_R1.conflict.fastq.gz"))
    expect_error(
        cleanRawReads(
            fastqDir  = NULL,
            fastq1    = "test_R1.fastq",
            postfix   = "conflict",
            overwrite = FALSE,
            workDir   = td
        ),
        "already exists"
    )
    file.remove(file.path(td, "test_R1.conflict.fastq.gz"))
})
