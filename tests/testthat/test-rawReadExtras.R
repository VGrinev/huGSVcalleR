test_that("findAdapters bins the packaged adapter", {
    ex <- system.file("extdata", package = "huGSVcalleR")
    skip_if(ex == "", "test data not available")
    hits <- findAdapters(
        fastqDir = NULL,
        fastq = "test_R1.fastq",
        adapters = "adapters.txt",
        workDir = ex
    )
    expect_true("TestAdapter" %in% colnames(hits))
    expect_true("positions" %in% colnames(hits))
    expect_true(sum(hits$TestAdapter) > 0)
})

test_that("findForeignSeqs returns a numeric bin for the contaminant", {
    ex <- system.file("extdata", package = "huGSVcalleR")
    skip_if(ex == "", "test data not available")
    hits <- findForeignSeqs(
        fastqDir = NULL,
        fastq = "test_R1.fastq",
        contaminants = "contaminants.txt",
        workDir = ex
    )
    expect_true("TestContaminant" %in% colnames(hits))
    expect_true(all(hits$TestContaminant >= 0))
})

test_that("cutLeftAdapter removes the leading adapter from read 1", {
    ex <- system.file("extdata", package = "huGSVcalleR")
    skip_if(ex == "", "test data not available")
    trimmed <- cutLeftAdapter(
        fastqDir = NULL,
        fastq = "test_R1.fastq",
        adapters = "adapters.txt",
        workDir = ex
    )
    expect_s4_class(trimmed, "ShortReadQ")
    expect_equal(as.integer(Biostrings::width(trimmed)[1]), 37L)
    expect_equal(length(trimmed), 10L)
})

test_that("cutRigthAdapter keeps the read count", {
    ex <- system.file("extdata", package = "huGSVcalleR")
    skip_if(ex == "", "test data not available")
    trimmed <- cutRigthAdapter(
        fastqDir = NULL,
        fastq = "test_R1.fastq",
        adapters = "adapters.txt",
        workDir = ex
    )
    expect_s4_class(trimmed, "ShortReadQ")
    expect_equal(length(trimmed), 10L)
    expect_true(all(Biostrings::width(trimmed) <= 50))
})

test_that("reportQAResults writes a non-empty metrics report", {
    ex <- system.file("extdata", package = "huGSVcalleR")
    skip_if(ex == "", "test data not available")
    qa <- assessQRawReads(
        fastqDir = NULL,
        fastq = "test_R1.fastq",
        adapters = NULL,
        contaminants = NULL,
        workDir = ex
    )
    td <- tempfile("qa-report-")
    dir.create(td)
    reportQAResults(x = qa, output = "example", workDir = td)
    metrics <- file.path(td, "QA metrics of example.pdf")
    plots <- file.path(td, "QA plots of example.pdf")
    fasta <- file.path(td, "Overrepresented reads of example.fasta")
    expect_gt(file.info(metrics)$size, 1000)
    expect_gt(file.info(plots)$size, 1000)
    seqs <- Biostrings::readDNAStringSet(fasta)
    expect_gt(length(seqs), 0)
})
