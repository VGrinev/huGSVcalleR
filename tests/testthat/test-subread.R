write_tiny_fasta <- function(dir) {
    dir.create(file.path(dir, "fasta"), showWarnings = FALSE)
    writeLines(
        c(">chr1", paste(rep("ACGT", 40), collapse = "")),
        file.path(dir, "fasta", "tiny.fa")
    )
}

test_that("buildIndexSubread writes an index for a tiny fasta", {
    td <- tempfile("index-")
    dir.create(td)
    write_tiny_fasta(td)
    buildIndexSubread(
        ref_genome = "index",
        ref_fasta = "fasta",
        index = "tiny",
        fa = "tiny.fa",
        memory = 200,
        workDir = td
    )
    index_files <- list.files(file.path(td, "index"), pattern = "^tiny")
    expect_gt(length(index_files), 0)
})

test_that("alignDNASubread writes a BAM with every input read", {
    ex <- system.file("extdata", "test_R1.fastq", package = "huGSVcalleR")
    skip_if(ex == "", "test FASTQ not available")
    td <- tempfile("align-")
    dir.create(td)
    write_tiny_fasta(td)
    buildIndexSubread(
        ref_genome = "index",
        ref_fasta = "fasta",
        index = "tiny",
        fa = "tiny.fa",
        memory = 200,
        workDir = td
    )
    reads <- readLines(ex)
    gz <- gzfile(file.path(td, "test_R1.fastq.gz"), "w")
    writeLines(reads, gz)
    close(gz)
    dir.create(file.path(td, "bam"))
    alignDNASubread(
        genome = "index/tiny",
        fastqDir = NULL,
        fastq1 = "test_R1.fastq.gz",
        bamDir = "bam",
        bamFile = "example.bam",
        threads = 1,
        workDir = td
    )
    counted <- Rsamtools::countBam(file.path(td, "bam", "example.bam"))
    expect_equal(counted$records, 10L)
})
