missing_gatk <- "huGSVcalleR-no-gatk"

test_that("indexGenomeGATK writes a fasta index without GATK", {
    td <- tempfile("fai-")
    dir.create(file.path(td, "fasta"), recursive = TRUE)
    writeLines(
        c(">chr1", paste(rep("ACGT", 40), collapse = "")),
        file.path(td, "fasta", "tiny.fa")
    )
    indexGenomeGATK(
        fastaDir = "fasta",
        faFile = "tiny.fa",
        fai = TRUE,
        dict = FALSE,
        gatk_path = missing_gatk,
        workDir = td
    )
    fai <- read.delim(
        file.path(td, "fasta", "tiny.fa.fai"),
        header = FALSE
    )
    expect_equal(fai[1, 1], "chr1")
    expect_equal(fai[1, 2], 160)
})

test_that("GATK wrappers return status 127 for a missing executable", {
    td <- tempfile("gatk-wrap-")
    dir.create(td)
    expect_identical(suppressWarnings(groupReadsGATK(
        bamFile = "test.bam",
        SE = FALSE,
        gatk_path = missing_gatk,
        workDir = td
    )), 127L)
    expect_identical(suppressWarnings(markDuplicatesGATK(
        bamFile = "test.bam",
        SE = FALSE,
        gatk_path = missing_gatk,
        workDir = td
    )), 127L)
    expect_identical(suppressWarnings(indexVcfGATK(
        vcfFile = "test.vcf",
        gatk_path = missing_gatk,
        workDir = td
    )), 127L)
    expect_identical(suppressWarnings(indexGenomeGATK(
        faFile = "tiny.fa",
        fai = FALSE,
        dict = TRUE,
        gatk_path = missing_gatk,
        workDir = td
    )), 127L)
    expect_null(suppressWarnings(recalibrateBamGATK(
        bamFile = "test.bam",
        vcfFile = "test.vcf",
        faFile = "tiny.fa",
        createBQSR = TRUE,
        recalibrate = FALSE,
        nameBQSR = "unused",
        gatk_path = missing_gatk,
        workDir = td
    )))
    expect_false(file.exists(file.path(td, "test.bqsr.table")))
})
