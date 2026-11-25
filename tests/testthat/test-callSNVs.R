read_vcf_body <- function(path) {
    lines <- readLines(path)
    lines <- lines[!startsWith(lines, "#")]
    read.delim(
        text = paste(lines, collapse = "\n"),
        header = FALSE,
        sep = "\t",
        col.names = c(
            "CHROM", "POS", "ID", "REF", "ALT",
            "QUAL", "FILTER", "INFO"
        )
    )
}

run_table_caller <- function(criteria, output, model = "rpart",
    pvalue = FALSE) {
    fr <- system.file("extdata", "test_pileup.csv", package = "huGSVcalleR")
    td <- tempfile(paste0("call-", criteria, "-"))
    dir.create(td)
    file.copy(fr, file.path(td, "pileup.csv"), overwrite = TRUE)
    ex_before <- list.files(system.file(
        "extdata",
        package = "huGSVcalleR"
    ))
    callSNVs(
        fr_data = "pileup.csv",
        criteria = criteria,
        model = model,
        pvalue = pvalue,
        vcfDir = NULL,
        outputVCF = output,
        workDir = td
    )
    ex_after <- list.files(system.file(
        "extdata",
        package = "huGSVcalleR"
    ))
    expect_identical(ex_after, ex_before)
    td
}

test_that("callSNVs fisher records the three pileup sites", {
    fr <- system.file("extdata", "test_pileup.csv", package = "huGSVcalleR")
    skip_if(fr == "", "test pileup not available")
    td <- run_table_caller("fisher", "sites")
    body <- read_vcf_body(file.path(td, "sites.fisher.vcf"))
    expect_equal(body$CHROM, rep("chr1", 3))
    expect_equal(body$POS, c(100, 200, 300))
    expect_equal(body$REF, c("A", "C", "G"))
    expect_equal(body$ALT, c("G", "G", "T"))
    expect_true(all(is.finite(body$QUAL)))
    expect_true(all(body$QUAL >= 0))
    strongest <- body$QUAL[body$POS == 200]
    expect_true(strongest >= max(body$QUAL[body$POS != 200]))
})

test_that("callSNVs binomial QUAL is a positive Phred score", {
    fr <- system.file("extdata", "test_pileup.csv", package = "huGSVcalleR")
    skip_if(fr == "", "test pileup not available")
    td <- run_table_caller("binomial", "sites")
    body <- read_vcf_body(file.path(td, "sites.binomial.vcf"))
    expect_equal(nrow(body), 3)
    expect_true(all(is.finite(body$QUAL)))
    expect_true(all(body$QUAL > 0))
})

test_that("callSNVs counting stores the reference-read ratio", {
    fr <- system.file("extdata", "test_pileup.csv", package = "huGSVcalleR")
    skip_if(fr == "", "test pileup not available")
    td <- run_table_caller("counting", "sites")
    body <- read_vcf_body(file.path(td, "sites.counting.vcf"))
    row <- body[body$POS == 100, ]
    expect_equal(row$REF, "A")
    expect_equal(row$QUAL, 0.8)
})

test_that("callSNVs entropy and poisson write one row per pileup site", {
    fr <- system.file("extdata", "test_pileup.csv", package = "huGSVcalleR")
    skip_if(fr == "", "test pileup not available")
    for (criteria in c("entropy", "poisson")) {
        td <- run_table_caller(criteria, "sites")
        body <- read_vcf_body(file.path(td, paste0("sites.", criteria, ".vcf")))
        expect_equal(body$POS, c(100, 200, 300))
        expect_true(all(body$QUAL >= 0))
    }
})

test_that("callSNVs rpart writes probabilities between 0 and 1", {
    fr <- system.file("extdata", "test_pileup.csv", package = "huGSVcalleR")
    skip_if(fr == "", "test pileup not available")
    skip_if_not_installed("rpart")
    td <- run_table_caller("ml", "sites", pvalue = TRUE)
    body <- read_vcf_body(file.path(td, "sites.ML_Prediction.vcf"))
    expect_true(nrow(body) > 0)
    expect_true(all(body$QUAL >= 0 & body$QUAL <= 1))
})

test_that("callSNVs stops without fr_data for table-based criteria", {
    expect_error(
        callSNVs(criteria = "fisher", outputVCF = "x"),
        "fr_data is required"
    )
})

test_that("callSNVs stops without bamFile for BAM-based criteria", {
    expect_error(
        callSNVs(criteria = "gatk", outputVCF = "x"),
        "bamFile is required"
    )
})

test_that("callSNVs stops when the GATK executable is missing", {
    td <- tempfile("gatk-call-")
    dir.create(td)
    expect_error(
        callSNVs(
            criteria = "gatk",
            bamFile = "test.bam",
            faFile = "tiny.fa",
            outputVCF = "sites",
            gatk_path = "huGSVcalleR-no-gatk",
            workDir = td
        ),
        "GATK"
    )
})

test_that("callSNVs stops on invalid criteria", {
    expect_error(
        callSNVs(fr_data = "x.csv", criteria = "unknown", outputVCF = "x"),
        "'arg' should be one of"
    )
})
