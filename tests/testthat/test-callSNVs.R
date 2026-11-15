test_that("callSNVs with fisher writes VCF", {
    fr <- system.file("extdata", "test_pileup.csv", package = "huGSVcalleR")
    skip_if(fr == "", "test pileup not available")
    td <- tempdir()
    callSNVs(
        fr_data   = basename(fr),
        criteria  = "fisher",
        vcfDir    = NULL,
        outputVCF = "test_fisher",
        workDir   = dirname(fr)
    )
    expect_true(file.exists(
        file.path(dirname(fr), "test_fisher.fisher.vcf")
    ))
})

test_that("callSNVs with binomial writes VCF", {
    fr <- system.file("extdata", "test_pileup.csv", package = "huGSVcalleR")
    skip_if(fr == "", "test pileup not available")
    callSNVs(
        fr_data   = basename(fr),
        criteria  = "binomial",
        outputVCF = "test_binom",
        workDir   = dirname(fr)
    )
    expect_true(file.exists(
        file.path(dirname(fr), "test_binom.binomial.vcf")
    ))
})

test_that("callSNVs with counting writes VCF", {
    fr <- system.file("extdata", "test_pileup.csv", package = "huGSVcalleR")
    skip_if(fr == "", "test pileup not available")
    callSNVs(
        fr_data   = basename(fr),
        criteria  = "counting",
        outputVCF = "test_count",
        workDir   = dirname(fr)
    )
    expect_true(file.exists(
        file.path(dirname(fr), "test_count.counting.vcf")
    ))
})

test_that("callSNVs with entropy writes VCF", {
    fr <- system.file("extdata", "test_pileup.csv", package = "huGSVcalleR")
    skip_if(fr == "", "test pileup not available")
    callSNVs(
        fr_data   = basename(fr),
        criteria  = "entropy",
        outputVCF = "test_ent",
        workDir   = dirname(fr)
    )
    expect_true(file.exists(
        file.path(dirname(fr), "test_ent.entropy.vcf")
    ))
})

test_that("callSNVs with poisson writes VCF", {
    fr <- system.file("extdata", "test_pileup.csv", package = "huGSVcalleR")
    skip_if(fr == "", "test pileup not available")
    callSNVs(
        fr_data   = basename(fr),
        criteria  = "poisson",
        outputVCF = "test_pois",
        workDir   = dirname(fr)
    )
    expect_true(file.exists(
        file.path(dirname(fr), "test_pois.poisson.vcf")
    ))
})

test_that("callSNVs with ml rpart writes VCF", {
    fr <- system.file("extdata", "test_pileup.csv", package = "huGSVcalleR")
    skip_if(fr == "", "test pileup not available")
    skip_if_not_installed("rpart")
    callSNVs(
        fr_data   = basename(fr),
        criteria  = "ml",
        model     = "rpart",
        pvalue    = TRUE,
        outputVCF = "test_ml_rpart",
        workDir   = dirname(fr)
    )
    expect_true(file.exists(
        file.path(dirname(fr), "test_ml_rpart.ML_Prediction.vcf")
    ))
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

test_that("callSNVs stops on invalid criteria", {
    expect_error(
        callSNVs(fr_data = "x.csv", criteria = "unknown", outputVCF = "x"),
        "'arg' should be one of"
    )
})
