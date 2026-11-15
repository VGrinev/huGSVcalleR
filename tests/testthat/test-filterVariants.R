test_that("filterVariants returns filtered CollapsedVCF", {
    vcf <- system.file("extdata", "test_truth.vcf", package = "huGSVcalleR")
    skip_if(vcf == "", "test VCF not available")
    td <- tempdir()
    file.copy(vcf, file.path(td, "test_truth.vcf"), overwrite = TRUE)
    res <- filterVariants(
        vcfDir = NULL, vcfFile = "test_truth.vcf",
        depth = 5, workDir = td
    )
    expect_s4_class(res, "CollapsedVCF")
    expect_true(file.exists(file.path(td, "test_truth.filtered.vcf")))
})

test_that("filterVariants applies QUAL threshold", {
    vcf <- system.file("extdata", "test_truth.vcf", package = "huGSVcalleR")
    skip_if(vcf == "", "test VCF not available")
    td <- tempdir()
    file.copy(vcf, file.path(td, "test_truth.vcf"), overwrite = TRUE)
    res <- filterVariants(
        vcfDir = NULL, vcfFile = "test_truth.vcf",
        score = 1, workDir = td
    )
    expect_s4_class(res, "CollapsedVCF")
})
