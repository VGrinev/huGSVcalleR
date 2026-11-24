test_that("filterVariants filters by depth", {
    vcf <- system.file("extdata", "test_truth.vcf", package = "huGSVcalleR")
    skip_if(vcf == "", "test VCF not available")
    td <- tempdir()
    file.copy(vcf, file.path(td, "test.vcf"), overwrite = TRUE)
    res <- filterVariants(
        vcfDir = NULL, vcfFile = "test.vcf",
        depth = 5, workDir = td
    )
    expect_s4_class(res, "CollapsedVCF")
    expect_equal(length(res), 2L)
})

test_that("filterVariants drops sites below the depth threshold", {
    vcf <- system.file("extdata", "test_truth.vcf", package = "huGSVcalleR")
    skip_if(vcf == "", "test VCF not available")
    td <- tempdir()
    file.copy(vcf, file.path(td, "test.vcf"), overwrite = TRUE)
    res <- filterVariants(
        vcfDir = NULL, vcfFile = "test.vcf",
        depth = 100, workDir = td
    )
    expect_equal(length(res), 0L)
})

test_that("filterVariants filters by region", {
    vcf <- system.file("extdata", "test_truth.vcf", package = "huGSVcalleR")
    reg <- system.file("extdata", "test_regions.txt", package = "huGSVcalleR")
    skip_if(vcf == "" || reg == "", "test data not available")
    td <- tempdir()
    file.copy(vcf, file.path(td, "test.vcf"), overwrite = TRUE)
    file.copy(reg, file.path(td, "regions.txt"), overwrite = TRUE)
    res <- filterVariants(
        vcfDir = NULL, vcfFile = "test.vcf",
        gr = "regions.txt", workDir = td
    )
    expect_s4_class(res, "CollapsedVCF")
})

test_that("filterVariants sets FILTER to PASS", {
    vcf <- system.file("extdata", "test_truth.vcf", package = "huGSVcalleR")
    skip_if(vcf == "", "test VCF not available")
    td <- tempdir()
    file.copy(vcf, file.path(td, "test.vcf"), overwrite = TRUE)
    res <- filterVariants(
        vcfDir = NULL, vcfFile = "test.vcf",
        depth = 1, workDir = td
    )
    filters <- as.character(fixed(res)$FILTER)
    expect_true(all(filters == "PASS" | filters == "."))
})

test_that("filterVariants writes filtered VCF", {
    vcf <- system.file("extdata", "test_truth.vcf", package = "huGSVcalleR")
    skip_if(vcf == "", "test VCF not available")
    td <- tempdir()
    file.copy(vcf, file.path(td, "test.vcf"), overwrite = TRUE)
    filterVariants(
        vcfDir = NULL, vcfFile = "test.vcf",
        depth = 1, workDir = td
    )
    expect_true(file.exists(file.path(td, "test.filtered.vcf")))
})
