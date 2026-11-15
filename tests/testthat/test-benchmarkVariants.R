test_that("benchmarkVariants returns metrics table", {
    ref <- system.file("extdata", "test_truth.vcf", package = "huGSVcalleR")
    qry <- system.file("extdata", "test_query.vcf", package = "huGSVcalleR")
    skip_if(ref == "" || qry == "", "test VCF not available")
    td <- tempdir()
    file.copy(ref, file.path(td, "ref.vcf"), overwrite = TRUE)
    file.copy(qry, file.path(td, "query.vcf"), overwrite = TRUE)
    metrics <- benchmarkVariants(
        vcfDir = NULL,
        ref_vcf = "ref.vcf",
        query_vcf = "query.vcf",
        workDir = td
    )
    expect_true(is.data.frame(metrics) || is.matrix(metrics))
    expect_true(all(c("Precision", "Recall", "F1") %in% colnames(metrics)))
    expect_true(all(c("SNPs", "INDELs", "ALL") %in% rownames(metrics)))
})
