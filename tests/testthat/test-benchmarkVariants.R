test_that("benchmarkVariants returns expected metrics", {
  ref <- system.file("extdata", "test_truth.vcf", package = "huGSVcalleR")
  qry <- system.file("extdata", "test_query.vcf", package = "huGSVcalleR")
  skip_if(ref == "" || qry == "", "test VCF not available")
  td <- tempdir()
  file.copy(ref, file.path(td, "ref.vcf"), overwrite = TRUE)
  file.copy(qry, file.path(td, "query.vcf"), overwrite = TRUE)
  metrics <- benchmarkVariants(
    vcfDir = NULL,
    ref_vcf = "ref.vcf", query_vcf = "query.vcf",
    workDir = td
  )
  # 1 TP (chr1:100), 1 FP (chr1:300), 1 FN (chr1:200)
  expect_equal(as.integer(metrics["ALL", "TP"]), 1)
  expect_equal(as.integer(metrics["ALL", "FP"]), 1)
  expect_equal(as.integer(metrics["ALL", "FN"]), 1)
})

test_that("benchmarkVariants writes feature table", {
  ref <- system.file("extdata", "test_truth.vcf", package = "huGSVcalleR")
  qry <- system.file("extdata", "test_query.vcf", package = "huGSVcalleR")
  skip_if(ref == "" || qry == "", "test VCF not available")
  td <- tempdir()
  file.copy(ref, file.path(td, "ref.vcf"), overwrite = TRUE)
  file.copy(qry, file.path(td, "query.vcf"), overwrite = TRUE)
  benchmarkVariants(
    vcfDir = NULL,
    ref_vcf = "ref.vcf", query_vcf = "query.vcf",
    featureTable = TRUE, workDir = td
  )
  outfile <- file.path(td, "query vs ref, benchmarking results.csv")
  expect_true(file.exists(outfile))
  features <- read.csv(outfile)
  expect_true(all(c("tag", "CHROM", "POS", "REF", "ALT") %in%
                    colnames(features)))
  expect_equal(sort(unique(features$tag)), c("FN", "FP", "TP"))
})
