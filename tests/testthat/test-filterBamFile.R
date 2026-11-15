test_that("filterBamFile filters BAM by mapq", {
  bam <- system.file("extdata", "test.bam", package = "huGSVcalleR")
  bai <- system.file("extdata", "test.bam.bai", package = "huGSVcalleR")
  skip_if(bam == "" || bai == "", "test BAM not available")
  td <- tempdir()
  file.copy(bam, file.path(td, "test.bam"), overwrite = TRUE)
  file.copy(bai, file.path(td, "test.bam.bai"), overwrite = TRUE)
  filterBamFile(
    bamDir    = NULL,
    bamFile   = "test.bam",
    mapq      = 30,
    overwrite = TRUE,
    workDir   = td
  )
  expect_true(file.exists(file.path(td, "test.filtered.bam")))
})

test_that("filterBamFile filters BAM by region", {
  bam <- system.file("extdata", "test.bam", package = "huGSVcalleR")
  bai <- system.file("extdata", "test.bam.bai", package = "huGSVcalleR")
  reg <- system.file("extdata", "test_regions.txt", package = "huGSVcalleR")
  skip_if(bam == "" || bai == "" || reg == "", "test data not available")
  td <- tempdir()
  file.copy(bam, file.path(td, "test.bam"), overwrite = TRUE)
  file.copy(bai, file.path(td, "test.bam.bai"), overwrite = TRUE)
  file.copy(reg, file.path(td, "test_regions.txt"), overwrite = TRUE)
  filterBamFile(
    bamDir    = NULL,
    bamFile   = "test.bam",
    gr        = "test_regions.txt",
    overwrite = TRUE,
    workDir   = td
  )
  expect_true(file.exists(file.path(td, "test.filtered.bam")))
})

test_that("filterBamFile stops if output exists and overwrite = FALSE", {
  bam <- system.file("extdata", "test.bam", package = "huGSVcalleR")
  bai <- system.file("extdata", "test.bam.bai", package = "huGSVcalleR")
  skip_if(bam == "" || bai == "", "test BAM not available")
  td <- tempdir()
  file.copy(bam, file.path(td, "test.bam"), overwrite = TRUE)
  file.copy(bai, file.path(td, "test.bam.bai"), overwrite = TRUE)
  writeLines("placeholder", file.path(td, "test.filtered.bam"))
  expect_error(
    filterBamFile(
      bamDir    = NULL,
      bamFile   = "test.bam",
      overwrite = FALSE,
      workDir   = td
    ),
    "already exists"
  )
  file.remove(file.path(td, "test.filtered.bam"))
})
