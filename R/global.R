# global.R
# Global settings and real function wrappers for huGSVcalleR

library(shiny)
library(shinythemes)
library(DT)
library(plotly)
library(ggplot2)
library(dplyr)
library(tidyr)
library(huGSVcalleR)

# Test data paths (adjust these to your actual test data locations)
TEST_PATHS <- list(
  workDir = "D:/Vasily Grinev",
  fastqDir = "Files_FASTQ",
  bamDir = "Files_BAM",
  vcfDir = "Files_VCF",
  referenceDir = "Reference_Genomes",
  fastaDir = "Files_FASTA",
  gtfDir = "Files_GTF"
)

# Real function wrapper for Module 1
run_module1_qc <- function(fastq1, fastq2 = NULL, workDir, ...) {
  setwd(workDir)

  results <- list()

  # QC for read 1
  qa1 <- huGSVcalleR::assessQRawReads(
    fastqDir = TEST_PATHS$fastqDir,
    fastq = fastq1,
    n = NULL,
    adapters = "adapters.txt",
    contaminants = "contaminants.txt",
    workDir = workDir
  )

  # Save results
  rds_file1 <- file.path(workDir, "test_seq_R1_QA.rds")
  saveRDS(qa1, rds_file1)

  # Generate report
  huGSVcalleR::reportQAResults(
    x = qa1,
    output = "test_seq_R1_QA_report",
    workDir = workDir
  )

  # QC for read 2 if provided
  if (!is.null(fastq2)) {
    qa2 <- huGSVcalleR::assessQRawReads(
      fastqDir = TEST_PATHS$fastqDir,
      fastq = fastq2,
      n = NULL,
      adapters = "adapters.txt",
      contaminants = "contaminants.txt",
      workDir = workDir
    )

    rds_file2 <- file.path(workDir, "test_seq_R2_QA.rds")
    saveRDS(qa2, rds_file2)

    huGSVcalleR::reportQAResults(
      x = qa2,
      output = "test_seq_R2_QA_report",
      workDir = workDir
    )

    results$qa2 <- qa2
    results$rds_file2 <- rds_file2
  }

  # Clean reads
  huGSVcalleR::cleanRawReads(
    fastqDir = TEST_PATHS$fastqDir,
    fastq1 = fastq1,
    fastq2 = fastq2,
    adapters = "adapters.txt",
    error = 0.2,
    min_match_flank = 3L,
    anchored = TRUE,
    indels = FALSE,
    tr_score = 20,
    tr_start = 10,
    tr_end = 30,
    phred = "phred_scores.txt",
    k = 3,
    halfwidth = NULL,
    successive = TRUE,
    readl = 100,
    readq = 20,
    readn = 3,
    dustScore = 875,
    batchSize = NA,
    postfix = "filtered",
    workDir = workDir
  )

  results$qa1 <- qa1
  results$rds_file1 <- rds_file1
  results$status <- "Completed"

  return(results)
}

# Real function wrapper for Module 2
run_module2_alignment <- function(filtered_fastq1, filtered_fastq2, workDir, ...) {
  setwd(workDir)

  results <- list()

  # Build index if needed
  index_path <- file.path(TEST_PATHS$referenceDir, "GRCh38")
  if (!file.exists(paste0(index_path, ".00.b.array"))) {
    huGSVcalleR::buildIndexSubread(
      ref_genome = dirname(index_path),
      ref_fasta = TEST_PATHS$fastaDir,
      index = basename(index_path),
      fa = "hg38.fa.gz",
      memory = 8000,
      workDir = workDir
    )
  }

  # Run alignment
  bam_file <- "test_seq.filtered.bam"
  huGSVcalleR::alignDNASubread(
    genome = index_path,
    fastqDir = TEST_PATHS$fastqDir,
    fastq1 = filtered_fastq1,
    fastq2 = filtered_fastq2,
    bamDir = TEST_PATHS$bamDir,
    bamFile = bam_file,
    orientation = "fr",
    threads = 4,
    SV = FALSE,
    workDir = workDir
  )

  # Sort and index
  huGSVcalleR::sortBamFile(
    bamDir = TEST_PATHS$bamDir,
    bamFile = bam_file,
    byQname = FALSE,
    workDir = workDir
  )

  # Assess alignment quality
  qal <- huGSVcalleR::assessQAlignedReads(
    bamDir = TEST_PATHS$bamDir,
    bamFile = bam_file,
    workDir = workDir
  )

  # Save and report
  rds_file <- file.path(workDir, "alignment_quality.rds")
  saveRDS(qal, rds_file)

  huGSVcalleR::reportQAlResults(
    x = qal,
    output = "alignment_report",
    workDir = workDir
  )

  results$alignment_stats <- qal
  results$bam_file <- file.path(TEST_PATHS$bamDir, bam_file)
  results$rds_file <- rds_file
  results$status <- "Completed"

  return(results)
}

# Real function wrapper for Module 4
run_module4_variant_calling <- function(bam_file, workDir, ...) {
  setwd(workDir)

  results <- list()

  # Call variants using Fisher's exact test
  vcf_result <- huGSVcalleR::fisherCalleR(
    bamDir = TEST_PATHS$bamDir,
    bamFile = bam_file,
    fastaDir = TEST_PATHS$fastaDir,
    faFile = "hg38.fa",
    baseq = 20,
    mindepth = 1,
    maxdepth = 1e6,
    qvalue = 12,
    trim = 0,
    threads = 3,
    workDir = workDir
  )

  vcf_file <- file.path(TEST_PATHS$vcfDir,
                        paste0(tools::file_path_sans_ext(bam_file),
                               "_fisherCalleR_SNVs.vcf"))

  results$vcf_file <- vcf_file
  results$status <- "Completed"

  return(results)
}

# Helper function to read and display VCF files
read_vcf_preview <- function(vcf_path, n_lines = 10) {
  if (!file.exists(vcf_path)) {
    return("File not found")
  }

  # Read header lines (start with ##)
  con <- file(vcf_path, "r")
  header_lines <- character()
  data_start <- 0

  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    if (grepl("^##", line)) {
      header_lines <- c(header_lines, line)
      data_start <- data_start + 1
    } else if (grepl("^#CHROM", line)) {
      header_lines <- c(header_lines, line)
      data_start <- data_start + 1
      break
    } else {
      break
    }
  }
  close(con)

  # Read data lines
  if (data_start > 0) {
    data_lines <- readLines(vcf_path)[(data_start + 1):min(data_start + n_lines,
                                                           length(readLines(vcf_path)))]
  } else {
    data_lines <- readLines(vcf_path)[1:min(n_lines, length(readLines(vcf_path)))]
  }

  return(list(header = header_lines, data = data_lines))
}
