#' Calculate quality scores for short reads
#' @description Calculation of the quality scores for short reads.
#' @param fastqDir character string specifying the name of the directory
#'     containing the FASTQ file(-s). NULL by default that mean the current
#'     working directory.
#' @param fastq character string specifying the name of the FASTQ file
#'     containing short reads. Allowed formats are "fastq" or "fastq.gz".
#' @param n integer value giving the number of reads loaded into the computer's
#'     RAM. The default value is NULL, at which all reads are loaded into RAM.
#' @param adapters character string specifying the name of the tab-delimited
#'     TXT file containing adapter sequences. The default value is NULL, at
#'     which the search for adaptor sequences in short reads is not carried out.
#'     If so, this file must contain the following three fields:
#'     i) adapter_id         - adapter ID;
#'     ii) adapter_name      - adapter name;
#'     iii) adapter_sequence - adapter sequence.
#' @param contaminants character string specifying the name of the
#'     tab-delimited TXT file containing foreign sequences. The default value
#'     is NULL, at which the search for foreign sequences in short reads is not
#'     carried out. If so, this file must contain the following three fields:
#'     i) foreign_id         - foreign sequence ID;
#'     ii) foreign_name      - foreign sequence name;
#'     iii) foreign_sequence - foreign sequence itself.
#' @param workDir character string specifying the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return object of class list containing all main quality scores for short
#'     reads.
#' @author Vasily V. Grinev
#' @examples
#' qa <- assessQRawReads(
#'     fastqDir = "Files_FASTQ",
#'     fastq = "2-Galkina-A_S2_L001_R1_001.fastq.gz",
#'     n = 1e5,
#'     adapters = NULL,
#'     contaminants = NULL,
#'     workDir = "D:/Vasily Grinev"
#' )
#' @export
#' @importFrom Biostrings DNAStringSet PhredQuality width letterFrequency
#' @importFrom ShortRead countFastq qa dustyScore
#' @importFrom data.table data.table rbindlist
#' @importFrom IRanges IntegerList
#' @importFrom stats density quantile
#' @importFrom methods as is
#' @importFrom DelayedArray mean
assessQRawReads <- function(fastqDir = NULL,
                            fastq,
                            n = NULL,
                            adapters = NULL,
                            contaminants = NULL,
                            workDir = NULL) {
  ### Setting the working directory.
  if (is.null(x = workDir)) {
    workDir <- getwd()
  }
  ### Full path to the FASTQ file(-s).
  if (is.null(x = fastqDir)) {
    path <- paste(workDir, fastqDir, sep = "")
  } else {
    path <- paste(workDir, fastqDir, sep = "/")
  }
  ### Retrieving and validation of the FASTQ file(-s) extension.
  frt <- tools::file_ext(x = fastq)
  if (!frt %in% c("fastq", "gz")) {
    stop("Invalid file format")
  }
  ### Collection of the quality metrics.
  ##  Sample name and file name.
  sample_name <- gsub(pattern = ".fastq.+", replacement = "", x = fastq)
  file_name <- fastq
  ##  Total number of reads and bases.
  TNBs <- countFastq(dirPath = path, pattern = fastq)
  TNRs <- TNBs[[1]]
  TNBs <- TNBs[[2]]
  ##  Distribution the length of reads.
  reads <- readFASTQ(
    fastqDir = fastqDir, fastq = fastq, n = n,
    workDir = workDir
  )
  RLs <- .computeReadLengths(reads)
  ##  Per entire reads quality.
  PERQs <- .computePerReadQuality(reads)
  ##  Per cycle quality of reads.
  if (is.null(x = n)) {
    QA <- qa(dirPath = path, pattern = fastq)
  } else {
    QA <- qa(dirPath = path, pattern = fastq, n = n)
  }
  PCQRs <- .computePerCycleQuality(QA)
  ##  Bases composition of the entire set of reads.
  BCs <- .computeBaseComposition(reads)
  ##  Per cycle content of bases.
  PCCBs <- .computePerCycleBases(QA)
  ##  GC composition of the entire set of reads.
  GCs <- .computeGcComposition(reads)
  ##  Summarization the low-complexity of reads.
  DSs <- .computeLowComplexity(reads)
  ##  Frequency of reads.
  FRs <- .computeReadFrequency(QA)
  ##  Overrepresented reads.
  ORs <- .computeOverrepresented(QA)
  ##  Identification of irrelevant sequences.
  ASs <- NA
  if (!is.null(x = adapters)) {
    ASs <- findAdapters(
      fastqDir = fastqDir, fastq = fastq, n = n,
      adapters = adapters, workDir = workDir
    )
  }
  FSs <- NA
  if (!is.null(x = contaminants)) {
    FSs <- findForeignSeqs(
      fastqDir = fastqDir, fastq = fastq, n = n,
      contaminants = contaminants, workDir = workDir
    )
  }
  ### Quality metrics.
  return(list(
    SN = sample_name,
    FN = file_name,
    TNRs = TNRs,
    TNBs = TNBs,
    RLs = RLs,
    PERQs = PERQs,
    PCQRs = PCQRs,
    BCs = BCs,
    PCCBs = PCCBs,
    GCs = GCs,
    DSs = DSs,
    FRs = FRs,
    ORs = ORs,
    ASs = ASs,
    FSs = FSs
  ))
}

#' Compute read length distribution
#' @description Internal helper for assessQRawReads(). Computes the
#'     summary statistics for the distribution of read lengths.
#' @param reads list with DNAStringSet object of read sequences.
#' @return numeric vector with summary statistics.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.computeReadLengths <- function(reads) {
  as.numeric(x = summary(object = width(x = reads$reads)))
}

#' Compute per-read quality metrics
#' @description Internal helper for assessQRawReads(). Computes the
#'     per-read average quality, its density, and quantiles.
#' @param reads list with DNAStringSet and PhredQuality objects.
#' @return list with means, density and metrics.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.computePerReadQuality <- function(reads) {
  PERQs <- mean(x = as(object = reads$QScores, Class = "IntegerList"))
  names(x = PERQs) <- NULL
  list(
    means = PERQs,
    density = do.call(
      what = cbind,
      args = density(x = PERQs)[seq_len(2)]
    ),
    metrics = quantile(
      x = PERQs,
      probs = c(0.1, 0.25, 0.5, 0.75, 0.9)
    )
  )
}

#' Compute per-cycle quality metrics
#' @description Internal helper for assessQRawReads(). Computes the
#'     per-cycle quality metrics from the ShortRead::qa() output.
#' @param QA ShortRead::qa() output.
#' @return data.frame with per-cycle quality metrics.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.computePerCycleQuality <- function(QA) {
  PCQRs <- do.call(
    what = rbind,
    args = lapply(
      X = split(
        x = QA[["perCycle"]]$quality,
        f = QA[["perCycle"]]$quality$Cycle
      ),
      FUN = function(y) {
        rep(x = y[, 3], times = y[, 4])
      }
    )
  )
  m <- list()
  for (i in seq(from = 1, to = nrow(x = PCQRs), by = 5)) {
    if (i + 4 < nrow(x = PCQRs)) {
      m[[i]] <- c(
        paste(i, i:i + 4, sep = "-"),
        quantile(
          x = PCQRs[i:i + 4, ],
          probs = c(0.1, 0.25, 0.5, 0.75, 0.9)
        )
      )
    } else {
      m[[i]] <- c(
        paste(i, nrow(x = PCQRs), sep = "-"),
        quantile(
          x = PCQRs[i:nrow(x = PCQRs), ],
          probs = c(0.1, 0.25, 0.5, 0.75, 0.9)
        )
      )
    }
  }
  PCQRs <- data.frame(do.call(args = m, what = rbind))
  colnames(x = PCQRs) <- c(
    "cycles", "minimum", "Q1", "median", "Q3", "maximum"
  )
  PCQRs[, -1] <- apply(X = PCQRs[, -1], MARGIN = 2, FUN = as.numeric)
  PCQRs
}

#' Compute base composition
#' @description Internal helper for assessQRawReads(). Computes the
#'     nucleotide frequency for each read and the summary statistics.
#' @param reads list with DNAStringSet object of read sequences.
#' @return list with frequency and metrics.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.computeBaseComposition <- function(reads) {
  BCs <- letterFrequency(
    x = reads$reads,
    letters = c("A", "C", "G", "T", "N")
  )
  BCs <- BCs / rowSums(x = BCs)
  metrics <- rbind(
    quantile(x = BCs[, 1], probs = c(0.1, 0.25, 0.5, 0.75, 0.9)),
    quantile(x = BCs[, 2], probs = c(0.1, 0.25, 0.5, 0.75, 0.9)),
    quantile(x = BCs[, 3], probs = c(0.1, 0.25, 0.5, 0.75, 0.9)),
    quantile(x = BCs[, 4], probs = c(0.1, 0.25, 0.5, 0.75, 0.9)),
    quantile(x = BCs[, 5], probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
  )
  metrics <- data.frame(cbind(colnames(x = BCs), metrics))
  colnames(x = metrics) <- c(
    "nucleotide", "minimum", "Q1", "median", "Q3", "maximum"
  )
  metrics[, -1] <- apply(X = metrics[, -1], MARGIN = 2, FUN = as.numeric)
  list(frequency = BCs, metrics = metrics)
}

#' Compute per-cycle base composition
#' @description Internal helper for assessQRawReads(). Computes the
#'     per-cycle nucleotide frequencies from the ShortRead::qa() output.
#' @param QA ShortRead::qa() output.
#' @return data.frame with per-cycle nucleotide frequencies.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.computePerCycleBases <- function(QA) {
  PCCBs <- QA[["perCycle"]]$baseCall[, -4]
  PCCBs$Base <- as.character(x = PCCBs$Base)
  PCCBs <- do.call(
    what = rbind,
    args = lapply(
      X = split(x = PCCBs, f = PCCBs$Cycle),
      FUN = function(y) {
        rep(x = y[, 2], times = y[, 3])
      }
    )
  )
  m <- data.frame(cycles = "", A = 0, C = 0, G = 0, T = 0, N = 0)
  for (i in seq(from = 1, to = nrow(x = PCCBs), by = 5)) {
    if (i + 4 < nrow(x = PCCBs)) {
      bases <- c(
        paste(i, i:i + 4, sep = "-"),
        table(PCCBs[i:i + 4, ])
      )
      m[i, 1] <- bases[1]
      m[i, -1] <- bases[match(colnames(x = m), names(x = bases))][-1]
    } else {
      bases <- c(
        paste(i, nrow(x = PCCBs), sep = "-"),
        table(PCCBs[i:nrow(x = PCCBs), ])
      )
      m[i, 1] <- bases[1]
      m[i, -1] <- bases[match(colnames(x = m), names(x = bases))][-1]
    }
  }
  PCCBs <- m[!is.na(x = m$cycles), ]
  PCCBs[is.na(x = PCCBs$N), ]$N <- 0
  PCCBs[, -1] <- apply(X = PCCBs[, -1], MARGIN = 2, FUN = as.numeric)
  rownames(x = PCCBs) <- NULL
  PCCBs[, -1] <- PCCBs[, -1] / rowSums(x = PCCBs[, -1])
  PCCBs
}

#' Compute GC composition
#' @description Internal helper for assessQRawReads(). Computes the GC
#'     content for each read and the summary statistics.
#' @param reads list with DNAStringSet object of read sequences.
#' @return list with frequency, density and metrics.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.computeGcComposition <- function(reads) {
  GCs <- letterFrequency(x = reads$reads, letters = c("G", "C"))
  GCs <- rowSums(x = GCs) / width(x = reads$reads)
  list(
    frequency = GCs,
    density = do.call(
      what = cbind,
      args = density(x = GCs)[seq_len(2)]
    ),
    metrics = quantile(x = GCs, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
  )
}

#' Compute low-complexity of reads
#' @description Internal helper for assessQRawReads(). Computes the dusty
#'     scores and their summary statistics.
#' @param reads list with DNAStringSet object of read sequences.
#' @return list with scores, density and metrics.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.computeLowComplexity <- function(reads) {
  DSs <- dustyScore(x = reads$reads)
  list(
    scores = DSs,
    density = do.call(
      what = cbind,
      args = density(x = log10(x = DSs))[seq_len(2)]
    ),
    metrics = quantile(x = DSs, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
  )
}

#' Compute read frequency
#' @description Internal helper for assessQRawReads(). Computes the
#'     frequency of each read occurrence.
#' @param QA ShortRead::qa() output.
#' @return data.frame with occurrence, reads and frequency.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.computeReadFrequency <- function(QA) {
  FRs <- QA[["sequenceDistribution"]][, seq_len(2)]
  colnames(x = FRs) <- c("occurrence", "reads")
  FRs$frequency <- FRs$reads / sum(x = FRs$reads)
  FRs
}

#' Compute overrepresented reads
#' @description Internal helper for assessQRawReads(). Extracts the
#'     overrepresented reads from the ShortRead::qa() output.
#' @param QA ShortRead::qa() output.
#' @return data.frame with sequence and count.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.computeOverrepresented <- function(QA) {
  ShortRead:::.freqSequences(QA, "read")[, seq_len(2)]
}


#' Read records from FASTQ-formatted file
#' @description Reading the records from FASTQ-formatted file.
#' @param fastqDir character string specifying the name of the directory
#'     containing the FASTQ file(s). NULL by default that mean the current
#'     working directory.
#' @param fastq character string specifying the name of the FASTQ file
#'     containing short reads. Allowed formats are "fastq" or "fastq.gz".
#' @param n integer value giving the number of reads loaded into the computer's
#'     RAM. The default value is NULL, at which all reads are loaded into RAM.
#' @param workDir character string specifying the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return object of class list containing DNAStringSet object of read
#'     sequences and PhredQuality object of read quality scores.
#' @author Vasily V. Grinev
#' @examples
#' fq <- readFASTQ(
#'     fastqDir = "Files_FASTQ",
#'     fastq = "example_seq.read1.fastq.gz",
#'     n = 1e5,
#'     workDir = "D:/Vasily Grinev"
#' )
#' @export
#' @importFrom Biostrings DNAStringSet PhredQuality
#' @importFrom tools file_ext
readFASTQ <- function(
    fastqDir = NULL,
    fastq,
    n = NULL,
    workDir = NULL) {
  ### Setting the working directory.
  if (is.null(x = workDir)) {
    workDir <- getwd()
  }
  ### Full path to the FASTQ file(-s).
  if (is.null(x = fastqDir)) {
    path <- paste(workDir, fastqDir, sep = "")
  } else {
    path <- paste(workDir, fastqDir, sep = "/")
  }
  ### Retrieving and validation of the FASTQ file(-s) extension.
  frt <- tools::file_ext(x = fastq)
  if (!frt %in% c("fastq", "gz")) {
    stop("Invalid file format")
  }
  ### Loading the FASTQ data.
  if (is.null(x = n)) {
    reads <- readLines(con = paste(path, fastq, sep = "/"))
    FQs <- .parseFastqLines(reads, n = NULL)
  } else {
    reads <- readLines(con = paste(path, fastq, sep = "/"), n = n * 4)
    FQs <- .parseFastqLines(reads, n = n)
  }
  rm(reads)
  ### Returning the final object.
  return(FQs)
}

#' Parse FASTQ lines into sequences and quality scores
#' @description Internal helper for readFASTQ(). Extracts the sequences and
#'     quality scores from the raw FASTQ lines.
#' @param reads character vector with raw FASTQ lines.
#' @param n integer value giving the number of reads, or NULL.
#' @return list with DNAStringSet and PhredQuality objects.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.parseFastqLines <- function(reads, n) {
  if (is.null(x = n)) {
    to <- length(x = reads)
  } else {
    to <- n * 4
  }
  FQs <- list(
    reads = DNAStringSet(x = reads[seq(from = 2, to = to, by = 4)]),
    QScores = PhredQuality(x = reads[seq(from = 4, to = to, by = 4)])
  )
  names(x = FQs$reads) <- reads[seq(from = 1, to = to, by = 4)]
  names(x = FQs$QScores) <- reads[seq(from = 1, to = to, by = 4)]
  FQs
}


#' Clean up low quality short reads
#' @description Generic function for trimming and filtering of low-quality
#'     DNA- and/or RNA-Seq short reads.
#' @param fastqDir character string giving the name (or path to and name) of
#'     directory with FASTQ file(-s). NULL by default, which means the current
#'     working directory will be used instead of the specified FASTQ directory.
#' @param fastq1 character string giving the name of FASTQ file containing
#'     short reads. For paired-end reads, this gives the name of FASTQ file for
#'     first reads in each sequencing library.
#' @param fastq2 character string including the name of FASTQ file containing
#'     second short reads. NULL by default, which means that the data will be
#'     processed as single-read ones.
#' @param adapters character string specifying the name of the tab-delimited
#'     TXT file containing adapter sequences. The default value is NULL, at
#'     which the trimming of adaptor sequence(-s) in short reads is not carried
#'     out. If so, this file must contain the following three fields:
#'     i) adapter_id         - adapter ID;
#'     ii) adapter_name      - adapter name;
#'     iii) adapter_sequence - adapter sequence.
#' @param error (when argument adapters is not NULL) numeric value in [0, 1]
#'     specifying error rate. The error rate is the proportion of mismatches
#'     allowed between the adapter and the aligned portion of the read. For a
#'     given adapter A, the number of allowed mismatches between each
#'     subsequence s of A and the read is computed as: error * L_s, where L_s
#'     is the length of the subsequence s. Default value is 0.2.
#' @param min_match_flank (when argument adapters is not NULL) integer value
#'     giving the number of nucleotides. Do not trim in flanks of the read if
#'     a match has min_match_flank of less length. Default value is 3L (trim
#'     only with >=4 coincidences in a flanking match).
#' @param anchored (when argument adapters is not NULL) logical, TRUE by
#'     default. Can the adapter or partial adapter be within (anchored=FALSE)
#'     or only in the terminal region(-s) (anchored=TRUE) of the read?
#' @param indels (when argument adapters is not NULL) logical, FALSE by
#'     default. If TRUE, indels are allowed in the alignments of the suffixes
#'     of adapter with the read, at its beginning.
#' @param tr_score integer value giving the quality score at or below of which
#'     a nucleotide will be marked as removable. NULL by default, which means
#'     that fixed (but nor dynamic) trimming will be applied to the low quality
#'     ends of short reads.
#' @param tr_start integer specifying the number of nucleotides to be removed
#'     from the left end of short reads. This argument is used only when
#'     tr_score=NULL and fixed trimming is carrying out. Default value is 10.
#' @param tr_end integer specifying the number of nucleotides to be removed
#'     from the right end of short reads. This argument is used only when
#'     tr_score=NULL and fixed trimming is carrying out. Default value is 30.
#' @param phred temporary argument for internal use.
#' @param k integer value describing the number of failing nucleotides required
#'     to trigger trimming. This argument is used only with dynamic trimming.
#'     Default value is 3.
#' @param halfwidth half width (cycles before or after the current nucleotide
#'     position) in which qualities are assessed. This argument is used only
#'     with dynamic trimming in the sliding window approach. Default value is
#'     NULL, what determines the choice of the successive approach in dynamic
#'     trimming of short reads.
#' @param successive logical value indicating whether failures can occur
#'     anywhere in the sequence, or must be successive. If successive=FALSE,
#'     then the k’th failed nucleotide and subsequent are removed. If
#'     successive=TRUE (by default), the first succession of k failed and
#'     subsequent nucleotides are removed. This argument is used only with
#'     dynamic trimming in the successive approach.
#' @param readl integer specifying the minimum length of reads to be saved.
#'     Zero is the default, which will allow reads of any length to be saved.
#' @param readq integer specifying the minimum quality of reads to be saved.
#'     Zero is the default, which will allow reads of any quality to be saved.
#' @param readn integer specifying the maximum number of ambiguous nucleotides
#'     allowed in a read to be saved. Zero is the default.
#' @param dustScore integer specifying the maximum dust score at which reads to
#'     be saved. By default this argument is NULL and such filter is not used.
#' @param batchSize NA (by default) or integer value indicating the maximum
#'     number of reads to be processed at any one time during filtration by
#'     dustScore filter. By default, all reads are processed simultaneously
#'     Smaller values use less memory but are computationally less efficient.
#' @param postfix character string that will be added to the output file name
#'     as postfix. "filtered" by default.
#' @param overwrite logical value indicating whether existing output files
#'     with the same names should be overwritten. FALSE by default, in which
#'     case the function stops with an error if an output file already exists.
#' @param workDir character string specifying the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return filtered FASTQ file(-s) will be written to the same directory as
#'     input FASTQ file(-s).
#' @author Vasily V. Grinev.
#' @examples
#' cleanRawReads(
#'     fastqDir = "Files_FASTQ",
#'     fastq1 = "example_seq.read1.fastq",
#'     fastq2 = "example_seq.read2.fastq",
#'     adapters = "adapters.txt",
#'     error = 0.2,
#'     min_match_flank = 3L,
#'     anchored = TRUE,
#'     indels = FALSE,
#'     tr_score = 20,
#'     tr_start = 10,
#'     tr_end = 30,
#'     phred = "phred_scores.txt",
#'     k = 3,
#'     halfwidth = 5,
#'     successive = TRUE,
#'     readl = 100,
#'     readq = 20,
#'     readn = 3,
#'     dustScore = 875,
#'     batchSize = NA,
#'     postfix = "filtered",
#'     overwrite = TRUE,
#'     workDir = "D:/Vasily Grinev"
#' )
#' @export
#' @importFrom ShortRead readFastq writeFastq trimEnds trimTailw trimTails id
#' @importFrom ShortRead polynFilter occurrenceFilter dustyFilter
#' @importFrom Biostrings width quality
#' @importFrom methods as
#' @importFrom utils read.table
#' @importFrom IRanges narrow
cleanRawReads <- function(
    fastqDir = NULL,
    fastq1,
    fastq2 = NULL,
    adapters = NULL,
    error = 0.2,
    min_match_flank = 3L,
    anchored = TRUE,
    indels = FALSE,
    tr_score = NULL,
    tr_start = 10,
    tr_end = 30,
    phred = NA,
    k = 3,
    halfwidth = NULL,
    successive = TRUE,
    readl = 0,
    readq = 0,
    readn = 0,
    dustScore = NULL,
    batchSize = NA,
    postfix = "filtered",
    overwrite = FALSE,
    workDir = NULL) {
  ### Setting the working directory.
  if (is.null(x = workDir)) {
    workDir <- getwd()
  }
  ### Full path to the FASTQ file(-s).
  if (is.null(x = fastqDir)) {
    path <- workDir
  } else {
    path <- paste(workDir, fastqDir, sep = "/")
  }
  ### Retrieving and validation of the FASTQ file(-s) extension.
  frt1 <- tools::file_ext(x = fastq1)
  if (!frt1 %in% c("fastq", "gz")) {
    stop("Invalid format of FASTQ file 1")
  }
  if (!is.null(x = fastq2)) {
    frt2 <- tools::file_ext(x = fastq2)
    if (!frt2 %in% c("fastq", "gz")) {
      stop("Invalid format of FASTQ file 2")
    }
  }
  ### Loading of the FASTQ file 1.
  FQ1 <- readFastq(dirPath = path, pattern = fastq1, withIds = TRUE)
  ### Raw reads processing.
  if (!is.null(x = adapters)) {
    FQ1 <- cutLeftAdapter(
      fastq = FQ1, adapters = adapters, error = error,
      min_match_flank = min_match_flank, anchored = anchored,
      indels = indels, workDir = workDir
    )
    FQ1 <- cutRigthAdapter(
      fastq = FQ1, adapters = adapters, error = error,
      min_match_flank = min_match_flank, anchored = anchored,
      indels = indels, workDir = workDir
    )
  }
  if (is.null(x = tr_score)) {
    FQ1 <- .trimFixedLowQuality(FQ1, tr_start, tr_end)
  } else {
    PhSs <- read.table(
      file = paste(workDir, phred, sep = "/"),
      sep = "\t", header = TRUE, quote = "",
      as.is = TRUE, comment.char = "z"
    )
    FQ1 <- .trimDynamicLowQuality(
      FQ1, PhSs, tr_score, k, halfwidth, successive
    )
  }
  FQ1 <- .applyReadFilters(FQ1, readl, readq, readn, dustScore, batchSize)
  if (!is.null(x = fastq2)) {
    FQ2 <- readFastq(dirPath = path, pattern = fastq2, withIds = TRUE)
    if (!is.null(x = adapters)) {
      FQ2 <- cutLeftAdapter(
        fastq = FQ2, adapters = adapters, error = error,
        min_match_flank = min_match_flank, anchored = anchored,
        indels = indels, workDir = workDir
      )
      FQ2 <- cutRigthAdapter(
        fastq = FQ2, adapters = adapters, error = error,
        min_match_flank = min_match_flank, anchored = anchored,
        indels = indels, workDir = workDir
      )
    }
    if (is.null(x = tr_score)) {
      FQ2 <- .trimFixedLowQuality(FQ2, tr_start, tr_end)
    } else {
      FQ2 <- .trimDynamicLowQuality(
        FQ2, PhSs, tr_score, k, halfwidth, successive
      )
    }
    FQ2 <- .applyReadFilters(
      FQ2, readl, readq, readn, dustScore, batchSize
    )
    synced <- .synchronizePairedReads(FQ1, FQ2)
    FQ1 <- synced$FQ1
    FQ2 <- synced$FQ2
  }
  ### Writing the filtered FASTQ file(-s).
  output1 <- paste(path,
                   gsub(
                     pattern = "fast.+",
                     replacement = paste(postfix, "fastq.gz", sep = "."),
                     x = fastq1
                   ),
                   sep = "/"
  )
  .checkOverwrite(output1, overwrite)
  writeFastq(object = FQ1, file = output1, compress = TRUE)
  if (!is.null(x = fastq2)) {
    output2 <- paste(path,
                     gsub(
                       pattern = "fast.+",
                       replacement = paste(postfix, "fastq.gz", sep = "."),
                       x = fastq2
                     ),
                     sep = "/"
    )
    .checkOverwrite(output2, overwrite)
    writeFastq(object = FQ2, file = output2)
  }
}

#' Trim low-quality ends with fixed trimming
#' @description Internal helper for cleanRawReads(). Removes a fixed number
#'     of nucleotides from each end of the reads.
#' @param FQ ShortReadQ object with reads.
#' @param tr_start integer, number of nucleotides to remove from the left.
#' @param tr_end integer, number of nucleotides to remove from the right.
#' @return ShortReadQ object with trimmed reads.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.trimFixedLowQuality <- function(FQ, tr_start, tr_end) {
  keep <- width(x = FQ) > (tr_start + tr_end)
  FQ <- FQ[keep]
  start_pos <- tr_start + 1
  end_pos <- if (tr_end > 0) -tr_end else -1L
  narrow(x = FQ, start = start_pos, end = end_pos)
}

#' Trim low-quality ends with dynamic trimming
#' @description Internal helper for cleanRawReads(). Trims the reads based
#'     on a quality threshold using either the sliding-window or the
#'     successive approach.
#' @param FQ ShortReadQ object with reads.
#' @param PhSs data.frame with ASCII quality scores.
#' @param tr_score integer, quality threshold.
#' @param k integer, number of failing nucleotides.
#' @param halfwidth integer, sliding-window half-width, or NULL.
#' @param successive logical, use successive approach?
#' @return ShortReadQ object with trimmed reads.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.trimDynamicLowQuality <- function(FQ, PhSs, tr_score, k, halfwidth,
                                   successive) {
  FQ <- trimEnds(
    object = FQ,
    a = PhSs[PhSs$score == tr_score, ]$ASCII,
    left = TRUE, right = FALSE,
    relation = "<="
  )
  if (!is.null(x = halfwidth)) {
    FQ <- trimTailw(
      object = FQ, k = k,
      a = PhSs[PhSs$score == tr_score, ]$ASCII,
      halfwidth = halfwidth
    )
  } else {
    FQ <- trimTails(
      object = FQ, k = k,
      a = PhSs[PhSs$score == tr_score, ]$ASCII,
      successive = successive
    )
  }
  FQ
}

#' Apply read-level filters
#' @description Internal helper for cleanRawReads(). Applies length,
#'     quality, ambiguous-nucleotide, duplicate and low-complexity filters.
#' @param FQ ShortReadQ object with reads.
#' @param readl integer, minimum read length.
#' @param readq integer, minimum read quality.
#' @param readn integer, maximum number of ambiguous nucleotides.
#' @param dustScore integer, maximum dust score, or NULL.
#' @param batchSize integer, batch size for dusty filter, or NA.
#' @return ShortReadQ object with filtered reads.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.applyReadFilters <- function(FQ, readl, readq, readn, dustScore,
                              batchSize) {
  FQ <- FQ[width(x = FQ) >= readl]
  FQ <- FQ[rowMeans(
    x = as(quality(x = FQ), "matrix"),
    na.rm = TRUE
  ) >= readq]
  FQ <- FQ[polynFilter(
    threshold = as.integer(x = readn),
    nuc = "other"
  )(FQ)]
  FQ <- FQ[occurrenceFilter(
    min = 1L, max = 1L,
    withSread = TRUE,
    duplicates = "sample"
  )(FQ)]
  if (!is.null(x = dustScore)) {
    FQ <- FQ[dustyFilter(
      threshold = dustScore,
      batchSize = batchSize
    )(FQ)]
  }
  FQ
}

#' Synchronize paired-end reads
#' @description Internal helper for cleanRawReads(). Keeps only the read
#'     pairs that are present in both FASTQ files after filtering.
#' @param FQ1 ShortReadQ object with first reads.
#' @param FQ2 ShortReadQ object with second reads.
#' @return list with synchronized FQ1 and FQ2.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.synchronizePairedReads <- function(FQ1, FQ2) {
  ID1 <- gsub(pattern = " 1.+", replacement = "", x = id(object = FQ1))
  ID2 <- gsub(pattern = " 2.+", replacement = "", x = id(object = FQ2))
  FQ1 <- FQ1[ID1 %in% ID2]
  FQ2 <- FQ2[ID2 %in% ID1]
  list(FQ1 = FQ1, FQ2 = FQ2)
}

#' Check overwrite flag for output file
#' @description Internal helper for cleanRawReads(). Removes the output file
#'     if overwrite = TRUE, otherwise stops with an error if the file exists.
#' @param output character string giving the name of the output file.
#' @param overwrite logical, whether to overwrite.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.checkOverwrite <- function(output, overwrite) {
  if (file.exists(output)) {
    if (isTRUE(overwrite)) {
      file.remove(output)
    } else {
      stop(
        "Output file already exists: ",
        output,
        ". Set overwrite = TRUE to replace it."
      )
    }
  }
  invisible(NULL)
}


#' Cut out adapter sequence(-s) from DNA- and/or RNA-Seq short reads
#' @description Generic function for trimming of adapter sequence(-s) from
#'     DNA- and/or RNA-Seq short reads. This function is based on low-level
#'     functions cutLseq() of R/Bioconductor package FastqCleaner and
#'     trimLRpattern() of R/Bioconductor package Biostrings.
#' @param fastqDir character string giving the name (or path to and name) of
#'     directory with FASTQ file(-s). NULL by default, which means the current
#'     working directory will be used instead of the specified FASTQ directory.
#' @param fastq character string giving the name of FASTQ file containing short
#'     reads of interest. Alternatively, fastq can be an internal R object with
#'     preloaded FASTQ data.
#' @param adapters character string specifying the name of the tab-delimited
#'     TXT file containing adapter sequences. This file must contain the
#'     following three fields:
#'     i) adapter_id         - adapter ID;
#'     ii) adapter_name      - adapter name;
#'     iii) adapter_sequence - adapter sequence.
#' @param error numeric value in [0, 1] specifying error rate. The error rate
#'     is the proportion of mismatches allowed between the adapter and the
#'     aligned portion of the read. For a given adapter A, the number of
#'     allowed mismatches between each subsequence s of A and the read is
#'     computed as: error * L_s, where L_s is the length of the subsequence s.
#'     Default value is 0.2.
#' @param min_match_flank integer value giving the number of nucleotides.
#'     Do not trim in flanks of the read if a match has min_match_flank of less
#'     length. Default value is 3L (trim only with >=4 coincidences in a
#'     flanking match).
#' @param anchored logical, TRUE by default. Can the adapter or partial adapter
#'     be within (anchored=FALSE) or only in the terminal region(-s)
#'     (anchored=TRUE) of the read?
#' @param indels logical, FALSE by default. If TRUE, indels are allowed in the
#'     alignments of the suffixes of adapter with the read, at its beginning.
#' @param workDir character string specifying the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return internal R object containing short reads with trimmed adapter
#'     sequence(-s).
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @examples
#' trimL <- cutLeftAdapter(
#'     fastqDir = "Files_FASTQ",
#'     fastq = "example_seq.read1.fastq",
#'     adapters = "adapters.txt",
#'     error = 0.2,
#'     min_match_flank = 3L,
#'     anchored = TRUE,
#'     indels = FALSE,
#'     workDir = "D:/Vasily Grinev"
#' )
#' @export
#' @importFrom ShortRead readFastq
#' @importFrom Biostrings trimLRPatterns width
#' @importFrom utils read.table
#' @importFrom methods is
cutLeftAdapter <- function(
    fastqDir = NULL,
    fastq,
    adapters,
    error = 0.2,
    min_match_flank = 3L,
    anchored = TRUE,
    indels = FALSE,
    workDir = NULL) {
  ### Setting the working directory.
  if (is.null(x = workDir)) {
    workDir <- getwd()
  }
  ### Retrieving of the FASTQ data.
  FQ <- .loadFastqObject(fastqDir, fastq, workDir)
  ### Loading of the adapter sequence(-s).
  ADs <- read.table(
    file = paste(workDir, adapters, sep = "/"),
    sep = "\t", header = TRUE, quote = "\"", as.is = TRUE
  )
  ### Left-end trimming of adapter sequence(-s).
  for (i in seq_len(nrow(x = ADs))) {
    FQ <- .trimLeftAdapterSingle(FQ, ADs[i, 3], error,
                                 min_match_flank, anchored, indels)
  }
  ### Returning the final object.
  return(FQ)
}

#' Load FASTQ object from file or return existing object
#' @description Internal helper for cutLeftAdapter() and cutRigthAdapter().
#'     Reads a FASTQ file if 'fastq' is a character string, or returns it
#'     as-is if it is already a ShortReadQ object.
#' @param fastqDir character string giving the directory name.
#' @param fastq character string with file name or ShortReadQ object.
#' @param workDir character string giving the path to working directory.
#' @return ShortReadQ object.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.loadFastqObject <- function(fastqDir, fastq, workDir) {
  if (is.object(x = fastq)) {
    return(fastq)
  }
  if (is.null(x = fastqDir)) {
    path <- paste(workDir, fastqDir, sep = "")
  } else {
    path <- paste(workDir, fastqDir, sep = "/")
  }
  frt <- tools::file_ext(x = fastq)
  if (!frt %in% c("fastq", "gz")) {
    stop("Invalid format of FASTQ file")
  }
  readFastq(dirPath = path, pattern = fastq, withIds = TRUE)
}

#' Trim a single adapter from the left end
#' @description Internal helper for cutLeftAdapter(). Trims a single
#'     adapter sequence from the 5' end of the reads.
#' @param FQ ShortReadQ object.
#' @param seq_l character string with adapter sequence.
#' @param error numeric error rate.
#' @param min_match_flank integer minimal match length.
#' @param anchored logical, anchored search?
#' @param indels logical, allow indels?
#' @return ShortReadQ object with trimmed reads.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.trimLeftAdapterSingle <- function(FQ, seq_l, error, min_match_flank,
                                   anchored, indels) {
  if (error > 0) {
    flank_seq <- as.integer(
      x = seq_len(length.out = nchar(x = seq_l)) * error
    )
  } else {
    flank_seq <- rep(
      x = 0,
      times = length(x = seq_len(length.out = nchar(x = seq_l)))
    )
  }
  if (min_match_flank >= 1L) {
    if (nchar(x = seq_l) > min_match_flank) {
      flank_seq[seq_len(length.out = min_match_flank)] <- -1
    } else {
      return(FQ)
    }
  }
  if (!anchored) {
    maxlen <- max(width(x = FQ)) - nchar(x = seq_l)
    if (maxlen > 0) {
      seq_l <- paste0(
        paste(rep(x = "N", times = maxlen), collapse = ""),
        seq_l
      )
    }
    flank_seq <- c(flank_seq, rep(x = 0, times = maxlen))
  }
  trimLRPatterns(
    Lpattern = seq_l, subject = FQ,
    max.Lmismatch = flank_seq, with.Lindels = indels
  )
}

#' Cut out adapter sequence(-s) from DNA- and/or RNA-Seq short reads
#' @description Generic function for trimming of adapter sequence(-s) from
#'     DNA- and/or RNA-Seq short reads. This function is based on low-level
#'     functions cutRseq() of R/Bioconductor package FastqCleaner and
#'     trimLRpattern() of R/Bioconductor package Biostrings.
#' @param fastqDir character string giving the name (or path to and name) of
#'     directory with FASTQ file(-s). NULL by default, which means the current
#'     working directory will be used instead of the specified FASTQ directory.
#' @param fastq character string giving the name of FASTQ file containing short
#'     reads of interest. Alternatively, fastq can be an internal R object with
#'     preloaded FASTQ data.
#' @param adapters character string specifying the name of the tab-delimited
#'     TXT file containing adapter sequences. This file must contain the
#'     following three fields:
#'     i) adapter_id         - adapter ID;
#'     ii) adapter_name      - adapter name;
#'     iii) adapter_sequence - adapter sequence.
#' @param error numeric value in [0, 1] specifying error rate. The error rate
#'     is the proportion of mismatches allowed between the adapter and the
#'     aligned portion of the read. For a given adapter A, the number of
#'     allowed mismatches between each subsequence s of A and the read is
#'     computed as: error * L_s, where L_s is the length of the subsequence s.
#'     Default value is 0.2.
#' @param min_match_flank integer value giving the number of nucleotides.
#'     Do not trim in flanks of the read if a match has min_match_flank of less
#'     length. Default value is 3L (trim only with >=4 coincidences in a
#'     flanking match).
#' @param anchored logical, TRUE by default. Can the adapter or partial adapter
#'     be within (anchored=FALSE) or only in the terminal region(-s)
#'     (anchored=TRUE) of the read?
#' @param indels logical, FALSE by default. If TRUE, indels are allowed in the
#'     alignments of the suffixes of adapter with the read, at its beginning.
#' @param workDir character string specifying the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return internal R object containing short reads with trimmed adapter
#'     sequence(-s).
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @examples
#' trimL <- cutRigthAdapter(
#'     fastqDir = "Files_FASTQ",
#'     fastq = "example_seq.read1.fastq",
#'     adapters = "adapters.txt",
#'     error = 0.2,
#'     min_match_flank = 3L,
#'     anchored = TRUE,
#'     indels = FALSE,
#'     workDir = "D:/Vasily Grinev"
#' )
#' @export
#' @importFrom ShortRead readFastq trimLRPatterns
#' @importFrom Biostrings width
#' @importFrom methods is
#' @importFrom utils read.table
cutRigthAdapter <- function(
    fastqDir = NULL,
    fastq,
    adapters,
    error = 0.2,
    min_match_flank = 3L,
    anchored = TRUE,
    indels = FALSE,
    workDir = NULL) {
  ### Setting the working directory.
  if (is.null(x = workDir)) {
    workDir <- getwd()
  }
  ### Retrieving of the FASTQ data.
  FQ <- .loadFastqObject(fastqDir, fastq, workDir)
  ### Loading of the adapter sequence(-s).
  ADs <- read.table(
    file = paste(workDir, adapters, sep = "/"),
    sep = "\t", header = TRUE, quote = "\"", as.is = TRUE
  )
  ### Right-end trimming of adapter sequence(-s).
  for (i in seq_len(nrow(x = ADs))) {
    FQ <- .trimRightAdapterSingle(FQ, ADs[i, 3], error,
                                  min_match_flank, anchored, indels)
  }
  ### Returning the final object.
  return(FQ)
}

#' Trim a single adapter from the right end
#' @description Internal helper for cutRigthAdapter(). Trims a single
#'     adapter sequence from the 3' end of the reads.
#' @param FQ ShortReadQ object.
#' @param seq_l character string with adapter sequence.
#' @param error numeric error rate.
#' @param min_match_flank integer minimal match length.
#' @param anchored logical, anchored search?
#' @param indels logical, allow indels?
#' @return ShortReadQ object with trimmed reads.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.trimRightAdapterSingle <- function(FQ, seq_l, error, min_match_flank,
                                    anchored, indels) {
  if (error > 0) {
    flank_seq <- as.integer(
      x = seq_len(length.out = nchar(x = seq_l)) * error
    )
  } else {
    flank_seq <- rep(
      x = 0,
      times = length(x = seq_len(length.out = nchar(x = seq_l)))
    )
  }
  if (min_match_flank >= 1L) {
    if (nchar(x = seq_l) > min_match_flank) {
      flank_seq[seq_len(length.out = min_match_flank)] <- -1
    } else {
      return(FQ)
    }
  }
  if (!anchored) {
    maxlen <- max(width(x = FQ)) - nchar(x = seq_l)
    if (maxlen > 0) {
      seq_l <- paste0(
        seq_l,
        paste(rep(x = "N", times = maxlen), collapse = "")
      )
    }
    flank_seq <- c(flank_seq, rep(x = 0, times = maxlen))
  }
  trimLRPatterns(
    Rpattern = seq_l, subject = FQ,
    max.Rmismatch = flank_seq, with.Rindels = indels
  )
}


#' Identify adapter sequences in short reads
#' @description Identification of the adapter sequences in short reads.
#' @param fastqDir character string specifying the name of the directory
#'     containing the FASTQ file(s). NULL by default that mean the current
#'     working directory.
#' @param fastq character string specifying the name of the FASTQ file
#'     containing short reads. Allowed formats are "fastq" or "fastq.gz".
#' @param n integer value giving the number of reads loaded into the computer's
#'     RAM. The default value is NULL, at which all reads are loaded into RAM.
#' @param adapters character string specifying the name of the tab-delimited
#'     TXT file containing adapter sequences. This file must contain the
#'     following three fields:
#'     i) adapter_id         - adapter ID;
#'     ii) adapter_name      - adapter name;
#'     iii) adapter_sequence - adapter sequence.
#' @param workDir character string specifying the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return object of class data frame with binned (by nucleotide positions)
#'     frequency data on the presence of adapter sequences in short reads.
#' @author Vasily V. Grinev
#' @examples
#' ad <- findAdapters(
#'     fastqDir = "Files_FASTQ",
#'     fastq = "2-Galkina-A_S2_L001_R1_001.fastq.gz",
#'     n = 1e5,
#'     adapters = "adapters.txt",
#'     workDir = "D:/Vasily Grinev"
#' )
#' @export
#' @importFrom Biostrings vmatchPattern width
#' @importFrom GenomicRanges GRanges intersect
#' @importFrom methods is
#' @importFrom utils read.table
#' @importFrom tools file_ext
findAdapters <- function(fastqDir = NULL,
                         fastq,
                         n = NULL,
                         adapters,
                         workDir = NULL) {
  ### Setting the working directory.
  if (is.null(x = workDir)) {
    workDir <- getwd()
  }
  ### Retrieving and validation of the FASTQ file(-s) extension.
  frt <- tools::file_ext(x = fastq)
  if (!frt %in% c("fastq", "gz")) {
    stop("Invalid file format")
  }
  ### Loading the FASTQ data.
  reads <- readFASTQ(
    fastqDir = fastqDir, fastq = fastq, n = n, workDir = workDir
  )[[1]]
  ### Loading the adapter sequences.
  ADs <- read.table(
    file = paste(workDir, adapters, sep = "/"),
    sep = "\t", header = TRUE, quote = "\"", as.is = TRUE
  )
  ### Computing the hits.
  hits <- lapply(
    seq_len(nrow(x = ADs)),
    function(i) .computeAdapterHits(
      reads, ADs$adapter_sequence[i]
    )
  )
  ### Building the output data.frame.
  positions <- .buildPositionBins(reads)
  hits <- data.frame(cbind(positions, do.call(what = cbind, args = hits)))
  colnames(x = hits) <- c("positions", ADs$adapter_name)
  if (ncol(hits) > 1) {
    hits[, -1] <- apply(
      X = hits[, -1, drop = FALSE],
      MARGIN = 2, FUN = as.numeric
    )
  }
  ### Returning the final object.
  return(hits)
}

#' Compute hits for a single adapter
#' @description Internal helper for findAdapters() and findForeignSeqs().
#'     Computes the binned frequency of matches for a single adapter.
#' @param reads DNAStringSet object with read sequences.
#' @param pattern character string with the adapter sequence.
#' @return numeric vector with binned frequencies.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.computeAdapterHits <- function(reads, pattern) {
  x <- unlist(x = vmatchPattern(
    pattern = pattern, max.mismatch = 2, subject = reads
  ))
  x <- GRanges(seqnames = names(x = x), ranges = x)
  x <- data.frame(intersect(x = range(x), y = x))[, 2:3]
  x <- unlist(x = apply(
    X = x, MARGIN = 1,
    FUN = function(y) {
      seq(from = y[1], to = y[2])
    }
  ))
  x <- cut(
    x = x,
    breaks = seq(from = 1, to = max(width(x = reads)) + 5, by = 5)
  )
  x <- data.frame(table(x))
  x$Freq / (length(x = reads) * 6)
}

#' Build position bins for adapter and contaminant hits
#' @description Internal helper for findAdapters() and findForeignSeqs().
#'     Builds the position labels for the binned frequency data.
#' @param reads DNAStringSet object with read sequences.
#' @return character vector with position labels.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.buildPositionBins <- function(reads) {
  if (max(width(x = reads)) %% 5 == 0) {
    paste(
      seq(from = 1, to = max(width(x = reads)), by = 5),
      seq(from = 5, to = max(width(x = reads)), by = 5),
      sep = "-"
    )
  } else {
    paste(
      seq(from = 1, to = max(width(x = reads)), by = 5),
      c(
        seq(from = 5, to = max(width(x = reads)), by = 5),
        max(width(x = reads))
      ),
      sep = "-"
    )
  }
}


#' Identify foreign sequences in short reads
#' @description Identification of the foreign sequences in short reads.
#' @param fastqDir character string specifying the name of the directory
#'     containing the FASTQ file(s). NULL by default that mean the current
#'     working directory.
#' @param fastq character string specifying the name of the FASTQ file
#'     containing short reads. Allowed formats are "fastq" or "fastq.gz".
#' @param n integer value giving the number of reads loaded into the computer's
#'     RAM. The default value is NULL, at which all reads are loaded into RAM.
#' @param contaminants character string specifying the name of the
#'     tab-delimited TXT file containing foreign sequences. The default value
#'     is NULL, at which the search for foreign sequences in short reads is not
#'     carried out. If so, this file must contain the following three fields:
#'     i) foreign_id         - foreign sequence ID;
#'     ii) foreign_name      - foreign sequence name;
#'     iii) foreign_sequence - foreign sequence itself.
#' @param workDir character string specifying the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return object of class data frame with binned (by nucleotide positions)
#'     frequency data on the presence of foreign sequences in short reads.
#' @author Vasily V. Grinev
#' @examples
#' conts <- findForeignSeqs(
#'     fastqDir = "Files_FASTQ",
#'     fastq = "2-Galkina-A_S2_L001_R1_001.fastq.gz",
#'     n = 1e5,
#'     contaminants = "contaminants.txt",
#'     workDir = "D:/Vasily Grinev"
#' )
#' @export
#' @importFrom Biostrings vmatchPattern width
#' @importFrom GenomicRanges GRanges intersect
#' @importFrom utils read.table
#' @importFrom methods is
#' @importFrom tools file_ext
findForeignSeqs <- function(
    fastqDir = NULL,
    fastq,
    n = NULL,
    contaminants,
    workDir = NULL) {
  ### Setting the working directory.
  if (is.null(x = workDir)) {
    workDir <- getwd()
  }
  ### Retrieving and validation of the FASTQ file(-s) extension.
  frt <- tools::file_ext(x = fastq)
  if (!frt %in% c("fastq", "gz")) {
    stop("Invalid file format")
  }
  ### Loading the FASTQ data.
  reads <- readFASTQ(
    fastqDir = fastqDir, fastq = fastq, n = n, workDir = workDir
  )[[1]]
  ### Loading the foreign sequences.
  FSs <- read.table(
    file = paste(workDir, contaminants, sep = "/"),
    sep = "\t", header = TRUE, quote = "\"", as.is = TRUE
  )
  ### Computing the hits.
  hits <- lapply(
    seq_len(nrow(x = FSs)),
    function(i) .computeAdapterHits(
      reads, FSs$foreign_sequence[i]
    )
  )
  ### Building the output data.frame.
  positions <- .buildPositionBins(reads)
  hits <- data.frame(cbind(positions, do.call(what = cbind, args = hits)))
  colnames(x = hits) <- c("positions", FSs$foreign_name)
  hits[, -1] <- apply(X = hits[, -1], MARGIN = 2, FUN = as.numeric)
  ### Returning the final object.
  return(hits)
}


#' Collect and present the results of the quality assessment of short reads
#' @description A high-level function for consolidating and presenting the
#'     results of quality assessment of short reads generated by the function
#'     assessQRawReads().
#' @param x character string specifying the name of object developed by the
#'     function assessQRawReads().
#' @param output character string that specifying a postfix for naming all
#'     output report files. The default value is NULL, at which the name of
#'     sample will be used for naming all output files.
#' @param workDir character string specifying the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return a series of text and graphic files.
#' @author Vasily V. Grinev
#' @examples
#' report <- reportQAResults(
#'     x = qa,
#'     output = NULL,
#'     workDir = "D:/Vasily Grinev"
#' )
#' @export
#' @importFrom Biostrings DNAStringSet writeXStringSet
#' @importFrom grid grid.draw textGrob gpar gTree gList
#' @importFrom gridExtra tableGrob ttheme_minimal
#' @importFrom vioplot vioplot
#' @importFrom graphics plot axis title rect legend mtext text barplot
#' @importFrom stats ks.test density
#' @importFrom grDevices pdf dev.off
reportQAResults <- function(x,
                            output = NULL,
                            workDir = NULL) {
  ### Textual report.
  .saveQATextReport(x, output, workDir)
  ### Saving of overrepresented reads.
  .saveQAOverrepresented(x, output, workDir)
  ### Graphical report.
  .saveQAGraphicalReport(x, output, workDir)
  ### Returning nothing.
  invisible(NULL)
}

#' Save textual QA report
#' @description Internal helper for reportQAResults(). Builds the summary
#'     table and writes it as a PDF file.
#' @param x list with QA metrics.
#' @param output character string with the output postfix, or NULL.
#' @param workDir character string giving the path to working directory.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.saveQATextReport <- function(x, output, workDir) {
  tab <- .buildQATable(x)
  file_output <- .resolveQAFileName("QA metrics of", x$SN, output, "pdf")
  pdf(
    file = paste(workDir, file_output, sep = "/"),
    height = nrow(x = tab) / 3
  )
  theme <- ttheme_minimal(
    core = list(
      fg_params = list(hjust = 0, x = 0.03, fontsize = 10),
      bg_params = list(
        fill = rep(c("azure1", "azure3"), each = 1),
        col = "black"
      )
    ),
    colhead = list(
      fg_params = list(col = "white"),
      bg_params = list(fill = "darkcyan", col = "black")
    )
  )
  table <- tableGrob(d = tab, rows = NULL, theme = theme)
  title1 <- textGrob(
    label = "Quality of the raw reads. Summary statistics",
    just = "centre", vjust = -60,
    gp = gpar(col = "black", font = 2, fontsize = 12)
  )
  title2 <- textGrob(
    label = sprintf("Sample: %s", x$SN),
    just = "centre", vjust = -63.5,
    gp = gpar(col = "black", fontsize = 11)
  )
  title3 <- textGrob(
    label = sprintf(paste("Date & Time:",
                          format(x = Sys.time(), "%Y-%m-%d %X"), sep = " ")),
    just = "centre", vjust = -67.8,
    gp = gpar(col = "black", fontsize = 10)
  )
  note <- textGrob(
    label = "Generated by assessQRawReads() & reportQAResults()",
    just = "centre", vjust = 88,
    gp = gpar(col = "black", fontsize = 8)
  )
  gt <- gTree(children = gList(table, title1, title2, title3, note))
  grid.draw(gt)
  invisible(dev.off())
  invisible(NULL)
}

#' Build QA summary table
#' @description Internal helper for reportQAResults(). Builds the data.frame
#'     with the 49 quality metrics.
#' @param x list with QA metrics.
#' @return data.frame with two columns: Parameter and Value.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.buildQATable <- function(x) {
  tab <- data.frame(
    Parameter = c(
      "Sample name", "File name",
      "Total number of reads", "Total number of bases",
      "Length of reads, min", "Length of reads, 1st quartile",
      "Length of reads, median", "Length of reads, 3rd quartile",
      "Length of reads, max",
      "Per entire reads quality, min",
      "Per entire reads quality, 1st quartile",
      "Per entire reads quality, median",
      "Per entire reads quality, 3rd quartile",
      "Per entire reads quality, max",
      "Frequency of A, min", "Frequency of A, 1st quartile",
      "Frequency of A, median", "Frequency of A, 3rd quartile",
      "Frequency of A, max",
      "Frequency of C, min", "Frequency of C, 1st quartile",
      "Frequency of C, median", "Frequency of C, 3rd quartile",
      "Frequency of C, max",
      "Frequency of G, min", "Frequency of G, 1st quartile",
      "Frequency of G, median", "Frequency of G, 3rd quartile",
      "Frequency of G, max",
      "Frequency of T, min", "Frequency of T, 1st quartile",
      "Frequency of T, median", "Frequency of T, 3rd quartile",
      "Frequency of T, max",
      "Frequency of N, min", "Frequency of N, 1st quartile",
      "Frequency of N, median", "Frequency of N, 3rd quartile",
      "Frequency of N, max",
      "Frequency of GC, min", "Frequency of GC, 1st quartile",
      "Frequency of GC, median", "Frequency of GC, 3rd quartile",
      "Frequency of GC, max",
      "Low-complexity of reads, min",
      "Low-complexity of reads, 1st quartile",
      "Low-complexity of reads, median",
      "Low-complexity of reads, 3rd quartile",
      "Low-complexity of reads, max"
    ),
    Value = 0
  )
  tab[1, 2] <- x$SN
  tab[2, 2] <- x$FN
  tab[3, 2] <- x$TNRs
  tab[4, 2] <- x$TNBs
  tab[5:9, 2] <- x$RLs[-4]
  tab[10:14, 2] <- sprintf(x = x$PERQs$metrics, fmt = "%.1f")
  tab[15:39, 2] <- formatC(
    x = c(t(x$BCs$metrics[, -1])), format = "e", digits = 2
  )
  tab[40:44, 2] <- formatC(x = x$GCs$metrics, format = "e", digits = 2)
  tab[45:49, 2] <- x$DSs$metrics
  tab
}

#' Save overrepresented reads
#' @description Internal helper for reportQAResults(). Writes the
#'     overrepresented reads to a FASTA file.
#' @param x list with QA metrics.
#' @param output character string with the output postfix, or NULL.
#' @param workDir character string giving the path to working directory.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.saveQAOverrepresented <- function(x, output, workDir) {
  ORs <- DNAStringSet(x = x$ORs$sequence)
  names(x = ORs) <- paste(
    paste("seq",
          formatC(x = seq_along(x = x$ORs$sequence),
                  width = 4, flag = "0"),
          sep = ""),
    x$ORs$count, sep = ", "
  )
  file_output <- .resolveQAFileName(
    "Overrepresented reads of", x$SN, output, "fasta"
  )
  writeXStringSet(
    x = ORs,
    filepath = paste(workDir, file_output, sep = "/"),
    format = "fasta"
  )
  invisible(NULL)
}

#' Resolve QA output file name
#' @description Internal helper for reportQAResults(). Builds the output file
#'     name based on the prefix, sample name or output, and extension.
#' @param prefix character string with the report prefix.
#' @param sample_name character string with the sample name.
#' @param output character string with the output postfix, or NULL.
#' @param ext character string with the file extension.
#' @return character string with the output file name.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.resolveQAFileName <- function(prefix, sample_name, output, ext) {
  if (is.null(x = output)) {
    paste(sprintf("%s %s", prefix, sample_name), ext, sep = ".")
  } else {
    paste(sprintf("%s %s", prefix, output), ext, sep = ".")
  }
}

#' Save graphical QA report
#' @description Internal helper for reportQAResults(). Writes all QA plots
#'     to a single PDF file.
#' @param x list with QA metrics.
#' @param output character string with the output postfix, or NULL.
#' @param workDir character string giving the path to working directory.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.saveQAGraphicalReport <- function(x, output, workDir) {
  file_output <- .resolveQAFileName("QA plots of", x$SN, output, "pdf")
  pdf(
    file = paste(workDir, file_output, sep = "/"),
    onefile = TRUE, paper = "a4", width = 8.27, height = 11.69
  )
  par(mfrow = c(3, 2), omi = c(0.79, 0.79, 1, 0.79))
  .plotReadQualityDistribution(x)
  .plotReadQualityThreshold(x)
  .plotPerCycleQuality(x)
  .plotBaseComposition(x)
  .plotPerCycleBaseComposition(x)
  .plotGcComposition(x)
  .plotLowComplexity(x)
  .plotReadFrequency(x)
  if (!is.null(x = nrow(x = x$ASs))) {
    .plotContaminations(x$ASs, "Adapter contaminations",
                        "Adapter frequency")
  }
  if (!is.null(x = nrow(x = x$FSs))) {
    .plotContaminations(x$FSs, "Foreign sequence contaminations",
                        "Sequence frequency")
  }
  invisible(dev.off())
  invisible(NULL)
}

#' Plot read quality distribution
#' @description Internal helper for reportQAResults(). Draws the violin plot
#'     of per-read quality.
#' @param x list with QA metrics.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.plotReadQualityDistribution <- function(x) {
  plot(0, type = "n", axes = FALSE, ann = FALSE, bty = "n",
       xlim = c(0, 2), ylim = c(0, 40))
  vioplot(
    x = x$PERQs$means[seq_len(10000)], add = TRUE,
    col = rgb(red = 86, green = 180, blue = 233, maxColorValue = 255),
    rectCol = rgb(red = 0, green = 158, blue = 115, maxColorValue = 255),
    pchMed = 0, colMed = NA, lwd = 0.9, wex = 1, frame.plot = FALSE
  )
  segments(
    x0 = 0.9, y0 = median(x = x$PERQs$means),
    x1 = 1.1, y1 = median(x$PERQs$means),
    col = rgb(red = 240, green = 228, blue = 66, maxColorValue = 255),
    lwd = 3, lend = 2
  )
  axis(side = 1, at = c(0, 1, 2), labels = FALSE, lwd = 0.8)
  axis(side = 2, at = c(0, 10, 20, 30, 40),
       labels = c(0, 10, 20, 30, 40), las = 1, lwd = 0.8)
  title(main = "Distribution of reads quality", cex.main = 1.5)
  title(xlab = "Standard vioplot", line = 1, cex.lab = 1.2)
  title(ylab = "Average read quality scores", line = 2.3, cex.lab = 1.2)
  invisible(NULL)
}

#' Plot read quality threshold
#' @description Internal helper for reportQAResults(). Draws the threshold
#'     plot of per-read quality.
#' @param x list with QA metrics.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.plotReadQualityThreshold <- function(x) {
  per <- data.frame(
    threshold = 0:ceiling(x = max(x = x$PERQs$means)),
    percent = 0
  )
  per[1, 2] <- 100
  L <- length(x = x$PERQs$means)
  for (i in seq_len(ceiling(x = max(x = x$PERQs$means)))) {
    p <- length(x = x$PERQs$means[x$PERQs$means > i]) / L * 100
    per[i + 1, 2] <- round(x = p, digits = 1)
  }
  plot(per, type = "l", bty = "n",
       xlim = c(0, 40), ylim = c(0, 100), lwd = 1,
       col = rgb(red = 0, green = 114, blue = 178, maxColorValue = 255),
       ann = FALSE, axes = FALSE)
  axis(side = 1, at = c(0, 10, 20, 30, 40),
       labels = c(0, 10, 20, 30, 40), las = 1, lwd = 0.8)
  axis(side = 2, at = c(0, 20, 40, 60, 80, 100),
       labels = c(0, 20, 40, 60, 80, 100), las = 1, lwd = 0.8)
  title(main = "Thresholding of reads quality", cex.main = 1.5)
  title(xlab = "Quality threshold", line = 2.5, cex.lab = 1.2)
  title(ylab = "Reads exceeding quality level, %", line = 2.8, cex.lab = 1.2)
  rect(xleft = 0, ybottom = 0, xright = 20, ytop = 100, border = NA,
       col = rgb(240, 228, 66, maxColorValue = 255, alpha = 40))
  rect(xleft = 20, ybottom = 0, xright = 30, ytop = 100, border = NA,
       col = rgb(230, 159, 0, maxColorValue = 255, alpha = 40))
  rect(xleft = 30, ybottom = 0, xright = 40, ytop = 100, border = NA,
       col = rgb(213, 94, 0, maxColorValue = 255, alpha = 70))
  text(x = 20, y = 20, srt = 90,
       labels = paste0(
         "Q20, ",
         round(length(x$PERQs$means[x$PERQs$means > 20]) / L * 100, 1),
         "%"
       ))
  text(x = 30, y = 20, srt = 90,
       labels = paste0(
         "Q30, ",
         round(length(x$PERQs$means[x$PERQs$means > 30]) / L * 100, 1),
         "%"
       ))
  invisible(NULL)
}

#' Plot per-cycle quality
#' @description Internal helper for reportQAResults(). Draws the boxplot of
#'     per-cycle quality.
#' @param x list with QA metrics.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.plotPerCycleQuality <- function(x) {
  bxp(
    z = list(
      stats = t(x$PCQRs[, -1]),
      n = rep(x = length(x = x$PERQs$means),
              times = nrow(x = x$PCQRs))
    ),
    ylim = c(0, 40),
    boxfill = rgb(240, 228, 66, maxColorValue = 255),
    frame = FALSE, axes = FALSE, lwd = 0.6,
    medlwd = 1, medcol = rgb(204, 121, 167, maxColorValue = 255),
    whisklty = "solid"
  )
  axis(side = 1, at = seq_along(x = x$PCQRs[, 1]),
       labels = FALSE, las = 2, lwd = 0.8)
  axis(side = 2, at = c(0, 10, 20, 30, 40),
       labels = c(0, 10, 20, 30, 40), las = 1, lwd = 0.8)
  text(
    x = seq(from = 1, to = length(x = x$PCQRs[, 1]), by = 2),
    par("usr")[3] - 3, offset = -0.5, cex = 0.9,
    labels = x$PCQRs[, 1][seq(from = 1,
                              to = length(x = x$PCQRs[, 1]), by = 2)],
    srt = 55, pos = 2, xpd = TRUE
  )
  title(main = "Per cycle quality of reads", cex.main = 1.5)
  title(xlab = "Cycles", line = 3.7, cex.lab = 1.2)
  title(ylab = "Quality scores across of all bases", line = 2.3,
        cex.lab = 1.2)
  rect(xleft = 0.5, ybottom = 0,
       xright = length(x = x$PCQRs[, 1]) + 0.5, ytop = 20, border = NA,
       col = rgb(240, 228, 66, maxColorValue = 255, alpha = 40))
  rect(xleft = 0.5, ybottom = 20,
       xright = length(x = x$PCQRs[, 1]) + 0.5, ytop = 30, border = NA,
       col = rgb(230, 159, 0, maxColorValue = 255, alpha = 40))
  rect(xleft = 0.5, ybottom = 30,
       xright = length(x = x$PCQRs[, 1]) + 0.5, ytop = 40, border = NA,
       col = rgb(213, 94, 0, maxColorValue = 255, alpha = 70))
  invisible(NULL)
}

#' Plot base composition
#' @description Internal helper for reportQAResults(). Draws the boxplot of
#'     per-read base composition.
#' @param x list with QA metrics.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.plotBaseComposition <- function(x) {
  colrs <- c(
    rgb(0, 114, 178, maxColorValue = 255),
    rgb(213, 94, 0, maxColorValue = 255),
    rgb(86, 180, 233, maxColorValue = 255),
    rgb(230, 159, 0, maxColorValue = 255),
    rgb(204, 121, 167, maxColorValue = 255)
  )
  bxp(
    z = list(
      stats = t(x$BCs$metrics[, -1]),
      n = rep(x = length(x = x$PERQs$means),
              times = nrow(x = x$BCs$metrics))
    ),
    ylim = c(0, 0.5), boxfill = colrs, frame = FALSE, axes = FALSE,
    lwd = 0.6, medlwd = 1, boxwex = 0.6, whisklty = "solid"
  )
  axis(side = 1, at = c(1, 2, 3, 4, 5),
       labels = x$BCs$metrics[, 1], las = 1, lwd = 0.8)
  axis(side = 2, at = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5),
       labels = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5), las = 2, lwd = 0.8)
  title(main = "Base composition of reads", cex.main = 1.5)
  title(xlab = "Base", line = 2.5, cex.lab = 1.2)
  title(ylab = "Frequency", line = 2.8, cex.lab = 1.2)
  invisible(NULL)
}

#' Plot per-cycle base composition
#' @description Internal helper for reportQAResults(). Draws the per-cycle
#'     base composition lines.
#' @param x list with QA metrics.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.plotPerCycleBaseComposition <- function(x) {
  colrs <- c(
    rgb(0, 114, 178, maxColorValue = 255),
    rgb(213, 94, 0, maxColorValue = 255),
    rgb(86, 180, 233, maxColorValue = 255),
    rgb(230, 159, 0, maxColorValue = 255),
    rgb(204, 121, 167, maxColorValue = 255)
  )
  matplot(x$PCCBs[, -1], type = "l", lty = 1, lwd = 1, col = colrs,
          ylim = c(0, 1), ann = FALSE, axes = FALSE)
  axis(side = 1, at = seq_along(x = x$PCCBs[, 1]),
       labels = FALSE, las = 2, lwd = 0.8)
  axis(side = 2, at = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
       labels = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), las = 1, lwd = 0.8)
  text(
    x = seq(from = 1, to = length(x = x$PCCBs[, 1]), by = 2),
    par("usr")[3] - 0.08, offset = -0.5, cex = 0.9,
    labels = x$PCCBs[, 1][seq(from = 1,
                              to = length(x = x$PCCBs[, 1]), by = 2)],
    srt = 55, pos = 2, xpd = TRUE
  )
  title(main = "Per cycle composition of reads", cex.main = 1.5)
  title(xlab = "Cycles", line = 3.8, cex.lab = 1.2)
  title(ylab = "Base frequency", line = 2.3, cex.lab = 1.2)
  legend(x = "topleft", legend = colnames(x = x$PCCBs)[-1],
         col = colrs, lty = 1, lwd = 1, bty = "n")
  invisible(NULL)
}

#' Plot GC composition
#' @description Internal helper for reportQAResults(). Draws the GC density
#'     plot with KS-test annotation.
#' @param x list with QA metrics.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.plotGcComposition <- function(x) {
  plot(x$GCs$density, bty = "n", type = "l", lwd = 1, xlim = c(0, 1),
       ann = FALSE, axes = FALSE,
       col = rgb(0, 114, 178, maxColorValue = 255))
  lines(density(x = rnorm(
    n = length(x = x$GCs$frequency),
    mean = mean(x = x$GCs$frequency),
    sd = sd(x = x$GCs$frequency)
  ))[seq_len(2)], xlim = c(0, 1), lwd = 1,
  col = rgb(213, 94, 0, maxColorValue = 255))
  axis(side = 1, at = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
       labels = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), lwd = 0.8)
  axis(side = 2,
       at = c(0.0, max(x = x$GCs$density[, 2]) / 4 * seq_len(4)),
       labels = round(c(0.0,
                        max(x = x$GCs$density[, 2]) / 4 * seq_len(4)),
                      digits = 2),
       las = 1, lwd = 0.8)
  title(main = "GC composition of reads", cex.main = 1.5)
  title(xlab = "GC frequency", line = 2.5, cex.lab = 1.2)
  title(ylab = "Density", line = 2.8, cex.lab = 1.2)
  legend(x = "topleft", legend = c("theoretical", "empirical"),
         col = c(rgb(213, 94, 0, maxColorValue = 255),
                 rgb(0, 114, 178, maxColorValue = 255)),
         lty = 1, lwd = 1, bty = "n")
  ks_p <- tryCatch(
    ks.test(x = x$GCs$frequency, y = "pnorm")$p.value,
    error = function(e) NA_real_,
    warning = function(w) NA_real_
  )
  legend(x = "topright",
         legend = paste("KS test: p = ", round(ks_p, 4), sep = ""),
         lty = 0, bty = "n")
  invisible(NULL)
}

#' Plot low-complexity
#' @description Internal helper for reportQAResults(). Draws the dusty-score
#'     density plot.
#' @param x list with QA metrics.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.plotLowComplexity <- function(x) {
  plot(x$DSs$density, bty = "n", type = "l", lwd = 1, xlim = c(2, 5),
       ann = FALSE, axes = FALSE,
       col = rgb(0, 114, 178, maxColorValue = 255))
  axis(side = 1, at = c(2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0),
       labels = c(2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0), lwd = 0.8)
  axis(side = 2,
       at = c(0.0, max(x = x$DSs$density[, 2]) / 4 * seq_len(4)),
       labels = round(c(0.0,
                        max(x = x$DSs$density[, 2]) / 4 * seq_len(4)),
                      digits = 2),
       las = 1, lwd = 0.8)
  title(main = "Low-complexity of reads", cex.main = 1.5)
  title(xlab = expression("Complexity, log"[10] * " (dusty scores)"),
        line = 2.75, cex.lab = 1.2)
  title(ylab = "Density", line = 3, cex.lab = 1.2)
  legend(x = 3.3, y = max(x = x$DSs$density[, 2]),
         legend = c("minimum:", "Q1:", "median:", "Q3:", "maximum:"),
         lty = 0, bty = "n")
  legend(x = 4.1, y = max(x = x$DSs$density[, 2]),
         legend = round(x = log10(x = x$DSs$metrics), digits = 2),
         lty = 0, bty = "n")
  rect(xleft = 2.0, ybottom = 0, xright = 3.0,
       ytop = max(x = x$DSs$density[, 2]),
       density = 20, border = NA, lwd = 0.3,
       col = rgb(0, 158, 115, maxColorValue = 255))
  invisible(NULL)
}

#' Plot read frequency
#' @description Internal helper for reportQAResults(). Draws the barplot of
#'     read occurrence frequency.
#' @param x list with QA metrics.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.plotReadFrequency <- function(x) {
  barplot(
    formula = x$FRs[, 3] ~ x$FRs[, 1], data = x$FRs,
    axes = FALSE, ann = FALSE, ylim = c(0, 1), border = NA,
    col = colorRampPalette(c(
      rgb(213, 94, 0, maxColorValue = 255),
      rgb(230, 159, 0, maxColorValue = 255)
    ))(nrow(x = x$FRs))
  )
  axis(side = 2, at = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
       labels = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), las = 1, lwd = 0.8)
  title(main = "Frequency of reads", cex.main = 1.5)
  title(xlab = "Occurrence of reads", line = 2.5, cex.lab = 1.2)
  title(ylab = "Frequency", line = 3, cex.lab = 1.2)
  invisible(NULL)
}

#' Plot adapter or foreign sequence contaminations
#' @description Internal helper for reportQAResults(). Draws the
#'     contamination frequency lines.
#' @param data data.frame with positions and one or more sequence columns.
#' @param main character string with the plot title.
#' @param ylab character string with the y-axis label.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.plotContaminations <- function(data, main, ylab) {
  colrs <- c(
    "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288",
    "#AA4499", "#44AA99", "#999933", "#882255", "#661100",
    "#6699CC", "#888888"
  )
  matplot(data[, -1], ann = FALSE, axes = FALSE, ylim = c(0, 0.35),
          type = "l", lty = 1, lwd = 1,
          col = colrs[seq_len(ncol(x = data[, -1]))])
  axis(side = 1, at = seq_along(x = data[, 1]), labels = FALSE,
       lwd = 0.8)
  axis(side = 2,
       at = c(0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35),
       labels = c(0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35),
       las = 1, lwd = 0.8)
  text(x = seq(from = 1, to = length(x = data[, 1]), by = 2),
       par("usr")[3] - 0.028, offset = -0.5,
       labels = data[, 1][seq(from = 1,
                              to = length(x = data[, 1]), by = 2)],
       cex = 0.9, srt = 55, pos = 2, xpd = TRUE)
  title(main = main, cex.main = 1.5)
  title(xlab = "Cycles", line = 3.8, cex.lab = 1.2)
  title(ylab = ylab, line = 3, cex.lab = 1.2)
  legend(x = "topleft", legend = colnames(x = data)[-1],
         cex = 0.8, lty = 1, lwd = 1,
         col = colrs[seq_len(ncol(x = data[, -1]))],
         bty = "n")
  invisible(NULL)
}
