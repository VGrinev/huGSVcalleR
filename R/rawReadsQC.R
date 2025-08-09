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
#' qa <- assessQRawReads(fastqDir="Files_FASTQ",
#'                       fastq="2-Galkina-A_S2_L001_R1_001.fastq.gz",
#'                       n=1e5,
#'                       adapters=NULL,
#'                       contaminants=NULL,
#'                       workDir="D:/Vasily Grinev")
#' @export
#' @importFrom Biostrings DNAStringSet PhredQuality width letterFrequency
#' @importFrom ShortRead countFastq qa dustyScore
#' @importFrom data.table data.table rbindlist
#' @importFrom IRanges IntegerList
#' @importFrom stats density quantile
#' @importFrom methods as is
#' @importFrom stats density quantile
#' @importFrom DelayedArray mean
assessQRawReads <- function(fastqDir=NULL,
                            fastq,
                            n=NULL,
                            adapters=NULL,
                            contaminants=NULL,
                            workDir=NULL){
  ### Setting the working directory.
  if (is.null(x=workDir)){
    workDir <- getwd()
  }
  ### Full path to the FASTQ file(-s).
  if (is.null(x=fastqDir)){
    path <- paste(workDir, fastqDir, sep="")
  }else{
    path <- paste(workDir, fastqDir, sep="/")
  }
  ### Retrieving and validation of the FASTQ file(-s) extension.
  frt <- tools::file_ext(x=fastq)
  if (!frt %in% c("fastq", "gz")){
    stop("Invalid file format")
  }
  ### Collection of the quality metrics.
  ##  Sample name.
  sample_name <- gsub(pattern=".fastq.+", replacement="", x=fastq)
  ##  File name.
  file_name <- fastq
  ##  Total number of reads.
  TNBs <- countFastq(dirPath=path, pattern=fastq)
  TNRs <- TNBs[[1]]
  ##  Total number of bases.
  TNBs <- TNBs[[2]]
  ##  Distribution the length of reads.
  reads <- readFASTQ(fastqDir=fastqDir, fastq=fastq, n=n, workDir=workDir)
  RLs <- as.numeric(x=summary(object=width(x=reads$reads)))
  ##  Per entire reads quality.
  PERQs <- mean(x=as(object=reads$QScores, Class="IntegerList"))
  names(x=PERQs) <- NULL
  PERQs <- list(means=PERQs,
                density=do.call(what=cbind, args=density(x=PERQs)[1:2]),
                metrics=quantile(x=PERQs,
                                 probs=c(0.1, 0.25, 0.5, 0.75, 0.9)))
  ##  Per cycle quality of reads.
  if (is.null(x=n)){
    QA <- qa(dirPath=path, pattern=fastq)
  }else{
    QA <- qa(dirPath=path, pattern=fastq, n=n)
  }
  PCQRs <- do.call(what=rbind,
                   args=lapply(X=split(x=QA[["perCycle"]]$quality,
                                       f=QA[["perCycle"]]$quality$Cycle),
                               FUN=function(y){rep(x=y[, 3], times=y[, 4])}))
  m <- list()
  for (i in seq(from=1, to=nrow(x=PCQRs), by=5)){
    if (i + 4 < nrow(x=PCQRs)){
      m[[i]] <- c(paste(i, i:i + 4, sep="-"),
                  quantile(x=PCQRs[i:i + 4, ],
                           probs=c(0.1, 0.25, 0.5, 0.75, 0.9)))
    }else{
      m[[i]] <- c(paste(i, nrow(x=PCQRs), sep="-"),
                  quantile(x=PCQRs[i:nrow(x=PCQRs), ],
                           probs=c(0.1, 0.25, 0.5, 0.75, 0.9)))

    }
  }
  PCQRs <- data.frame(do.call(args=m, what=rbind))
  colnames(x=PCQRs) <- c("cycles", "minimum", "Q1", "median", "Q3", "maximum")
  PCQRs[, -1] <- apply(X=PCQRs[, -1], MARGIN=2, FUN=as.numeric)
  ##  Bases composition of the entire set of reads.
  BCs <- letterFrequency(x=reads$reads,
                         letters=c("A", "C", "G", "T", "N"))
  BCs <- BCs/rowSums(x=BCs)
  metrics=rbind(quantile(x=BCs[, 1], probs=c(0.1, 0.25, 0.5, 0.75, 0.9)),
                quantile(x=BCs[, 2], probs=c(0.1, 0.25, 0.5, 0.75, 0.9)),
                quantile(x=BCs[, 3], probs=c(0.1, 0.25, 0.5, 0.75, 0.9)),
                quantile(x=BCs[, 4], probs=c(0.1, 0.25, 0.5, 0.75, 0.9)),
                quantile(x=BCs[, 5], probs=c(0.1, 0.25, 0.5, 0.75, 0.9)))
  metrics <- data.frame(cbind(colnames(x=BCs), metrics))
  colnames(x=metrics) <- c("nucleotide",
                           "minimum", "Q1", "median", "Q3", "maximum")
  metrics[, -1] <- apply(X=metrics[, -1], MARGIN=2, FUN=as.numeric)
  BCs <- list(frequency=BCs, metrics=metrics)
  ##  Per cycle content of bases.
  PCCBs <- QA[["perCycle"]]$baseCall[, -4]
  PCCBs$Base <- as.character(x=PCCBs$Base)
  PCCBs <- do.call(what=rbind,
                   args=lapply(X=split(x=PCCBs,
                                       f=PCCBs$Cycle),
                               FUN=function(y){rep(x=y[, 2], times=y[, 3])}))
  m <- data.frame(cycles="", A=0, C=0, G=0, T=0, N=0)
  for (i in seq(from=1, to=nrow(x=PCCBs), by=5)){
    if (i + 4 < nrow(x=PCCBs)){
      bases <- c(paste(i, i:i + 4, sep="-"),
                 table(PCCBs[i:i + 4, ]))
      m[i, 1] <- bases[1]
      m[i, -1] <- bases[match(colnames(x=m), names(x=bases))][-1]
    }else{
      bases <- c(paste(i, nrow(x=PCCBs), sep="-"),
                 table(PCCBs[i:nrow(x=PCCBs), ]))
      m[i, 1] <- bases[1]
      m[i, -1] <- bases[match(colnames(x=m), names(x=bases))][-1]
    }
  }
  PCCBs <- m[!is.na(x=m$cycles), ]
  PCCBs[is.na(x=PCCBs$N), ]$N <- 0
  PCCBs[, -1] <- apply(X=PCCBs[, -1], MARGIN=2, FUN=as.numeric)
  rownames(x=PCCBs) <- NULL
  PCCBs[, -1] <- PCCBs[, -1]/rowSums(x=PCCBs[, -1])
  ##  GC composition of the entire set of reads.
  GCs <- letterFrequency(x=reads$reads, letters=c("G", "C"))
  GCs <- rowSums(x=GCs)/width(x=reads$reads)
  GCs <- list(frequency=GCs,
              density=do.call(what=cbind, args=density(x=GCs)[1:2]),
              metrics=quantile(x=GCs, probs=c(0.1, 0.25, 0.5, 0.75, 0.9)))
  ##  Summarization the low-complexity of reads (using dusty scores).
  DSs <- dustyScore(x=reads$reads)
  DSs <- list(scores=DSs,
              density=do.call(what=cbind, args=density(x=log10(x=DSs))[1:2]),
              metrics=quantile(x=DSs, probs=c(0.1, 0.25, 0.5, 0.75, 0.9)))
  ##  Frequency of reads.
  FRs <- QA[["sequenceDistribution"]][, 1:2]
  colnames(x=FRs) <- c("occurrence", "reads")
  FRs$frequency <- FRs$reads/sum(x=FRs$reads)
  ##  Overrepresented reads.
  ORs <- ShortRead:::.freqSequences(QA, "read")[, 1:2]
  ##  Identification of irrelevant sequences.
  #   Identification of adapters.
  ASs <- NA
  if (!is.null(x=adapters)){
    ASs <- findAdapters(fastqDir=fastqDir,
                        fastq=fastq,
                        n=n,
                        adapters=adapters,
                        workDir=workDir)
  }
  #   Identification of foreign sequences.
  FSs <- NA
  if (!is.null(x=contaminants)){
    FSs <- findForeignSeqs(fastqDir=fastqDir,
                           fastq=fastq,
                           n=n,
                           contaminants=contaminants,
                           workDir=workDir)
  }
  ### Quality metrics.
  return(list(SN=sample_name,
              FN=file_name,
              TNRs=TNRs,
              TNBs=TNBs,
              RLs=RLs,
              PERQs=PERQs,
              PCQRs=PCQRs,
              BCs=BCs,
              PCCBs=PCCBs,
              GCs=GCs,
              DSs=DSs,
              FRs=FRs,
              ORs=ORs,
              ASs=ASs,
              FSs=FSs)
  )
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
#' fq <- readFASTQ(fastqDir="Files_FASTQ",
#'                 fastq="example_seq.read1.fastq.gz",
#'                 n=1e5,
#'                 workDir="D:/Vasily Grinev")
#' @export
#' @importFrom Biostrings DNAStringSet PhredQuality
#' @importFrom tools file_ext
readFASTQ <- function(fastqDir=NULL,
                      fastq,
                      n=NULL,
                      workDir=NULL){
  ### Setting the working directory.
  if (is.null(x=workDir)){
    workDir <- getwd()
  }
  ### Full path to the FASTQ file(-s).
  if (is.null(x=fastqDir)){
    path <- paste(workDir, fastqDir, sep="")
  }else{
    path <- paste(workDir, fastqDir, sep="/")
  }
  ### Retrieving and validation of the FASTQ file(-s) extension.
  frt <- tools::file_ext(x=fastq)
  if (!frt %in% c("fastq", "gz")){
    stop("Invalid file format")
  }
  ### Loading the FASTQ data.
  if (is.null(x=n)){
    reads <- readLines(con=paste(path, fastq, sep="/"))
    FQs <- list(reads=DNAStringSet(x=reads[seq(from=2,
                                               to=length(x=reads),
                                               by=4)]),
                QScores=PhredQuality(x=reads[seq(from=4,
                                                 to=length(x=reads),
                                                 by=4)]))
    names(x=FQs$reads) <- reads[seq(from=1, to=length(x=reads), by=4)]
    names(x=FQs$QScores) <- reads[seq(from=1, to=length(x=reads), by=4)]
  }else{
    reads <- readLines(con=paste(path, fastq, sep="/"), n=n * 4)
    FQs <- list(reads=DNAStringSet(x=reads[seq(from=2, to=n * 4, by=4)]),
                QScores=PhredQuality(x=reads[seq(from=4, to=n * 4, by=4)]))
    names(x=FQs$reads) <- reads[seq(from=1, to=n * 4, by=4)]
    names(x=FQs$QScores) <- reads[seq(from=1, to=n * 4, by=4)]
  }
  rm(reads)
  ### Returning the final object.
  return(FQs)
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
#' @param dustScore integer specifying the maximum dust score at which reads to
#'     be saved. By default this argument is NULL and such filter is not used.
#' @param batchSize NA (by default) or integer value indicating the maximum
#'     number of reads to be processed at any one time during filtration by
#'     dustScore filter. By default, all reads are processed simultaneously
#'     Smaller values use less memory but are computationally less efficient.
#' @param postfix character string that will be added to the output file name
#'     as postfix. "filtered" by default.
#' @param workDir character string specifying the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return filtered FASTQ file(-s) will be written to the same directory as
#'     input FASTQ file(-s).
#' @author Vasily V. Grinev.
#' @examples
#' cleanRawReads(fastqDir="Files_FASTQ",
#'               fastq1="example_seq.read1.fastq",
#'               fastq2="example_seq.read2.fastq",
#'               adapters="adapters.txt",
#'               error=0.2,
#'               min_match_flank=3L,
#'               anchored=TRUE,
#'               indels=FALSE,
#'               tr_score=20,
#'               tr_start=10,
#'               tr_end=30,
#'               phred="phred_scores.txt",
#'               k=3,
#'               halfwidth=5,
#'               successive=TRUE,
#'               readl=100,
#'               readq=20,
#'               readn=3,
#'               dustScore=875,
#'               batchSize=NA,
#'               postfix="filtered",
#'               workDir="D:/Vasily Grinev")
#' @export
#' @importFrom ShortRead readFastq writeFastq trimEnds trimTailw trimTails id
#' @importFrom ShortRead polynFilter occurrenceFilter dustyFilter
#' @importFrom Biostrings width
#' @importFrom methods as
#' @importFrom utils read.table
cleanRawReads <- function(fastqDir=NULL,
                          fastq1,
                          fastq2=NULL,
                          adapters=NULL,
                          error=0.2,
                          min_match_flank=3L,
                          anchored=TRUE,
                          indels=FALSE,
                          tr_score=NULL,
                          tr_start=10,
                          tr_end=30,
                          phred=NA,
                          k=3,
                          halfwidth=NULL,
                          successive=TRUE,
                          readl=0,
                          readq=0,
                          readn=0,
                          dustScore=NULL,
                          batchSize=NA,
                          postfix="filtered",
                          workDir=NULL){
  ### Setting the working directory.
  if (is.null(x=workDir)){
    workDir <- getwd()
  }
  ### Full path to the FASTQ file(-s).
  if (is.null(x=fastqDir)){
    path <- workDir
  }else{
    path <- paste(workDir, fastqDir, sep="/")
  }
  ### Retrieving and validation of the FASTQ file(-s) extension.
  frt1 <- tools::file_ext(x=fastq1)
  if (!frt1 %in% c("fastq", "gz")){
    stop("Invalid format of FASTQ file 1")
  }
  if (!is.null(x=fastq2)){
    frt2 <- tools::file_ext(x=fastq2)
    if (!frt2 %in% c("fastq", "gz")){
      stop("Invalid format of FASTQ file 2")
    }
  }
  ### Loading of the FASTQ file 1.
  FQ1 <- readFastq(dirPath=path, pattern=fastq1, withIds=TRUE)
  ### Raw reads processing.
  ##  Step 1, trimming of the adapter sequences.
  if (!is.null(x=adapters)){
    FQ1 <- cutLeftAdapter(fastq=FQ1,
                          adapters=adapters,
                          error=error,
                          min_match_flank=min_match_flank,
                          anchored=anchored,
                          indels=indels,
                          workDir=workDir)
    FQ1 <- cutRigthAdapter(fastq=FQ1,
                           adapters=adapters,
                           error=error,
                           min_match_flank=min_match_flank,
                           anchored=anchored,
                           indels=indels,
                           workDir=workDir)
  }
  ##  Step 2, trimming of the low quality ends.
  #   Fixed trimming.
  if (is.null(x=tr_score)){
    FQ1 <- narrow(x=FQ1,
                  start=tr_start + 1,
                  end=max(x=width(x=FQ1)) - tr_end)
  }else{
    #   Dynamic trimming.
    PhSs <- read.table(file=paste(workDir, phred, sep="/"),
                       sep="\t",
                       header=TRUE,
                       quote="",
                       as.is=TRUE,
                       comment.char="z")
    FQ1 <- trimEnds(object=FQ1,
                    a=PhSs[PhSs$score == tr_score, ]$ASCII,
                    left=TRUE, right=FALSE,
                    relation="<=")
    if (!is.null(x=halfwidth)){
      FQ1 <- trimTailw(object=FQ1,
                       k=k,
                       a=PhSs[PhSs$score == tr_score, ]$ASCII,
                       halfwidth=halfwidth)
    }else{
      FQ1 <- trimTails(object=FQ1,
                       k=k,
                       a=PhSs[PhSs$score == tr_score, ]$ASCII,
                       successive=successive)
    }
  }
  ##  Step 3, removing of the too short reads.
  FQ1 <- FQ1[width(x=FQ1) >= readl]
  ##  Step 4, removing of the too low quality reads.
  FQ1 <- FQ1[rowMeans(x=as(quality(x=FQ1), "matrix"), na.rm=TRUE) >= readq]
  ##  Step 5, removing of the reads with ambiguous nucleotides.
  FQ1 <- FQ1[polynFilter(threshold=as.integer(x=readn), nuc="other")(FQ1)]
  ##  Step 6, removing of the duplicated reads.
  FQ1 <- FQ1[occurrenceFilter(min=1L, max=1L,
                              withSread=TRUE,
                              duplicates="sample")(FQ1)]
  if (!is.null(x=dustScore)){
    ##  Step 7, removing of the too low-complexity reads.
    FQ1 <- FQ1[dustyFilter(threshold=dustScore, batchSize=batchSize)(FQ1)]
  }
  if (!is.null(x=fastq2)){
    ### Loading of the FASTQ file 2.
    FQ2 <- readFastq(dirPath=path, pattern=fastq2, withIds=TRUE)
    ### Raw reads processing.
    ##  Step 1, trimming of the adapter sequences.
    if (!is.null(x=adapters)){
      FQ2 <- cutLeftAdapter(fastq=FQ2,
                            adapters=adapters,
                            error=error,
                            min_match_flank=min_match_flank,
                            anchored=anchored,
                            indels=indels,
                            workDir=workDir)
      FQ2 <- cutRigthAdapter(fastq=FQ2,
                             adapters=adapters,
                             error=error,
                             min_match_flank=min_match_flank,
                             anchored=anchored,
                             indels=indels,
                             workDir=workDir)
    }
    ##  Step 2, trimming of the low quality ends.
    #   Fixed trimming.
    if (is.null(x=tr_score)){
      FQ2 <- narrow(x=FQ2,
                    start=tr_start + 1,
                    end=max(x=width(x=FQ2)) - tr_end)
    }else{
      #   Dynamic trimming.
      FQ2 <- trimEnds(object=FQ2,
                      a=PhSs[PhSs$score == tr_score, ]$ASCII,
                      left=TRUE, right=FALSE,
                      relation="<=")
      if (!is.null(x=halfwidth)){
        FQ2 <- trimTailw(object=FQ2,
                         k=k,
                         a=PhSs[PhSs$score == tr_score, ]$ASCII,
                         halfwidth=halfwidth)
      }else{
        FQ2 <- trimTails(object=FQ2,
                         k=k,
                         a=PhSs[PhSs$score == tr_score, ]$ASCII,
                         successive=successive)
      }
    }
    ##  Step 3, removing of the too short reads.
    FQ2 <- FQ2[width(x=FQ2) >= readl]
    ##  Step 4, removing of the too low quality reads.
    FQ2 <- FQ2[rowMeans(x=as(quality(x=FQ2), "matrix"),
                        na.rm=TRUE) >= readq]
    ##  Step 5, removing of the reads with ambiguous nucleotides.
    FQ2 <- FQ2[polynFilter(threshold=as.integer(x=readn),
                           nuc="other")(FQ2)]
    ##  Step 6, removing of the duplicated reads.
    FQ2 <- FQ2[occurrenceFilter(min=1L, max=1L,
                                withSread=TRUE,
                                duplicates="sample")(FQ2)]
    if (!is.null(x=dustScore)){
      ##  Step 7, removing of the too low-complexity reads.
      FQ2 <- FQ2[dustyFilter(threshold=dustScore,
                             batchSize=batchSize)(FQ2)]
    }
    ID1 <- gsub(pattern=" 1.+", replacement="", x=id(object=FQ1))
    ID2 <- gsub(pattern=" 2.+", replacement="", x=id(object=FQ2))
    FQ1 <- FQ1[ID1 %in% ID2]
    FQ2 <- FQ2[ID2 %in% ID1]
  }
  ### Writing the filtered FASTQ file(-s).
  output1 <- paste(path,
                   gsub(pattern="fast.+",
                        replacement=paste(postfix, "fastq.gz", sep="."),
                        x=fastq1),
                   sep="/")
  writeFastq(object=FQ1, file=output1, compress=TRUE)
  if (!is.null(x=fastq2)){
    output2 <- paste(path,
                     gsub(pattern="fast.+",
                          replacement=paste(postfix, "fastq.gz", sep="."),
                          x=fastq2),
                     sep="/")
    writeFastq(object=FQ2, file=output2)
  }
  rm(list=ls())
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
#' trimL <- cutLeftAdapter(fastqDir="Files_FASTQ",
#'                         fastq="example_seq.read1.fastq",
#'                         adapters="adapters.txt",
#'                         error=0.2,
#'                         min_match_flank=3L,
#'                         anchored=TRUE,
#'                         indels=FALSE,
#'                         workDir="D:/Vasily Grinev")
#' @export
#' @importFrom ShortRead readFastq
#' @importFrom Biostrings trimLRPatterns width
#' @importFrom utils read.table
#' @importFrom methods is
cutLeftAdapter <- function(fastqDir=NULL,
                           fastq,
                           adapters,
                           error=0.2,
                           min_match_flank=3L,
                           anchored=TRUE,
                           indels=FALSE,
                           workDir=NULL){
  ### Setting the working directory.
  if (is.null(x=workDir)){
    workDir <- getwd()
  }
  ### Retrieving of the FASTQ data.
  if (is.object(x=fastq)){
    FQ=fastq
  }else{
    ### Full path to the FASTQ file.
    if (is.null(x=fastqDir)){
      path <- paste(workDir, fastqDir, sep="")
    }else{
      path <- paste(workDir, fastqDir, sep="/")
    }
    ### Retrieving and validation of the FASTQ file extension.
    frt <- tools::file_ext(x=fastq)
    if (!frt %in% c("fastq", "gz")){
      stop("Invalid format of FASTQ file")
    }
    FQ <- readFastq(dirPath=path, pattern=fastq, withIds=TRUE)
  }
  ### Loading of the adapter sequence(-s).
  ADs <- read.table(file=paste(workDir, adapters, sep="/"),
                    sep="\t",
                    header=TRUE,
                    quote="\"",
                    as.is=TRUE)
  ### Left-end trimming of adapter sequence(-s).
  for (i in 1:nrow(x=ADs)){
    seq_l <- ADs[i, 3]
    if (error > 0){
      flank_seq <- as.integer(x=seq_len(length.out=nchar(x=seq_l)) *
                                error)
    }else{
      flank_seq <- rep(x=0,
                       times=length(x=seq_len(length.out=nchar(x=seq_l))))
    }
    if (min_match_flank >= 1L){
      if (nchar(x=seq_l) > min_match_flank){
        flank_seq[seq_len(length.out=min_match_flank)] <- -1
      }else{
        return(FQ)
      }
    }
    if (!anchored){
      maxlen <-  max(width(x=FQ)) - nchar(x=seq_l)
      if (maxlen > 0){
        seq_l <- paste0(paste(rep(x="N", times=maxlen),
                              collapse=""),
                        seq_l)
      }
      flank_seq <- c(flank_seq, rep(x=0, times=maxlen))
    }
    FQ <- trimLRPatterns(Lpattern=seq_l,
                         subject=FQ,
                         max.Lmismatch=flank_seq,
                         with.Lindels=indels)
  }
  ### Returning the final object.
  return(FQ)
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
#' trimL <- cutRigthAdapter(fastqDir="Files_FASTQ",
#'                          fastq="example_seq.read1.fastq",
#'                          adapters="adapters.txt",
#'                          error=0.2,
#'                          min_match_flank=3L,
#'                          anchored=TRUE,
#'                          indels=FALSE,
#'                          workDir="D:/Vasily Grinev")
#' @export
#' @importFrom ShortRead readFastq trimLRPatterns
#' @importFrom Biostrings width
#' @importFrom methods is
#' @importFrom utils read.table
cutRigthAdapter <- function(fastqDir=NULL,
                            fastq,
                            adapters,
                            error=0.2,
                            min_match_flank=3L,
                            anchored=TRUE,
                            indels=FALSE,
                            workDir=NULL){
  ### Setting the working directory.
  if (is.null(x=workDir)){
    workDir <- getwd()
  }
  ### Retrieving of the FASTQ data.
  if (is.object(x=fastq)){
    FQ=fastq
  }else{
    ### Full path to the FASTQ file.
    if (is.null(x=fastqDir)){
      path <- paste(workDir, fastqDir, sep="")
    }else{
      path <- paste(workDir, fastqDir, sep="/")
    }
    ### Retrieving and validation of the FASTQ file extension.
    frt <- tools::file_ext(x=fastq)
    if (!frt %in% c("fastq", "gz")){
      stop("Invalid format of FASTQ file")
    }
    FQ <- readFastq(dirPath=path, pattern=fastq, withIds=TRUE)
  }
  ### Loading of the adapter sequence(-s).
  ADs <- read.table(file=paste(workDir, adapters, sep="/"),
                    sep="\t",
                    header=TRUE,
                    quote="\"",
                    as.is=TRUE)
  ### Rigth-end trimming of adapter sequence(-s).
  for (i in 1:nrow(x=ADs)){
    seq_l <- ADs[i, 3]
    if (error > 0){
      flank_seq <- as.integer(x=seq_len(length.out=nchar(x=seq_l)) *
                                error)
    }else{
      flank_seq <- rep(x=0,
                       times=length(x=seq_len(length.out=nchar(x=seq_l))))
    }
    if (min_match_flank >= 1L){
      if (nchar(x=seq_l) > min_match_flank){
        flank_seq[seq_len(length.out=min_match_flank)] <- -1
      }else{
        return(FQ)
      }
    }
    if (!anchored){
      maxlen <-  max(width(x=FQ)) - nchar(x=seq_l)
      if (maxlen > 0){
        seq_l <- paste0(seq_l, paste(rep(x="N", times=maxlen),
                                     collapse=""))
      }
      flank_seq <- c(flank_seq, rep(x=0, times=maxlen))
    }
    FQ <- trimLRPatterns(Rpattern=seq_l,
                         subject=FQ,
                         max.Rmismatch=flank_seq,
                         with.Rindels=indels)
  }
  ### Returning the final object.
  return(FQ)
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
#' ad <- findAdapters(fastqDir="Files_FASTQ",
#'                    fastq="2-Galkina-A_S2_L001_R1_001.fastq.gz",
#'                    n=1e5,
#'                    adapters="adapters.txt",
#'                    workDir="D:/Vasily Grinev")
#' @export
#' @importFrom Biostrings vmatchPattern width
#' @importFrom GenomicRanges GRanges intersect
#' @importFrom methods is
#' @importFrom utils read.table
#' @importFrom tools file_ext
findAdapters <- function(fastqDir=NULL,
                         fastq,
                         n=NULL,
                         adapters,
                         workDir=NULL){
  ### Setting the working directory.
  if (is.null(x=workDir)){
    workDir <- getwd()
  }
  ### Full path to the FASTQ file(-s).
  if (is.null(x=fastqDir)){
    path <- paste(workDir, fastqDir, sep="")
  }else{
    path <- paste(workDir, fastqDir, sep="/")
  }
  ### Retrieving and validation of the FASTQ file(-s) extension.
  frt <- tools::file_ext(x=fastq)
  if (!frt %in% c("fastq", "gz")){
    stop("Invalid file format")
  }
  ### Loading the FASTQ data.
  reads <- readFASTQ(fastqDir=fastqDir,
                     fastq=fastq,
                     n=n,
                     workDir=workDir)[[1]]
  ### Loading the adapter sequences.
  ADs <- read.table(file=paste(workDir, adapters, sep="/"),
                    sep="\t",
                    header=TRUE,
                    quote="\"",
                    as.is=TRUE)
  hits <- list()
  for (i in 1:nrow(x=ADs)){
    x <- unlist(x=vmatchPattern(pattern=ADs$adapter_sequence[i],
                                max.mismatch=2, subject=reads))
    x <- GRanges(seqnames=names(x=x), ranges=x)
    x <- data.frame(intersect(x=range(x), y=x))[, 2:3]
    x <- unlist(x=apply(X=x,
                        MARGIN=1,
                        FUN=function(y){seq(from=y[1], to=y[2])}))
    x <- cut(x=x,
             breaks=seq(from=1,
                        to=max(width(x=reads)) + 5,
                        by=5))
    x <- data.frame(table(x))
    hits[[i]] <- x$Freq/(length(x=reads) * 6)
  }
  if (max(width(x=reads)) %% 5 == 0){
    positions <- paste(seq(from=1, to=max(width(x=reads)), by=5),
                       seq(from=5, to=max(width(x=reads)), by=5),
                       sep="-")
  }else{
    positions <- paste(seq(from=1, to=max(width(x=reads)), by=5),
                       c(seq(from=5, to=max(width(x=reads)), by=5),
                         max(width(x=reads))),
                       sep="-")
  }
  hits <- data.frame(cbind(positions, do.call(what=cbind, args=hits)))
  colnames(x=hits) <- c("positions", ADs$adapter_name)
  hits[, -1] <- apply(X=hits[, -1], MARGIN=2, FUN=as.numeric)
  ### Returning the final object.
  return(hits)
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
#' conts <- findForeignSeqs(fastqDir="Files_FASTQ",
#'                          fastq="2-Galkina-A_S2_L001_R1_001.fastq.gz",
#'                          n=1e5,
#'                          contaminants="contaminants.txt",
#'                          workDir="D:/Vasily Grinev")
#' @export
#' @importFrom Biostrings vmatchPattern width
#' @importFrom GenomicRanges GRanges intersect
#' @importFrom utils read.table
#' @importFrom methods is
#' @importFrom tools file_ext
findForeignSeqs <- function(fastqDir=NULL,
                            fastq,
                            n=NULL,
                            contaminants,
                            workDir=NULL){
  ### Setting the working directory.
  if (is.null(x=workDir)){
    workDir <- getwd()
  }
  ### Full path to the FASTQ file(-s).
  if (is.null(x=fastqDir)){
    path <- paste(workDir, fastqDir, sep="")
  }else{
    path <- paste(workDir, fastqDir, sep="/")
  }
  ### Retrieving and validation of the FASTQ file(-s) extension.
  frt <- tools::file_ext(x=fastq)
  if (!frt %in% c("fastq", "gz")){
    stop("Invalid file format")
  }
  ### Loading the FASTQ data.
  reads <- readFASTQ(fastqDir=fastqDir,
                     fastq=fastq,
                     n=n,
                     workDir=workDir)[[1]]
  ### Loading the foreign sequences.
  FSs <- read.table(file=paste(workDir, contaminants, sep="/"),
                    sep="\t",
                    header=TRUE,
                    quote="\"",
                    as.is=TRUE)
  hits <- list()
  for (i in 1:nrow(x=FSs)){
    x <- unlist(x=vmatchPattern(pattern=FSs$foreign_sequence[i],
                                max.mismatch=2,subject=reads))
    x <- GRanges(seqnames=names(x=x), ranges=x)
    x <- data.frame(intersect(x=range(x), y=x))[, 2:3]
    x <- unlist(x=apply(X=x,
                        MARGIN=1,
                        FUN=function(y){seq(from=y[1], to=y[2])}))
    x <- cut(x=x,
             breaks=seq(from=1,
                        to=max(width(x=reads)) + 5,
                        by=5))
    x <- data.frame(table(x))
    hits[[i]] <- x$Freq/(length(x=reads) * 6)
  }
  if (max(width(x=reads)) %% 5 == 0){
    positions <- paste(seq(from=1, to=max(width(x=reads)), by=5),
                       seq(from=5, to=max(width(x=reads)), by=5),
                       sep="-")
  }else{
    positions <- paste(seq(from=1, to=max(width(x=reads)), by=5),
                       c(seq(from=5, to=max(width(x=reads)), by=5),
                         max(width(x=reads))),
                       sep="-")
  }
  hits <- data.frame(cbind(positions, do.call(what=cbind, args=hits)))
  colnames(x=hits) <- c("positions", FSs$foreign_name)
  hits[, -1] <- apply(X=hits[, -1], MARGIN=2, FUN=as.numeric)
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
#' report <- reportQAResults(x=qa,
#'                           output=NULL,
#'                           workDir="D:/Vasily Grinev")
#' @export
#' @importFrom Biostrings DNAStringSet writeXStringSet
#' @importFrom grid grid.draw textGrob gpar gTree gList
#' @importFrom gridExtra tableGrob ttheme_minimal
#' @importFrom vioplot vioplot
#' @importFrom graphics plot axis title rect legend mtext text barplot
#' @importFrom stats ks.test density
#' @importFrom grDevices pdf dev.off
reportQAResults <- function(x,
                            output=NULL,
                            workDir=NULL){
  ### Textual reports.
  ##  Template of textual table.
  tab <- data.frame(Parameter=c("Sample name",
                                "File name",
                                "Total number of reads",
                                "Total number of bases",
                                "Length of reads, min",
                                "Length of reads, 1st quartile",
                                "Length of reads, median",
                                "Length of reads, 3rd quartile",
                                "Length of reads, max",
                                "Per entire reads quality, min",
                                "Per entire reads quality, 1st quartile",
                                "Per entire reads quality, median",
                                "Per entire reads quality, 3rd quartile",
                                "Per entire reads quality, max",
                                "Frequency of A, min",
                                "Frequency of A, 1st quartile",
                                "Frequency of A, median",
                                "Frequency of A, 3rd quartile",
                                "Frequency of A, max",
                                "Frequency of C, min",
                                "Frequency of C, 1st quartile",
                                "Frequency of C, median",
                                "Frequency of C, 3rd quartile",
                                "Frequency of C, max",
                                "Frequency of G, min",
                                "Frequency of G, 1st quartile",
                                "Frequency of G, median",
                                "Frequency of G, 3rd quartile",
                                "Frequency of G, max",
                                "Frequency of T, min",
                                "Frequency of T, 1st quartile",
                                "Frequency of T, median",
                                "Frequency of T, 3rd quartile",
                                "Frequency of T, max",
                                "Frequency of N, min",
                                "Frequency of N, 1st quartile",
                                "Frequency of N, median",
                                "Frequency of N, 3rd quartile",
                                "Frequency of N, max",
                                "Frequency of GC, min",
                                "Frequency of GC, 1st quartile",
                                "Frequency of GC, median",
                                "Frequency of GC, 3rd quartile",
                                "Frequency of GC, max",
                                "Low-complexity of reads, min",
                                "Low-complexity of reads, 1st quartile",
                                "Low-complexity of reads, median",
                                "Low-complexity of reads, 3rd quartile",
                                "Low-complexity of reads, max"),
                    Value=0)
  ##  Values of textual table.
  tab[1, 2] <- x$SN
  tab[2, 2] <- x$FN
  tab[3, 2] <- x$TNRs
  tab[4, 2] <- x$TNBs
  tab[5:9, 2] <- x$RLs[-4]
  tab[10:14, 2] <- sprintf(x=x$PERQs$metrics, fmt="%.1f")
  tab[15:39, 2] <- formatC(x=c(t(x$BCs$metrics[, -1])), format="e", digits=2)
  tab[40:44, 2] <- formatC(x=x$GCs$metrics, format="e", digits=2)
  tab[45:49, 2] <- x$DSs$metrics
  ##  Saving of textual table as *.PDF file.
  if (is.null(x=output)){
    file_output <- paste(sprintf("QA metrics of %s", x$SN),
                         "pdf",
                         sep=".")
  }else{
    file_output <- paste(sprintf("QA metrics of %s", output),
                         "pdf",
                         sep=".")
  }
  pdf(file=paste(workDir, file_output, sep="/"), height=nrow(x=tab)/3)
  theme=ttheme_minimal(core=list(fg_params=list(hjust=0,
                                                x=0.03,
                                                fontsize=10),
                                 bg_params=list(fill=rep(c("azure1",
                                                           "azure3"),
                                                         each = 1),
                                                col="black")),
                       colhead=list(fg_params=list(col="white"),
                                    bg_params=list(fill="darkcyan",
                                                   col="black")))
  table <- tableGrob(d=tab, rows=NULL, theme=theme)
  title1 <- textGrob(label="Quality of the raw reads. Summary statistics",
                     just="centre",
                     vjust=-60,
                     gp=gpar(col="black",
                             font=2,
                             fontsize=12))
  title2 <- textGrob(label=sprintf("Sample: %s", x$SN),
                     just="centre",
                     vjust=-63.5,
                     gp=gpar(col="black",
                             fontsize=11))
  title3 <- textGrob(label=sprintf(paste("Date & Time:",
                                         format(x=Sys.time(), "%b %d %Y %X"),
                                         sep=" ")),
                     just="centre",
                     vjust=-67.8,
                     gp=gpar(col="black",
                             fontsize=10))
  note <- "Generated by assessQRawReads() & reportQAResults()"
  note <- textGrob(label=sprintf(note),
                   just="centre",
                   vjust=88,
                   gp=gpar(col="black",
                           fontsize=8))
  gt <- gTree(children=gList(table, title1, title2, title3, note))
  grid.draw(gt)
  suppressMessages(expr=dev.off())
  ##  Saving of overrepresented reads.
  ORs <- DNAStringSet(x=x$ORs$sequence)
  names(x=ORs) <- paste(paste("seq",
                              formatC(x=1:length(x=x$ORs$sequence),
                                      width=4, flag="0"), sep=""),
                        x$ORs$count, sep=", ")
  if (is.null(x=output)){
    file_output <- paste(sprintf("Overrepresented reads of %s", x$SN),
                         "fasta",
                         sep=".")
  }else{
    file_output <- paste(sprintf("Overrepresented reads of %s", output),
                         "fasta",
                         sep=".")
  }
  writeXStringSet(x=ORs,
                  filepath=paste(workDir, file_output, sep="/"),
                  format="fasta")
  ### Graphical reports.
  ##  Per entire reads quality plots.
  if (is.null(x=output)){
    file_output <- paste(sprintf("QA plots of %s", x$SN),
                         "pdf",
                         sep=".")
  }else{
    file_output <- paste(sprintf("QA plots of %s", output),
                         "pdf",
                         sep=".")
  }
  pdf(file=paste(workDir, file_output, sep="/"),
      onefile=TRUE,
      paper="a4",
      width=8.27, height=11.69)
  par(mfrow=c(3, 2),
      omi=c(0.79, 0.79, 1, 0.79))
  plot(0, type="n", axes=FALSE, ann=FALSE, bty="n",
       xlim=c(0, 2), ylim=c(0, 40))
  vioplot(x=x$PERQs$means[1:10000],
          add=TRUE,
          col=rgb(red=86, green=180, blue=233, maxColorValue=255),
          rectCol=rgb(red=0, green=158, blue=115, maxColorValue=255),
          pchMed=0, colMed=NA,
          lwd=0.9,
          wex=1,
          frame.plot=FALSE)
  segments(x0=0.9, y0=median(x=x$PERQs$means),
           x1=1.1, y1=median(x$PERQs$means),
           col=rgb(red=240, green=228, blue=66, maxColorValue=255),
           lwd=3, lend=2)
  axis(side=1, at=c(0, 1, 2), labels=FALSE, lwd=0.8)
  axis(side=2, at=c(0, 10, 20, 30, 40), labels=c(0, 10, 20, 30, 40),
       las=1, lwd=0.8)
  title(main="Distribution of reads quality", cex.main=1.5)
  title(xlab="Standard vioplot", line=1, cex.lab=1.2)
  title(ylab="Average read quality scores", line=2.3, cex.lab=1.2)
  per <- data.frame(threshold=0:ceiling(x=max(x=x$PERQs$means)), percent=0)
  per[1, 2] <- 100
  L <- length(x=x$PERQs$means)
  for (i in 1:ceiling(x=max(x=x$PERQs$means))){
    p <- length(x=x$PERQs$means[x$PERQs$means > i])/L * 100
    per[i + 1, 2] <- round(x=p, digits=1)
  }
  plot(per, type="l", bty="n",
       xlim=c(0, 40), ylim=c(0, 100),
       lwd=1,
       col=rgb(red=0, green=114, blue=178, maxColorValue=255),
       ann=FALSE, axes=FALSE)
  axis(side=1,
       at=c(0, 10, 20, 30, 40),
       labels=c(0, 10, 20, 30, 40),
       las=1, lwd=0.8)
  axis(side=2,
       at=c(0, 20, 40, 60, 80, 100),
       labels=c(0, 20, 40, 60, 80, 100),
       las=1, lwd=0.8)
  title(main="Thresholding of reads quality", cex.main=1.5)
  title(xlab="Quality threshold", line=2.5, cex.lab=1.2)
  title(ylab="Reads exceeding quality level, %", line=2.8, cex.lab=1.2)
  rect(xleft=0, ybottom=0, xright=20, ytop=100, border=NA,
       col=rgb(red=240, green=228, blue=66, maxColorValue=255, alpha=40))
  rect(xleft=20, ybottom=0, xright=30, ytop=100, border=NA,
       col=rgb(red=230, green=159, blue=0, maxColorValue=255, alpha=40))
  rect(xleft=30, ybottom=0, xright=40, ytop=100, border=NA,
       col=rgb(red=213, green=94, blue=0, maxColorValue=255, alpha=70))
  text(x=20, y=20, srt=90,
       labels=paste(paste("Q20",
                          round(x=length(x=x$PERQs$means[x$PERQs$means >
                                                           20])/L * 100, digits=1),
                          sep=", "),
                    "%",
                    sep=""))
  text(x=30, y=20, srt=90,
       labels=paste(paste("Q30",
                          round(x=length(x=x$PERQs$means[x$PERQs$means >
                                                           30])/L * 100, digits=1),
                          sep=", "),
                    "%",
                    sep=""))
  ##  Per cycle reads quality plots.
  bxp(z=list(stats=t(x$PCQRs[, -1]),
             n=rep(x=length(x=x$PERQs$means), times=nrow(x=x$PCQRs))),
      ylim=c(0, 40),
      boxfill=rgb(red=240, green=228, blue=66, maxColorValue=255),
      frame=FALSE,
      axes=FALSE,
      lwd=0.6,
      medlwd=1, medcol=rgb(red=204, green=121, blue=167, maxColorValue=255),
      whisklty="solid")
  axis(side=1,
       at=1:length(x=x$PCQRs[, 1]),
       labels=FALSE,
       las=2, lwd=0.8)
  axis(side=2,
       at=c(0, 10, 20, 30, 40),
       labels=c(0, 10, 20, 30, 40),
       las=1, lwd=0.8)
  text(x=seq(from=1, to=length(x=x$PCQRs[, 1]), by=2),
       par("usr")[3] - 3, offset=-0.5, cex=0.9,
       labels=x$PCQRs[, 1][seq(from=1, to=length(x=x$PCQRs[, 1]), by=2)],
       srt=55, pos=2, xpd=TRUE)
  title(main="Per cycle quality of reads", cex.main=1.5)
  title(xlab="Cycles", line=3.7, cex.lab=1.2)
  title(ylab="Quality scores across of all bases", line=2.3, cex.lab=1.2)
  rect(xleft=0.5, ybottom=0, xright=length(x=x$PCQRs[, 1]) + 0.5, ytop=20,
       border=NA,
       col=rgb(red=240, green=228, blue=66, maxColorValue=255, alpha=40))
  rect(xleft=0.5, ybottom=20, xright=length(x=x$PCQRs[, 1]) + 0.5, ytop=30,
       border=NA,
       col=rgb(red=230, green=159, blue=0, maxColorValue=255, alpha=40))
  rect(xleft=0.5, ybottom=30, xright=length(x=x$PCQRs[, 1]) + 0.5, ytop=40,
       border=NA,
       col=rgb(red=213, green=94, blue=0, maxColorValue=255, alpha=70))
  ##  Base composition of the entire set of reads.
  colrs <- c(rgb(red=0, green=114, blue=178, maxColorValue=255),
             rgb(red=213, green=94, blue=0, maxColorValue=255),
             rgb(red=86, green=180, blue=233, maxColorValue=255),
             rgb(red=230, green=159, blue=0, maxColorValue=255),
             rgb(red=204, green=121, blue=167, maxColorValue=255))
  bxp(z=list(stats=t(x$BCs$metrics[, -1]),
             n=rep(x=length(x=x$PERQs$means), times=nrow(x=x$BCs$metrics))),
      ylim=c(0, 0.5),
      boxfill=colrs,
      frame=FALSE,
      axes=FALSE,
      lwd=0.6,
      medlwd=1,
      boxwex=0.6,
      whisklty="solid")
  axis(side=1,
       at=c(1, 2, 3, 4, 5),
       labels=x$BCs$metrics[, 1],
       las=1, lwd=0.8)
  axis(side=2,
       at=c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5),
       labels=c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5),
       las=2, lwd=0.8)
  title(main="Base composition of reads", cex.main=1.5)
  title(xlab="Base", line=2.5, cex.lab=1.2)
  title(ylab="Frequency", line=2.8, cex.lab=1.2)
  ##  Per cycle base composition of reads.
  matplot(x$PCCBs[, -1],
          type="l", lty=1, lwd=1, col=colrs,
          ylim=c(0, 1),
          ann=FALSE, axes=FALSE)
  axis(side=1,
       at=1:length(x=x$PCCBs[, 1]),
       labels=FALSE,
       las=2, lwd=0.8)
  axis(side=2,
       at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
       labels=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
       las=1, lwd=0.8)
  text(x=seq(from=1, to=length(x=x$PCCBs[, 1]), by=2),
       par("usr")[3] - 0.08, offset=-0.5, cex=0.9,
       labels=x$PCCBs[, 1][seq(from=1, to=length(x=x$PCCBs[, 1]), by=2)],
       srt=55, pos=2, xpd=TRUE)
  title(main="Per cycle composition of reads", cex.main=1.5)
  title(xlab="Cycles", line=3.8, cex.lab=1.2)
  title(ylab="Base frequency", line=2.3, cex.lab=1.2)
  legend(x="topleft",
         legend=colnames(x=x$PCCBs)[-1],
         col=colrs, lty=1, lwd=1, bty="n")
  ##  GC composition of the entire set of reads.
  plot(x$GCs$density,
       bty="n",
       type="l", lwd=1,
       xlim=c(0, 1),
       ann=FALSE, axes=FALSE,
       col=rgb(red=0, green=114, blue=178, maxColorValue=255))
  lines(density(x=rnorm(n=length(x=x$GCs$frequency),
                        mean=mean(x=x$GCs$frequency),
                        sd=sd(x=x$GCs$frequency)))[1:2],
        xlim=c(0, 1),
        lwd=1,
        col=rgb(red=213, green=94, blue=0, maxColorValue=255))
  axis(side=1,
       at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
       labels=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
       lwd=0.8)
  axis(side=2,
       at=c(0.0, max(x=x$GCs$density[, 2])/4 * 1:4),
       labels=round(x=c(0.0, max(x=x$GCs$density[, 2])/4 * 1:4), digits=2),
       las=1, lwd=0.8)
  title(main="GC composition of reads", cex.main=1.5)
  title(xlab="GC frequency", line=2.5, cex.lab=1.2)
  title(ylab="Density", line=2.8, cex.lab=1.2)
  legend(x="topleft",
         legend=c("theoretical", "empirical"),
         col=c(rgb(red=213, green=94, blue=0, maxColorValue=255),
               rgb(red=0, green=114, blue=178, maxColorValue=255)),
         lty=1, lwd=1, bty="n")
  suppressWarnings(expr=legend(x="topright",
                               legend=paste("KS test: p = ",
                                            ks.test(x=x$GCs$frequency, y="pnorm")$p.value,
                                            sep=""),
                               lty=0, bty="n"))
  ##  Title, subtitle, sub-subtitle & footnote.
  mtext(text="Quality of the raw reads. Summary plots",
        side=3,
        line=3.5,
        cex=1.1,
        col="black",
        font=2,
        outer=TRUE)
  mtext(text=sprintf("Sample: %s", x$SN),
        side=3,
        line=2,
        cex=0.95,
        col="black",
        outer=TRUE)
  mtext(text=sprintf(paste("Date & Time:",
                           format(x=Sys.time(), "%b %d %Y %X"),
                           sep=" ")),
        side=3,
        line=0.5,
        cex=0.8,
        col="black",
        outer=TRUE)
  note <- "Generated by assessQRawReads() & reportQAResults()"
  mtext(text=sprintf(note),
        side=1,
        line=2,
        cex=0.7,
        col="black",
        outer=TRUE)
  ##  Summarization the low-complexity of reads (using dusty scores).
  plot(x$DSs$density,
       bty="n",
       type="l", lwd=1,
       xlim=c(2, 5),
       ann=FALSE, axes=FALSE,
       col=rgb(red=0, green=114, blue=178, maxColorValue=255))
  axis(side=1,
       at=c(2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0),
       labels=c(2.0, 2.5, 3.0, 3.5,  4.0, 4.5, 5.0),
       lwd=0.8)
  axis(side=2,
       at=c(0.0, max(x=x$DSs$density[, 2])/4 * 1:4),
       labels=round(x=c(0.0, max(x=x$DSs$density[, 2])/4 * 1:4), digits=2),
       las=1, lwd=0.8)
  title(main="Low-complexity of reads", cex.main=1.5)
  title(xlab=expression("Complexity, log"[10] * " (dusty scores)"),
        line=2.75, cex.lab=1.2)
  title(ylab="Density", line=3, cex.lab=1.2)
  legend(x=3.3, y=max(x=x$DSs$density[, 2]),
         legend=c("minimum:", "Q1:", "median:", "Q3:", "maximum:"),
         lty=0, bty="n")
  legend(x=4.1, y=max(x=x$DSs$density[, 2]),
         legend=round(x=log10(x=x$DSs$metrics), digits=2),
         lty=0, bty="n")
  rect(xleft=2.0, ybottom=0, xright=3.0, ytop=max(x=x$DSs$density[, 2]),
       density=20, border=NA, lwd=0.3,
       col=rgb(red=0, green=158, blue=115, maxColorValue=255))
  ##  Frequency of reads.
  barplot(formula=x$FRs[, 3] ~ x$FRs[, 1], data=x$FRs,
          axes=FALSE, ann=FALSE, ylim=c(0, 1), border=NA,
          col=colorRampPalette(c(rgb(red=213,
                                     green=94,
                                     blue=0,
                                     maxColorValue=255),
                                 rgb(red=230,
                                     green=159,
                                     blue=0,
                                     maxColorValue=255)))(nrow(x=x$FRs)))
  axis(side=2,
       at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
       labels=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
       las=1, lwd=0.8)
  title(main="Frequency of reads", cex.main=1.5)
  title(xlab="Occurrence of reads",
        line=2.5, cex.lab=1.2)
  title(ylab="Frequency", line=3, cex.lab=1.2)
  ##  Contamination of reads with adapters.
  if (!is.null(x=nrow(x=x$ASs))){
    colrs <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288",
               "#AA4499", "#44AA99", "#999933", "#882255", "#661100",
               "#6699CC", "#888888")
    matplot(x$ASs[, -1],
            ann=FALSE,
            axes=FALSE,
            ylim=c(0, 0.35),
            type="l", lty=1, lwd=1,
            col=colrs[1:ncol(x=x$ASs[, -1])])
    axis(side=1,
         at=1:length(x=x$ASs[, 1]),
         labels=FALSE,
         lwd=0.8)
    axis(side=2,
         at=c(0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35),
         labels=c(0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35), las=1,
         lwd=0.8)
    text(x=seq(from=1, to=length(x=x$ASs[, 1]), by=2),
         par("usr")[3] - 0.028, offset=-0.5,
         labels=x$ASs[, 1][seq(from=1, to=length(x=x$ASs[, 1]), by=2)],
         cex=0.9, srt=55, pos=2, xpd=TRUE)
    title(main="Adapter contaminations", cex.main=1.5)
    title(xlab="Cycles", line=3.8, cex.lab=1.2)
    title(ylab="Adapter frequency", line=3, cex.lab=1.2)
    legend(x="topleft",
           legend=colnames(x=x$ASs)[-1],
           cex=0.8,
           lty=1, lwd=1,
           col=colrs[1:ncol(x=x$ASs[, -1])],
           bty="n")
  }
  ##  Contamination of reads with foreign sequences.
  if (!is.null(x=nrow(x=x$FSs))){
    colrs <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288",
               "#AA4499", "#44AA99", "#999933", "#882255", "#661100",
               "#6699CC", "#888888")
    matplot(x$FSs[, -1],
            ann=FALSE,
            axes=FALSE,
            ylim=c(0, 0.35),
            type="l", lty=1, lwd=1,
            col=colrs[1:ncol(x=x$FSs[, -1])])
    axis(side=1,
         at=1:length(x=x$FSs[, 1]),
         labels=FALSE,
         lwd=0.8)
    axis(side=2,
         at=c(0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35),
         labels=c(0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35), las=1,
         lwd=0.8)
    text(x=seq(from=1, to=length(x=x$FSs[, 1]), by=2),
         par("usr")[3] - 0.028, offset=-0.5,
         labels=x$FSs[, 1][seq(from=1, to=length(x=x$FSs[, 1]), by=2)],
         cex=0.9, srt=55, pos=2, xpd=TRUE)
    title(main="Foreign sequence contaminations", cex.main=1.5)
    title(xlab="Cycles", line=3.8, cex.lab=1.2)
    title(ylab="Sequence frequency", line=3, cex.lab=1.2)
    legend(x="topleft",
           legend=colnames(x=x$FSs)[-1],
           cex=0.8,
           lty=1, lwd=1,
           col=colrs[1:ncol(x=x$FSs[, -1])],
           bty="n")
  }
  ##  Title, subtitle, sub-subtitle & footnote.
  mtext(text="Quality of the raw reads. Summary plots (continuation)",
        side=3,
        line=3.5,
        cex=1.1,
        col="black",
        font=2,
        outer=TRUE)
  mtext(text=sprintf("Sample: %s", x$SN),
        side=3,
        line=2,
        cex=0.95,
        col="black",
        outer=TRUE)
  mtext(text=sprintf(paste("Date & Time:",
                           format(x=Sys.time(), "%b %d %Y %X"),
                           sep=" ")),
        side=3,
        line=0.5,
        cex=0.8,
        col="black",
        outer=TRUE)
  note <- "Generated by assessQRawReads() & reportQAResults()"
  mtext(text=sprintf(note),
        side=1,
        line=2,
        cex=0.7,
        col="black",
        outer=TRUE)
  suppressMessages(expr=dev.off())
}


