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
#' Last updated: July 25, 2025.

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
    ### Loading the required package.
    #   This code was successfully tested with the package ShortRead v.1.67.0.
    suppressMessages(expr=library(package=ShortRead))
    ### Loading the required auxiliary functions.
    source(file=paste(workDir, "cutLeftAdapter.r", sep="/"))
    source(file=paste(workDir, "cutRigthAdapter.r", sep="/"))
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
