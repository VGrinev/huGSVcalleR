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
#  Last updated: July 25, 2025.

cutLeftAdapter <- function(fastqDir=NULL,
                           fastq,
                           adapters,
                           error=0.2,
                           min_match_flank=3L,
                           anchored=TRUE,
                           indels=FALSE, 
                           workDir=NULL){
    ### Loading the required package.
    #   This code was successfully tested with the package ShortRead v.1.67.0.
    suppressMessages(expr=library(package=ShortRead))
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
