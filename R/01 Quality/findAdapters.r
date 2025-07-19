#' Identify adapter sequences in short reads
#' @description Identification of the adapter sequences in short reads.
#' @param fastqDir a character string specifying the name of the directory
#'     containing the FASTQ file(s). NULL by default that mean the current
#'     working directory.
#' @param fastq a character string specifying the name of the FASTQ file
#'     containing short reads. Allowed formats are "fastq" or "fastq.gz".
#' @param n integer value giving the number of reads loaded into the computer's
#'     RAM. The default value is NULL, at which all reads are loaded into RAM.
#' @param adapters a character string specifying the name of the tab-delimited
#'     TXT file containing adapter sequences. The default value is NULL, at
#'     which the search for adaptor sequences in short reads is not carried out.
#'     If so, this file must contain the following three fields:
#'     i) adapter_id         - adapter ID;
#'     ii) adapter_name      - adapter name;
#'     iii) adapter_sequence - adapter sequence.
#' @param workDir character string specifying the path to and name of working
#'     directory. NULL by default that mean the current working directory.
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
#' Last updated: July 13, 2025.

findAdapters <- function(fastqDir=NULL,
                         fastq,
                         n=NULL,
                         adapters,
                         workDir=NULL){
    ### Loading the required packages.
    #   This code was successfully tested with the packages Biostrings v.2.72.1
    #   and GenomicRanges v.1.56.1.
    suppressMessages(expr=library(package=Biostrings))
    suppressMessages(expr=library(package=GenomicRanges))
    ### Setting the working directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Loading the required auxiliary function.
    source(file=paste(workDir, "readFASTQ.r", sep="/"))
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

