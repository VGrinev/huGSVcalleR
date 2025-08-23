#' Filter a BAM file
#' @description Multiparameter BAM file filtering.
#' @param bamDir character string giving the name (or path to and name) of
#'     directory with BAM file(s). NULL by default, which means the current
#'     working directory will be used instead of the specified BAM directory.
#' @param bamFile character string giving the name of input BAM file.
#' @param gr character string specifying the name of the tab-delimited TXT file
#'     containing coordinates of genomic location(-s) of interest. The default
#'     value is NULL. If so, this file must contains the following four fields:
#'     i) seqnames - the name of chromosome or scaffold with prefix "chr";
#'     ii) start   - start coordinate of genomic interval of interest;
#'     iii) end    - end coordinate of genomic interval of interest;
#'     iv) strand  - (optionally) strand information about genomic interval
#'                   of interest.
#'     TXT file with genomic location(-s) must be in working directory.
#' @param flag logical vector used to filter reads based on their flag entry.
#'     This is most easily created with the helper function scanBamFlag() from
#'     package Rsamtools. The function scanBamFlag() recognizes the following
#'     flags: isPaired, isProperPair, isUnmappedQuery, hasUnmappedMate,
#'     isMinusStrand, isMateMinusStrand, isFirstMateRead, isSecondMateRead,
#'     isNotPrimaryRead, isSecondaryAlignment, isNotPassingQualityControls,
#'     isDuplicate and isSupplementaryAlignment. The chosen flag must have
#'     value TRUE to be involved in filtration process and only such flag
#'     should be submitted to the function. Filtering by flag is not performed
#'     by default.
#' @param tlen integer value limiting the maximum template length. NULL by
#'     default that mean no any limitations.
#' @param mapq integer value limiting the minimum mapping quality. NULL by
#'     default that mean no any limitations.
#' @param index logical, FALSE by default. If TRUE, it allows to sort and index
#'     a new filtered BAM file.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @return filtered BAM file (optionally, sorted and indexed if index=TRUE)
#'     with postfix "filtered".
#' @author Vasily V. Grinev
#' @examples
#' filterBamFile(bamDir="Files_BAM",
#'               bamFile="SRR21721515.bam",
#'               gr="genomic_intervals.txt",
#'               flag=scanBamFlag(isPaired=TRUE, isProperPair=TRUE),
#'               tlen=300,
#'               mapq=20,
#'               index=TRUE,
#'               workDir="D:/Vasily Grinev")
#' @export
#  Last updated: July 23, 2025.

filterBamFile <- function(bamDir=NULL,
                          bamFile,
                          gr=NULL,
                          flag=scanBamFlag(),
                          tlen=NULL,
                          mapq=NULL,
                          index=FALSE,
                          workDir=NULL){
    ### Loading the required packages.
    #   This code was successfully tested with the packages parallel and
    #   Rsamtools v.2.20.0.
    suppressMessages(expr=library(package=parallel))
    suppressMessages(expr=library(package=Rsamtools))
    ### Setting the working directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Full path to the BAM file.
    if (is.null(x=bamDir)){
        path <- workDir
    }else{
        path <- paste(workDir, bamDir, sep="/")
    }
    ### Development an object of class IRangesList to filter the BAM file of
    #   interest in given genomic interval(-s) (if so).
    if (!is.null(x=gr)){
        GRs <- read.table(file=paste(workDir, gr, sep="/"),
                          sep="\t",
                          header=TRUE,
                          quote="\"",
                          as.is=TRUE)
        GRs$split <- GRs$seqnames
        GRs <- makeGRangesListFromDataFrame(df=GRs, split.field="split")
        GRs <- IRangesList(GRs)
    }
    ### Filtering the BAM file of interest.
    param <- ScanBamParam(flag=flag)
    if (!is.null(x=gr)){
        param@which <- GRs
    }
    if (!is.null(x=tlen)){
        param@what <- "isize"
    }
    if (!is.null(x=mapq)){
        param@mapqFilter <- as.integer(x=mapq)
    }
    if (!is.null(x=tlen)){
        filter <- FilterRules(list(function(x){x$isize <= tlen}))
    }else{
        filter <- FilterRules()
    }
    dest <- gsub(pattern="bam",
                 replacement="filtered.bam",
                 x=paste(path, bamFile, sep="/"))
    f <- filterBam(file=paste(path, bamFile, sep="/"),
                   index=paste(path, bamFile, sep="/"),
                   filter=filter,
                   param=param,
                   destination=dest,
                   indexDestination=FALSE)
    if (isTRUE(x=index)){
        s <- sortBam(file=f,
                     destination=gsub(pattern=".bam", replacement="", x=f),
                     nThreads=detectCores(logical=TRUE) - 1)
        indexBam(files=s)
    }
}
