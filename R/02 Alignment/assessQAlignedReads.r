#' Evaluate alignment quality metrics
#' @description Evaluation of the alignment quality metrics for aligned DAN-Seq
#'     and/or RNA-Seq short reads.
#' @param bamDir character vector giving the path to and name of directory
#'     with BAM file(s).
#' @param bamFile character vector giving the name of BAM file to be evaluated.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return object of class list containing all main quality metrics for aligned
#'     short reads.
#' @author Vasily V. Grinev
#' @examples
#' qal <- assessQAlignedReads(bamDir="Files_BAM",
#'                            bamFile="SRR21721515.bam",
#'                            workDir="D:/Vasily Grinev")
#' @export
#' Last updated: July 18, 2025.

assessQAlignedReads <- function(bamDir, bamFile, workDir=NULL){
    ### Loading the required package.
    #   This code was successfully tested with the package Rsamtools v.2.20.0.
    suppressMessages(expr=library(package=Rsamtools))
    ### Setting the work directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Full path to the BAM file.
    if (is.null(x=bamDir)){
        path <- paste(workDir, bamDir, sep="")
    }else{
        path <- paste(workDir, bamDir, sep="/")
    }
    ### Collection of the quality metrics.
    name <- gsub(pattern=".bam", replacement="", x=bamFile)
    qa <- scanBam(file=paste(path, bamFile, sep="/"),
                  index=paste(paste(path, bamFile, sep="/"), "bai", sep="."),
                  param=ScanBamParam(what=c("flag", "mapq", "isize")))
    numbers <- length(x=qa[[1]]$flag)
    mapq <- qa[[1]]$mapq
    mapq <- mapq[!is.na(x=mapq)]
    isize <- abs(x=qa[[1]]$isize)
    isize <- isize[!is.na(x=isize)]
    flag <- colSums(x=bamFlagAsBitMatrix(flag=qa[[1]]$flag))
    rm("qa")
    ### Quality metrics.
    return(list(SN=name,
                FN=bamFile,
                TNRs=numbers,
                MAPQs=mapq,
                TLENs=isize,
                FLAGs=flag)
           )
}
