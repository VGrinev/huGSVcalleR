#' Sort and index BAM file(-s)
#' @description This is a wrapper function that uses the low-level function
#'     sortBam() from package Rsamtools to sort and index BAM file(-s).
#' @param bamDir character vector giving the path to and name of directory
#'     with BAM file(s).
#' @param bamFile character vector giving the name(-s) of BAM file(-s) to be
#'     sorted and indexed.
#' @param byQname logical. It indicats whether the sorted destination file
#'     should be sorted by query-name (TRUE) or by mapping position (FALSE).
#' @param workDir character string giving the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return sorted and indexed BAM file(-s).
#' @author Vasily V. Grinev
#' @examples
#' sortBamFile(bamDir="Files_BAM",
#'             bamFile="SRR21721515.bam",
#'             byQname=TRUE,
#'             workDir="D:/Vasily Grinev")
#' @export
#' Last updated: July 23, 2025.

sortBamFile <- function(bamDir, bamFile, byQname, workDir=NULL){
    ### Loading the required packages.
    #   This code was successfully tested with the packages parallel and
    #   Rsamtools v.2.20.0.
    suppressMessages(expr=library(package=parallel))
    suppressMessages(expr=library(package=Rsamtools))
    ### Sorting and indexing of generated BAM file(-s).
    for (i in 1:length(x=bamFile)){
    sortedBam <- sortBam(file=paste(paste(workDir, bamDir, sep="/"),
                                    bamFile[i],
                                    sep="/"),
                         destination=sub(pattern=".bam",
                                         replacement="",
                                         x=paste(paste(workDir,
                                                       bamDir,
                                                       sep="/"),
                                                 bamFile[i],
                                                 sep="/")),
                         byQname=byQname,
                         nThreads=detectCores(logical=TRUE) - 1)
    suppressMessages(expr=indexBam(files=sortedBam))
    }
}
