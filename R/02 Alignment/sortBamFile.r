#' Sort and index BAM file(-s)
#' @description This is a wrapper function that uses the low-level function
#'     sortBam() from package Rsamtools to sort and index BAM file(-s).
#' @param bamDir character vector giving the path to and name of directory
#'     with BAM file(s).
#' @param bamFile character vector giving the name(-s) of BAM file(-s) to be
#'     sorted and indexed.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return sorted and indexed BAM file(-s).
#' @author Vasily V. Grinev
#' @examples
#' sortBamFile(bamDir="Files_BAM",
#'             bamFile=c("example_seq.bam"),
#'             workDir="/mnt/data/grinev")
#' @export
#' Last updated: July 18, 2025.

sortBamFile <- function(bamDir, bamFile, workDir=NULL){
    ### Loading the required package.
    #   This code was successfully tested with the package Rsamtools v.2.20.0.
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
                         byQname=FALSE,
                         maxMemory=40000)
    suppressMessages(expr=indexBam(files=sortedBam))
    }
}
