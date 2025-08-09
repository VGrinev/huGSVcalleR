#' Detect and mark duplicated reads
#' @description This is a wrapper function that uses the low-level function
#'     MarkDuplicates() from package GATK to find out and mark duplicated reads
#'     in BAM file(-s) of interest.
#' @param bamDir character string giving the name (or path to and name) of
#'     directory with BAM file(s). NULL by default, which means the current
#'     working directory will be used instead of the specified BAM directory.
#' @param bamFile character string giving the name of input BAM file.
#' @param postfix character string that will be added to the output file name
#'     as postfix. "markedDuplicates" by default.
#' @param removeDupls a way to deal with detected duplicated reads. With "no"
#'     all duplicated reads will be kept in output BAM file, with "yes" all
#'     duplicated reads will be removed, and with "seq" only optical and other
#'     sequencing duplicates will be removed.
#' @param SE if TRUE (by default), standart error determined by function
#'     MarkDuplicates() will be writen in file STDERR.MARKDUP.TXT.
#' @param gatk_path character string giving the path to internal or external
#'     (system) GATK executive script. This path can be obtained by running
#'     Find.Gatk() or settled manually. Default value is "gatk", i.e. function
#'     will uses GATK systemly from $PATH.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @return a new BAM file with processed duplicated reads, a file with metrics
#'     of input BAM file and (optionally) a file with standart error log.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @examples
#' ...
#' @export
markDuplicatesGATK <- function(bamDir=NULL,
                               bamFile,
                               postfix="markedDuplicates",
                               removeDupls="no",
                               SE=TRUE,
                               gatk_path="gatk",
                               workDir=NULL){
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
    ### Settings of run.
    input <- paste("-I",
                   paste(path, bamFile, sep="/"),
                   sep=" ")
    output <- paste("-O",
                    paste(path,
                          paste(gsub(pattern=".bam",
                                     replacement="",
                                     x=bamFile),
                                postfix,
                                "bam",
                                sep="."),
                          sep="/"),
                    sep=" ")
    metrics <- paste("-M",
                     paste(path,
                           paste(gsub(pattern=".bam",
                                      replacement="",
                                      x=bamFile),
                                 postfix,
                                 "metrics.txt",
                                 sep="."),
                           sep="/"),
                     sep=" ")
    ### Settings for removing duplicated reads.
    if (removeDupls == "no"){
        remove.d <- paste("--REMOVE_DUPLICATES", "FALSE")
    }else{
        if (removeDupls == "yes"){
            remove.d <- paste("--REMOVE_DUPLICATES", "TRUE")
        }else{
            if (removeDupls == "seq"){
                remove.d <- paste("--REMOVE_DUPLICATES", "FALSE",
                                  "--REMOVE_SEQUENCING_DUPLICATES", "TRUE")
            }
        }
    }
    ### Settings for saving of standart error.
    se <- ""
    if (SE == TRUE){
        se <- paste(gsub(pattern=".bam", replacement="", x=bamFile),
                    "STDERR.MARKDUP.TXT",
                    sep=".")
    }
    ### Running of command.
    system2(command=gatk_path,
            args=c("MarkDuplicates", input, output, metrics, remove.d),
            if (se == ""){
                stderr=se
            }else{
                stderr=paste(path, se, sep="/")
            })
}
