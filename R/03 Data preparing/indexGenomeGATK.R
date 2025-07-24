#' Index a reference genome
#' @description This is a wrapper function that uses the low-level function
#'     CreateSequenceDictionary() from package GATK to built .dict and/or .fai
#'     indices for reference genome of interest.
#' @param fastaDir character string giving the name (or path to and name) of
#'     folder for storing of reference genome.
#' @param faFile charater string giving the name of file containing reference
#'     genome in format *.FASTA or *.FA.
#' @param fai logical. If TRUE (by default), .fai index of FASTA file will be
#'     created in working diectory.
#' @param dict logical. If TRUE (by default), .dict index of FASTA file will be
#'     created in working diectory.
#' @param gatk_path character string giving the path to internal or external
#'     (system) GATK executive script. This path can be obtained by running
#'     Find.Gatk() or settled manualy. Default value is "gatk", i.e. function
#'     will uses GATK systemly from $PATH.
#' @param workDir character string giving the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return index files in working directory.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @examples
#' @export
#' Last updated: July 23, 2025.

indexGenomeGATK <- function(fastaDir=NULL,
                            faFile,
                            fai=TRUE,
                            dict=TRUE,
                            gatk_path="gatk",
                            workDir=NULL){
    ### Loading the required package.
    #   This code was successfully tested with the package Rsamtools v.2.20.0.
    suppressMessages(expr=library(package=Rsamtools))
    ### Setting the working directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Full path to the FASTA/FA file.
    if (is.null(x=fastaDir)){
        path <- workDir
    }else{
        path <- paste(workDir, fastaDir, sep="/")
    }
    ### Building the .fai index.
    if (fai == TRUE){
        indexFa(file=paste(path, faFile, sep="/")
    }
    ### Building the .dict index.
    if (dict == TRUE){
        system2(command=gatk_path,
                args=c("CreateSequenceDictionary",
                       "-R",
                       paste(path, faFile, sep="/")))
    }
}