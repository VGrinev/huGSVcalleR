#' Group short reads in BAM file
#' @description This is a wrapper function that uses the low-level function
#'     AddOrReplaceReadGroup() from package GATK to add read group tags to BAM
#'     file of interest. This stage could be not nessesery, if BAM file already
#'     contains read group tags.
#' @param bamDir character string giving the name (or path to and name) of
#'     directory with BAM file(s). NULL by default, which means the current
#'     working directory will be used instead of the specified BAM directory.
#' @param bamFile character string giving the name of input BAM file.
#' @param postfix character string that will be added to the output file name
#'     as postfix. "RG" by default.
#' @param SE if TRUE (by default), standart error determined by function
#'     AddOrReplaceReadGroup() will be writen in file STDERR.RG.TXT.
#' @param RGID RGID value for appropriate BAM field.
#' @param RGLB RGLB value for appropriate BAM field.
#' @param RGPL RGPL value for appropriate BAM field.
#' @param RGPU RGPU value for appropriate BAM field.
#' @param RGSM RGSM value for appropriate BAM field.    
#' @param gatk_path character string giving the path to internal or external
#'     (system) GATK executive script. This path can be obtained by running
#'     Find.Gatk() or settled manually. Default value is "gatk", i.e. function
#'     will uses GATK systemly from $PATH.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @return a new BAM file with additional information about read groups and
#'     (optionally) auxillary file with standart error log.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @examples
#' @export
#' Last updated: July 23, 2025.

groupReadsGATK <- function(bamDir=NULL,
                           bamFile,
                           postfix="RG",
                           SE=TRUE,
                           RGID="Bam.Recal.RG",
                           RGLB="default.LIBRARY.ID",
                           RGPL="default.ILLUMINA",
                           RGPU="default.PLATFORM",
                           RGSM="default.SAMPLE",
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
    ### Naming of output file.
    output <- paste(gsub(pattern=".bam",
                         replacement="",
                         x=bamFile),
                    postfix,
                    "bam",
                    sep=".")
    ### Settings for saving of standart error.
    se <- ""
    if (SE == TRUE){
        se <- paste(gsub(pattern=".bam", replacement="", x=bamFile),
                    "STDERR.RG.TXT",
                    sep=".")
    }
    ### Running of command.
    system2(command=gatk_path,
            args=c("AddOrReplaceReadGroups",
                   paste("-I", paste(path, bamFile, sep="/"), sep=" "),
                   paste("-O", paste(path, output, sep="/"), sep=" "),
                   paste("--RGID", RGID, "--RGLB", RGLB, "--RGPL", RGPL,
                         "--RGPU", RGPU, "--RGSM", RGSM, sep=" ")),
            if (se == ""){
                stderr=se
            }else{
                stderr=paste(path, se, sep="/")
            })
}
