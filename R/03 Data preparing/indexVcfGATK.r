#' Index a VCF file with known variants
#' @description This is a wrapper function that uses the low-level function
#'     IndexFeatureFile() from package GATK to built a .tbi index for VCF file
#'     with known variants.
#' @param vcfDir character string giving the name (or path to and name) of
#'     directory with VCF file(s). NULL by default, which means that the
#'     current working directory will be used instead of the specified VCF
#'     directory.
#' @param vcfFile character string giving the name of input VCF file.
#' @param output_name character string giving the name of index. NULL by
#'     default, which means that the name of input VCF file will be uased.
#'     Important note: to be used correctly by the GATK, the VCF file and its
#'     index must have the same name.
#' @param gatk_path character string giving the path to internal or external
#'     (system) GATK executive script. This path can be obtained by running
#'     Find.Gatk() or settled manually. Default value is "gatk", i.e. function
#'     will uses GATK systemly from $PATH.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @return .tbi file in working directory.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @examples
#' @export
#  Last updated: July 23, 2025.

indexVcfGATK <- function(vcfDir=NULL,
                         vcfFile=NULL,
                         output_name=NULL,
                         gatk_path="gatk",
                         workDir=NULL){
    ### Setting the working directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Full path to the VCF file.
    if (is.null(x=vcfDir)){
        path <- workDir
    }else{
        path <- paste(workDir, vcfDir, sep="/")
    }
    ### Running of command.
    if (is.null(x=output_name)){
        system2(command=gatk_path,
                args=c("IndexFeatureFile -I",
                        paste(path, vcfFile, sep="/")))
    }else{
        system2(command=gatk_path,
                args=c("IndexFeatureFile -I",
                        paste(path, vcfFile, sep="/"),
                       "-O", output_name))
    }
}
