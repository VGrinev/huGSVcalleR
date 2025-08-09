#' Recalibrate a BAM file
#' @description This is a wrapper function that uses the low-level functions
#'     BaseRecalibrator() and ApplyBQSR() from package GATK to recalibrate
#'     sequencing quality scores of bases in BAM file of interest.
#' @param bamDir character string giving the name (or path to and name) of
#'     directory with BAM file(s). NULL by default, which means the current
#'     working directory will be used instead of the specified BAM directory.
#' @param bamFile character string giving the name of input BAM file.
#' @param vcfDir character string giving the name (or path to and name) of
#'     directory with VCF file(s). NULL by default, which means that the
#'     current working directory will be used instead of the specified VCF
#'     directory.
#' @param vcfFile character string giving the name of input VCF file.
#' @param fastaDir character string giving the name (or path to and name) of
#'     folder for storing of reference genome.
#' @param faFile charater string giving the name of file containing reference
#'     genome in format *.FASTA or *.FA.
#' @param createBQSR logical, TRUE by default. It allows to create a new BQSR
#'     table with provided BAM file, known variation sites and reference genome.
#' @param recalibrate logical, TRUE by default. It allows to run a recalibration
#'     process. For this purpose, the BQSR table provided trough nameBQSR or
#'     created by function itself if createBQSR=TRUE, will be used.
#' @param nameBQSR character string giving the name of pre-computed BQSR table.
#'     NULL by default.
#' @param postfix character string that will be added to the output file name
#'     as postfix. "Recalibrated" by default.
#' @param SE if TRUE (by default), standart error determined by function
#'     MarkDuplicates() will be writen in file STDERR.RECALIBRATION.TXT.
#' @param gatk_path character string giving the path to internal or external
#'     (system) GATK executive script. This path can be obtained by running
#'     Find.Gatk() or settled manually. Default value is "gatk", i.e. function
#'     will uses GATK systemly from $PATH.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @return if createBQSR=TRUE, a new BQSR table in working directory;
#'     if recalibrate=TRUE, a new BAM file with recalibrated sequencing quality
#'     scores of bases in working directory; if SE=TRUE, a file with standart
#'     error log in working directory.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @examples
#' ...
#' @export
recalibrateBamGATK <- function(bamDir=NULL,
                               bamFile,
                               vcfDir=NULL,
                               vcfFile,
                               fastaDir=NULL,
                               faFile,
                               createBQSR=TRUE,
                               recalibrate=TRUE,
                               nameBQSR,
                               postfix="Recalibrated",
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
  ### Full path to the VCF file.
  if (is.null(x=vcfDir)){
    path2 <- workDir
  }else{
    path2 <- paste(workDir, vcfDir, sep="/")
  }
  ### Full path to the FASTA/FA file.
  if (is.null(x=fastaDir)){
    path3 <- workDir
  }else{
    path3 <- paste(workDir, fastaDir, sep="/")
  }
  ### Generation of BQSR table.
  if (isTRUE(x=createBQSR)){
    bqsr_name <- paste(gsub(pattern=".bam", replacement="", x=bamFile),
                       "bqsr.table",
                       sep=".")
    system2(command=gatk_path,
            args=c("BaseRecalibrator",
                   paste("-I",
                         paste(path, bamFile, sep="/"),
                         sep = " "),
                   paste("-O",
                         paste(path, bqsr_name, sep="/"),
                         sep = " "),
                   paste("--known-sites",
                         paste(path2, vcfFile, sep="/"),
                         sep = " "),
                   paste("-R",
                         paste(path3, faFile, sep="/"),
                         sep = " ")))
  }
  ### Recalibration process.
  if (isTRUE(x=recalibrate)){
    if (isTRUE(x=createBQSR)){
      bqsr_name <- paste(gsub(pattern=".bam", replacement="", x=bamFile),
                         "bqsr.table",
                         sep=".")
    }else{
      bqsr_name <- nameBQSR
    }
    ### Settings for saving of standart error.
    se <- ""
    if (SE == TRUE){
      se <- paste(gsub(pattern=".bam", replacement="", x=bamFile),
                  "STDERR.RECALIBRATION.TXT",
                  sep=".")
    }
    ### Running of command.
    system2(command=gatk_path,
            args=c("ApplyBQSR",
                   paste("-bqsr", paste(path, bqsr_name, sep="/"), sep=" "),
                   paste("-I", paste(path, bamFile, sep="/"), sep=" "),
                   paste("-O", paste(path,
                                     paste(gsub(pattern=".bam",
                                                replacement="",
                                                x=bamFile),
                                           postfix,
                                           "bam",
                                           sep="."),
                                     sep="/"), sep=" ")))
    if (se == ""){
      stderr=se
    }else{
      stderr=paste(path, se, sep="/")
    }
  }
}
