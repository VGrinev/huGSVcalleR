#' Call of single nucleotide variations using Fisher’s exact test
#' @description This is a wrapper function that uses the low-level function
#'     exactSNP() from package Rsubread for accurate and efficient SNVs calling.
#' @param bamDir character string giving the name (or path to and name) of
#'     directory with BAM file(-s). NULL by default, which means the current
#'     working directory will be used instead of the specified BAM directory.
#' @param bamFile character string giving the name of input BAM file.
#' @param fastaDir character string giving the name (or path to and name) of
#'     folder for storing of reference genome.
#' @param faFile charater string giving the name of file containing reference
#'     genome in format *.FASTA or *.FA.
#' @param baseq integer value, minimal "base quality" for each nucleotide in an
#'     alignment. Default value is 20 that means read nucleotides with quality
#'     scores less than 20 will be excluded from analysis.
#' @param mindepth integer value giving the minimal coverage of location
#'     containing the variation. Default value is 1.
#' @param maxdepth integer value giving the maximal coverage of location
#'     containing the variation. Default value is 1e6. This option is useful
#'     for removing PCR artifacts.
#' @param qvalue numeric value giving the q-value cutoff for variations calling
#'     at sequencing depth 50X. Default value is 12. The q-value is calcuated
#'     as -log10(p), where p is a p-value yielded from the Fisher’s exact test.
#'     The function exactSNP() automatically adjusts the q-value cutoff for
#'     each chromosomal location according to its sequencing depth, based on
#'     this cutoff.
#' @param trim integer value giving the number of nucleotides trimmed off from
#'     each end of the read. Default value is 0.
#' @param threads numeric value giving the number of threads/CPUs used. Default
#'     value is 1.
#' @param workDir character string giving the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return file of VCF format is written to the current working directory.
#' @author Vasily V. Grinev.
#' @examples
#' res <- fisherCalleR(bamDir="Files_BAM",
#'                     bamFile="SRR21721515.bam",
#'                     fastaDir="Files_FASTA",
#'                     faFile="hg38.fa",
#'                     baseq=20,
#'                     mindepth=1,
#'                     maxdepth=1e6,
#'                     qvalue=12,
#'                     trim=0,
#'                     threads=3,
#'                     workDir="D:/Vasily Grinev")
#' @export
#  Last updated: August 19, 2025.

fisherCalleR <- function(bamDir=NULL,
                         bamFile,
                         fastaDir=NULL,
                         faFile,
                         baseq=20,
                         mindepth=1, maxdepth=1e6,
                         qvalue=12,
                         trim=0,
                         threads=1,
                         workDir=NULL){
    ### Loading the required package.
    #   This code was successfully tested with the package Rsubread v.2.32.2.
    suppressMessages(expr=library(package=Rsubread))
    ### Setting the working directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Full path to the BAM file.
    if (is.null(x=bamDir)){
        path1 <- workDir
    }else{
        path1 <- paste(workDir, bamDir, sep="/")
    }
    ### Full path to the reference genome.
    if (is.null(x=fastaDir)){
        path2 <- workDir
    }else{
        path2 <- paste(workDir, fastaDir, sep="/")
    }
    ### Setting the BAM file.
    bam <- paste(path1, bamFile, sep="/")
    ### Setting the FASTA file.
    fa <- paste(path2, faFile, sep="/")
    ### Setting the output VCF file.
    vcf <- gsub(pattern=".bam",
                replacement=", fisherCalleR, SNVs.vcf",
                x=bam)
    ### Detection of SNVs.
    SNVs <- exactSNP(readFile=bam,
                     isBAM=TRUE,
                     refGenomeFile=fa,
                     minBaseQuality=baseq,
                     minReads=mindepth,
                     maxReads=maxdepth,
                     qvalueCutoff=qvalue,
                     nTrimmedBases=trim,
                     nthreads=threads,
                     outputFile=vcf)
    ### Returning a final object.
    return(SNVs)
}
