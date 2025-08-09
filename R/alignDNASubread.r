#' Local alignment of short DNA-Seq reads by seed-and-vote algorithm
#' @description This is a wrapper function that uses the low-level function
#'     align() from package Rsubread to align the short DNA-Seq reads by
#'     seed-and-vote algorithm.
#' @param genome character vector giving the path to and basename of index file.
#' @param fastqDir character vector giving the path to and name of directory
#'     with FASTQ file(-s).
#' @param fastq1 character vector including the name of FASTQ file containing
#'     DNA-Seq reads to be aligned. For paired-end reads, this gives the name
#'     of FASTQ file for first reads in each DNA-Seq library.
#' @param fastq2 character vector including the name of FASTQ file containing
#'     second DNA-Seq reads to be aligned. NULL by default.
#' @param bamDir character vector giving the path to and name of output
#'     directory for BAM file.
#' @param bamFile character vector giving the name of output BAM file.
#' @param orientation character string giving the orientation of the two reads
#'     from the same pair. Default value is "fr" (forward for the first read
#'     and reverse for second one).
#' @param threads numeric value giving the number of threads used for mapping.
#'     Default value is 1.
#' @param SV logical value indicating if structural variants will be detected
#'     during read mapping. Default value is FALSE.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @return BAM and VCF (if SV argument is TRUE) file(-s) are written to the
#'     bamDir directory.
#' @author Vasily V. Grinev
#' @examples
#' genome <- "Reference_Genomes/GRCh38/GRCh38"
#' fastq1 <- c("example_seq.read1.fastq.gz")
#' fastq2 <- c("example_seq.read1.fastq.gz")
#' fastqDir <- "Files_FASTQ"
#' bamDir <- "Files_BAM"
#' bam <- "example_seq.bam"
#' workDir <- "/mnt/data/grinev"
#' al <- alignDNASubread(genome=genome,
#'                       fastqDir=fastqDir, fastq1=fastq1, fastq2=fastq2,
#'                       bamDir=bamDir,
#'                       bamFile=bam,
#'                       orientation="fr",
#'                       threads=28,
#'                       SV=FALSE,
#'                       workDir=workDir)
#' @export
#' @importFrom Rsubread align
#' @importFrom Rsamtools indexBam
alignDNASubread <- function(genome,
                            fastqDir, fastq1, fastq2=NULL,
                            bamDir,
                            bamFile,
                            orientation="fr",
                            threads=1,
                            SV=FALSE,
                            workDir=NULL){
    ### Assignment the basename of index file.
    ref <- paste(workDir, genome, sep="/")
    ### Assignment the path to and name of FASTQ file(-s) including first reads
    #   in DNA-Seq library(-ies).
    fq1 <- paste(paste(workDir, fastqDir, sep="/"), fastq1, sep="/")
    ### Assignment the path to and name of FASTQ file(-s) including second reads
    #   in DNA-Seq library(-ies).
    if (is.null(x=fastq2) == FALSE){
        fq2 <- paste(paste(workDir, fastqDir, sep="/"), fastq2, sep="/")
    }else{
        fq2 <- NULL
    }
    ### Alignment of short DNA-Seq reads against reference genome.
    align <- align(index=ref,
                   readfile1=fq1,
                   readfile2=fq2,
                   type="dna",
                   input_format="gzFASTQ",
                   output_format="BAM",
                   output_file=paste(paste(workDir, bamDir, sep="/"),
                                     bamFile,
                                     sep="/"),
                   PE_orientation=orientation,
                   nthreads=threads,
                   unique=TRUE,
                   detectSV=SV)
}
