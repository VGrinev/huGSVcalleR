#' Call of single nucleotide variations using haplotype caller GATK
#' @description This is a wrapper function to run GATK within the R environment
#'     for accurate and efficient SNVs calling.
#' @param bamDir character string giving the name (or path to and name) of
#'     directory with BAM file(s). NULL by default, which means the current
#'     working directory will be used instead of the specified BAM directory.
#' @param bamFile character string giving the name of input BAM file.
#' @param postfix character string that will be added to the output VCF file
#'     name as postfix. NA by default at which the postfix "GATK based SNVs"
#'     will be used.
#' @param fastaDir character string giving the name (or path to and name) of
#'     folder for storing of fai-indexed reference genome.
#' @param faFile charater string giving the name of file containing reference
#'     genome in format *.FASTA or *.FA.
#' @param intervals character string specifying the name of the tab-delimited
#'     TXT file containing coordinates of genomic location(-s) of interest. The
#'     default value is NULL. If so, this file must contains the following four
#'     fields:
#'     i) seqnames - the name of chromosome or scaffold with prefix "chr";
#'     ii) start   - start coordinate of genomic interval of interest;
#'     iii) end    - end coordinate of genomic interval of interest;
#'     iv) strand  - (optionally) strand information about genomic interval
#'                   of interest.
#'     TXT file with genomic location(-s) must be in working directory.
#' @param pcr_model GATK argument --pcr-indel-model. Possible values: NONE,
#'     AGRESSIVE, HOSTILE, CONSERVATIVE. See HaplotypeCaller
#'     at https://gatk.broadinstitute.org for further details.
#' @param soft_clip GATK argument --dont-use-soft-clipped. TRUE by default.
#' @param ERC GATK ERC mode. "NONE" by default. Possible values: NONE, GVCF or
#'     BP_RESOLUTION. See HaplotypeCaller at https://gatk.broadinstitute.org
#'     for further details.
#' @param cores number of CPU cores to use for Hidden Markov Models. If "MAX"
#'     (by default), all available cores will be detected and used.
#' @param SE if TRUE (by default), standart error determined by GATK will be
#'     writen in file STDERR.HC.TXT.
#' @param gatk_path character string giving the path to internal or external
#'     (system) GATK executive script. This path can be obtained by running
#'     Find.Gatk() or settled manually. Default value is "gatk", i.e. function
#'     will uses GATK systemly from $PATH.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @details GATK performs SNVs and short indels calling by local de novo
#'     assembly of mapped reads
#' @return a new VCF file with detected SNVs.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @examples
#' @export
#' Last updated: July 26, 2025.

gatkCalleR <- function(bamDir=NULL,
                       bamFile,
                       postfix=NA,
                       fastaDir=NULL,
                       faFile,
                       intervals=NULL,
                       pcr_model,
                       soft_clip=TRUE,
                       ERC="NONE",
                       cores="MAX",
                       SE=TRUE,
                       gatk_path=NULL,
                       workDir=NULL){
    ### Loading the required package.
    #   This code was successfully tested with the package parallel.
    suppressMessages(expr=library(package=parallel))
    ### Validation of GATK path.
    if(is.null(x=gatk_path)){
        return("Cannot run GATK without path. Please provide full path to GATK
                folder and execution script! Format of 
                path: /home/Bioinformatics/Soft/gatk-4.2.6.1/gatk. Or use a
                just *gatk*, if it can be called directly from bash console")
    }
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
    ### Settings of run.
    #   Settings for input BAM file.
    input <- paste("-I",
                   paste(path1, bamFile, sep="/"),
                   sep=" ")
    #   Settings for output VCF file.
    if (!is.na(x=postfix)){
        output <- paste("-O",
                        paste(path1,
                              paste(gsub(pattern=".bam",
                                         replacement="",
                                         x=bamFile),
                                    postfix,
                                    "vcf",
                                    sep="."),
                              sep="/"),
                        sep=" ")
    }else{
        output <- paste("-O",
                        paste(path1,
                              paste(gsub(pattern=".bam",
                                         replacement="",
                                         x=bamFile),
                                    "GATK based SNVs",
                                    "vcf",
                                    sep="."),
                              sep="/"),
                        sep=" ")
    }
    #   Settings for fai-indexed reference genome.
    fa <- paste("-R", paste(path2, faFile, sep="/"), sep=" ")
    #   Settings for intervals.
    if (!is.null(x=intervals)){
        GRs <- read.table(file=paste(workDir, intervals, sep="/"),
                          sep="\t",
                          header=TRUE,
                          quote="\"",
                          as.is=TRUE)
       GRs <- paste(paste(GRs[, 1], GRs[, 2], sep=":"), GRs[, 3], sep="-")
       GRs <- paste("--intervals", GRs)
       fa <- paste(fa, GRs)
    }  
    #   Settings for PCR indel model.
    pcr <- paste("--pcr-indel-model", pcr_model)
    #   Settings for using for discovery the soft clipped regions
    softclip <- paste("--dont-use-soft-clipped-bases", soft_clip)
    # Settings for ERC output.
    erc <- paste("-ERC", ERC)
    #   Setting the multithreads process.
    if (cores == "MAX"){
        cores <- detectCores()
    }
    cores <- paste("--native-pair-hmm-threads", cores)
    ### Settings for saving of standart error.
    se <- ""
    if (SE == TRUE){
        se <- paste(gsub(pattern=".bam", replacement="", x=bamFile),
                    "STDERR.HC.TXT",
                    sep=".")
    }
    ### Running of command.
    system2(command=gatk_path,
            args=c("HaplotypeCaller", input, output, fa, pcr, softclip, erc, cores),
            if (se == ""){
                stderr=se
            }else{
                stderr=paste(path1, se, sep="/")
            })
}
