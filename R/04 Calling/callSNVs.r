#' Call single nucleotide variations using different statistical tests
#' @description This function is a wrapper that allows to call and control
#'     various low-level functions for identification of single nucleotide
#'     variations from NGS data.
#' @param fr_data either a name of CSV file containing pileup statistics (as
#'     returned by collectBasesPileup() function) or an object of class
#'     data.frame with frequancy table.
#' @param criteria character vector specifying which test(-s) to be performed.
#'     The following values are possible: "binomial", "counting", "entropy"
#'     (by default), "fisher", "fisherSmyth", "gatk" and "poisson".
#' @param vcfDir character string giving the name (or path to and name) of
#'     directory for storing of generated VCF file(-s). NULL by default, which
#'     means the current working directory will be used instead of the
#'     specified BAM directory.
#' @param outputVCF character string giving the name of output VCF file.
#' @param bamDir character string giving the name (or path to and name) of
#'     directory with BAM file(-s). NULL by default, which means the current
#'     working directory will be used instead of the specified BAM directory.
#' @param bamFile character string giving the name of input BAM file.
#' @param fastaDir character string giving the name (or path to and name) of
#'     folder for storing of reference genome (with .fai indices if GATK).
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
#' @param postfix character string that will be added to the output VCF file
#'     name as postfix. NA by default at which the postfix "GATK based SNVs"
#'     will be used. Only with GATK.
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
#' @return VCF file(-s) with detected single nucleotide variations.
#' @author Dzianis D. Sarnatski, Mikalai M. Yatskou, Vasily V. Grinev
#' @examples
#' @export
#' Last updated: August 14, 2025.

callSNVs <- function(fr_data, 
                     criteria="entropy",
                     vcfDir=NULL,
                     outputVCF,
                     #common arguments to be passed to "fisherSmyth" and "gatk"
                     bamDir=NULL,
                     bamFile=NA,
                     fastaDir=NULL,
                     faFile=NA,
                     #specific arguments to be passed to "fisherSmyth"
                     baseq=20,
                     mindepth=1,
                     maxdepth=1e6,
                     qvalue=12,
                     trim=0,
                     threads=1,
                     #specific arguments to be passed to "gatk"
                     postfix=NA,
                     intervals=NULL,
                     pcr_model=NA,
                     soft_clip=TRUE,
                     ERC="NONE",
                     cores="MAX",
                     SE=TRUE,
                     gatk_path=NULL,
                     workDir=NULL){
    ### Loading the required package.
    #   This code was successfully tested with the package data.table v.1.15.4.
    suppressMessages(expr=library(package=data.table))
    ### Validation of the input criteria.
    if (length(x=criteria) == 0){
        stop("No test(-s) will be performed. Please specify valid criteria")
    }
    invalid_criteria <- setdiff(x=criteria, y=SUPPORTED_CLASSICAL_TESTS)
    if (length(x=invalid_criteria) > 0){
        stop(paste("Invalid criteria specified:",
                   paste(invalid_criteria, collapse=", ")))
    }
    ### Setting the working directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Classical tests.
    if (criteria %in% c("binomial", "counting", "entropy", "fisher", "poisson")){
        ##  Retrieving of the frequancy table.
        if (is.object(x=fr_data)){
            fr_table=fr_data
        }else{
        #   Full path to the CSV file.
            path <- paste(workDir, fr_data, sep="/")
        #   Retrieving and validation of the file extension.
            frt <- tools::file_ext(x=path)
            if (!frt %in% "csv"){
                stop("Invalid format of CSV file")
            }
        #   Does the file actually exist in the specified location?
            if (!file.exists(path)){
                stop(paste("File not found:", path))
            }
        #   Loading of the CSV file as an object of class data.frame.
            fr_table <- read.csv(file=path, header=TRUE)
        }
        ##  Initialization.
        results <- list()
        ##  Perform tests based on the specified criteria.
        tests <- list(counting=snp_counting_cpp_parallel,
                      binomial=snp_binomial_cpp_parallel,
                      entropy=snp_entropy_cpp_parallel,
                      fisher=snp_fisher_cpp_parallel,
                      poisson=snp_poisson_cpp_parallel)
        results <- lapply(X=criteria,
                          FUN=function(test){tests[[test]](fr_table)})
        names(x=results) <- criteria
        ##  Transformation and writing of the results to VCF file.
        input <- data.table(cbind(fr_table, results))
        source(file=paste(workDir, "generateVCF.r", sep="/"))
        generateVCF(input=input,
                    version="4.0",
                    vcfDir=vcfDir,
                    outputVCF=outputVCF,
                    workDir=workDir)
    }
    ### Gordon K. Smyth's based approach.
    if (criteria == "fisherSmyth"){
        source(file=paste(workDir, "fisherCalleR.r", sep="/"))
        fisherSmythSNVs <- fisherCalleR(bamDir=bamDir,
                                        bamFile=bamFile,
                                        fastaDir=fastaDir,
                                        faFile=faFile,
                                        baseq=baseq,
                                        mindepth=mindepth,
                                        maxdepth=maxdepth,
                                        qvalue=qvalue,
                                        trim=trim,
                                        threads=threads,
                                        workDir=workDir)
    }
    ### GATK approach.
    if (criteria == "gatk"){
        source(file=paste(workDir, "gatkCalleR.r", sep="/"))
        gatkSNVs <- gatkCalleR(bamDir=bamDir,
                               bamFile=bamFile,
                               postfix=postfix,
                               fastaDir=fastaDir,
                               faFile=faFile,
                               intervals=NULL,
                               pcr_model,
                               soft_clip=TRUE,
                               ERC="NONE",
                               cores="MAX",
                               SE=TRUE,
                               gatk_path=NULL,
                               workDir=workDir)
    }
}
