#' Benchmark called variants against reference variants
#' @description This function performs benchmarking of called variants against
#'    reference variants based on chromosome, position and allele identity.
#' @param vcfDir character string giving the name (or path to and name) of
#'     directory with VCF file(-s). NULL by default, which means the current
#'     working directory will be used instead of the specified VCF directory.
#' @param ref_vcf character string giving the name of VCF file with reference
#'    variants. Valid data format is *.vcf.
#' @param query_vcf character string giving the name of VCF file with called
#'    variants. Valid data format is *.vcf.
#' @param featureTable logical argument, whether to create the feature table
#'    that links each variant to its predicted classification outcome
#'    (TP - true positive, FP - false positive, FN - false negative). Default
#'    value is FALSE.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return consolidated data frame containing values of benchmarking metrics;
#'     in addition, if featureTable=TRUE, in the working directory is writing
#'     CSV file containing a detailed report of the benchmarking results.
#' @author Liudmila S. Varaniuk.
#' @examples
#' res <- benchmarkVariants(vcfDir="Files_VCF",
#'                          ref_vcf="test_truth.vcf",
#'                          query_vcf="test_query.vcf",
#'                          featureTable=FALSE,
#'                          workDir="D:/Vasily Grinev")
#' @export
#  Last updated: August 22, 2025.

benchmarkVariants <- function(vcfDir,
                              ref_vcf,
                              query_vcf,
                              featureTable=FALSE,
                              workDir=NULL){
    ### Loading the required package.
    #   This code was successfully tested with the package data.table v.1.15.4.
    suppressMessages(expr=library(package=data.table))
    ### Setting the working directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Setting the VCF directory.
    if (is.null(x=vcfDir)){
        vcfDir <- workDir
    }else{
        vcfDir <- paste(workDir, vcfDir, sep="/")
    }
    ### Defining of auxiliary functions.
    ##  Auxiliary function for reading VCF data into an object of class
    #   data.table.
    read_vcf <- function(input){
        if (is.character(x=input) && file.exists(x=input)){
            vcf <- fread(input, skip="#CHROM", header=TRUE, data.table=TRUE)
            colnames(x=vcf)[1:8] <- c("CHROM", "POS", "ID", "REF",
                                      "ALT", "QUAL", "FILTER", "INFO")
        }else if (is.data.frame(x=input)){
            vcf <- as.data.table(x=input)
        }else{
            stop("Input data must be a path to a file or data frame.")
        }
        vcf <- vcf[, .(ALT=unlist(x=strsplit(x=ALT, split=","))),
                     by=.(CHROM, POS=as.integer(x=POS), REF)]
        vcf[, type := ifelse(nchar(x=REF) == 1 & nchar(x=ALT) == 1,
              "SNP", "INDEL")]
        vcf[, key := paste0(CHROM, ":", POS, "_", REF, ">", ALT)]
        return(vcf)
    }
    ##  Auxiliary function for collecting of key metrics.
    get_metrics <- function(truth_keys, test_keys){
        truth_keys <- unique(truth_keys)
        test_keys  <- unique(test_keys)
        TP <- sum(test_keys %in% truth_keys)
        FP <- sum(!test_keys %in% truth_keys)
        FN <- sum(!truth_keys %in% test_keys)
        precision <- if ((TP + FP) == 0) NA else TP/(TP + FP)
        recall <- if ((TP + FN) == 0) NA else TP/(TP + FN)
        f1 <- if (is.na(x=precision) || is.na(x=recall) || (precision + recall) == 0)
              NA else 2 * precision * recall/(precision + recall)
        return(data.frame(TP=as.integer(x=TP),
                          FP=as.integer(x=FP),
                          FN=as.integer(x=FN),
                          Precision=round(x=precision, digits=4),
                          Recall=round(x=recall, digits=4),
                          F1=round(x=f1, digits=4)))
    }
    ### Loading of VCF data.
    truth <- read_vcf(paste(vcfDir, ref_vcf, sep="/"))
    test  <- read_vcf(paste(vcfDir, query_vcf, sep="/"))
    ### Developing a featureTable.
    if (isTRUE(x=featureTable)){
        merged <- merge(truth[, .(key, CHROM, POS, REF, ALT, type)],
                        test[, .(key, CHROM, POS, REF, ALT, type)],
                        by="key", all=TRUE, suffixes=c(".truth", ".test"))
        merged[, tag := fifelse(!is.na(x=CHROM.truth) & !is.na(x=CHROM.test), "TP",
                        fifelse(is.na(x=CHROM.truth) & !is.na(x=CHROM.test), "FP", "FN"))]
        clean <- data.frame(CHROM=ifelse(is.na(x=merged$CHROM.test),
                                         merged$CHROM.truth,
                                         merged$CHROM.test),
                            POS=ifelse(is.na(x=merged$POS.test),
                                       merged$POS.truth,
                                       merged$POS.test),
                            REF=merged$REF.test,
                            REF.truth=merged$REF.truth,
                            ALT=merged$ALT.test,
                            type=ifelse(is.na(x=merged$type.test),
                                        merged$type.truth,
                                        merged$type.test),
                            tag=merged$tag)
        outfile <- paste0(sub(pattern=".vcf", replacement="", x=query_vcf),
                          " vs ",
                          sub(pattern=".vcf", replacement="", x=ref_vcf),
                          ", benchmarking results.csv")
        write.csv(x=clean,
                  file=paste(workDir, outfile, sep="/"),
                  row.names=FALSE)
    }
    ### Development and return of the final object.
    metrics <- rbind(SNPs=get_metrics(truth_keys=truth[type == "SNP", key],  
                                      test_keys=test[type == "SNP", key]),
                     INDELs=get_metrics(truth_keys=truth[type == "INDEL", key],
                                        test_keys=test[type == "INDEL", key]),
                     ALL=get_metrics(truth_keys=truth$key, test_keys=test$key))
    return(metrics)
}
