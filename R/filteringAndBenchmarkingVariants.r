#' Filter of single nucleotide variations
#' @description This function performs multi-stage filtering of single
#'     nucleotide variations given in the standard VCF format.
#' @param vcfDir character string giving the name (or path to and name) of
#'     directory with VCF file(-s). NULL by default, which means the current
#'     working directory will be used instead of the specified VCF directory.
#' @param vcfFile character string giving the name of input VCF file.
#' @param gr character string specifying the name of the tab-delimited TXT file
#'     containing coordinates of genomic location(-s) of interest. The default
#'     value is NULL. If so, this file must contains the following four fields:
#'     i) seqnames - the name of chromosome or scaffold with prefix "chr";
#'     ii) start   - start coordinate of genomic interval of interest;
#'     iii) end    - end coordinate of genomic interval of interest;
#'     iv) strand  - (optionally) strand information about genomic interval
#'                   of interest.
#'     TXT file with genomic location(-s) must be in working directory.
#' @param depth integer value giving the threshold for the minimum sequencing
#'     depth. Default value is 1.
#' @param score numeric value that specifies the threshold for the QUAL field.
#'     Default value is NULL.
#' @param heterozyg numeric value that specifies the threshold for the total
#'     heterozygosity calculated by the fields MMsum and DP. Default value is 0.
#' @param heterozyg1 numeric value that specifies the threshold for the
#'     heterozygosity of first alternative allele calculated by the fields MM
#'     and DP. Default value is 0.
#' @param index logical agrument, whether to bgzip the output file and generate
#'     a tabix index.
#' @param workDir character string giving the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return filtered VCF file with postfix "filtered".
#' @author Vasily V. Grinev.
#' @examples
#' SNVs_filter <- filterVariants(vcfDir="Files_VCF",
#'                               vcfFile="example_seq.vcf",
#'                               gr="genomic_intervals.txt",
#'                               depth=10,
#'                               score=5,
#'                               heterozyg=0.2,
#'                               heterozyg1=0.2,
#'                               index=TRUE,
#'                               workDir="D:/Vasily Grinev")
#' @export
#' @importFrom VariantAnnotation readVcf writeVcf
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom methods setClassUnion
filterVariants <- function(vcfDir=NULL,
                           vcfFile,
                           gr=NULL,
                           depth=1,
                           score=NULL,
                           heterozyg=0,
                           heterozyg1=0,
                           index=FALSE,
                           workDir=NULL){
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
    ### Loading of the single nucleotide variations as an object of class
    #   CollapsedVCF.
    SNVs <- readVcf(file=paste(vcfDir, vcfFile, sep="/"))
    ### Multi-stage Filtering.
    if (!is.null(x=gr)){
        ##  Stage 1: Filter by genomic location(-s).
        #   Loading of the genomic location(-s) as an object of class GRanges.
        intervals <- read.table(file=paste(workDir, gr, sep="/"),
                                sep="\t",
                                header=TRUE,
                                quote="\"",
                                as.is=TRUE)
        intervals <- makeGRangesFromDataFrame(df=intervals)
        #   Filtering at stage 1.
        hits <- findOverlaps(query=rowRanges(x=SNVs),
                             subject=intervals,
                             type="within",
                             ignore.strand=TRUE)
        SNVs <- SNVs[queryHits(x=hits)]
    }
    if (!is.null(x=score)){
        ##  Stage 2: Filter by score.
        SNVs <- SNVs[rowRanges(x=SNVs)$QUAL >= score]
    }
    if (depth > 1){
        ##  Stage 3: Filter by sequencing depth.
        SNVs <- SNVs[info(x=SNVs)$DP >= depth]
    }
    if (heterozyg > 0){
        ##  Stage 4: Filter by heterozygosity.
        if ("MMsum" %in% colnames(info(x=SNVs))){
            SNVs <- SNVs[info(x=SNVs)$MMsum/info(x=SNVs)$DP >= heterozyg]
        }
    }
    if (heterozyg1 > 0){
        ##  Stage 5: Filter by heterozygosity of the first alternative allele.
        MM <- suppressWarnings(expr=as.numeric(x=info(x=SNVs)$MM))
        MM_na <- info(x=SNVs)$MM[is.na(x=MM)]
        MM_na <- as.numeric(x=unlist(x=lapply(X=MM_na,
                       FUN=function(y){max(x=strsplit(x=y, split=",")[[1]])})))
        MM[is.na(x=MM)] <- MM_na
        SNVs <- SNVs[MM/info(x=SNVs)$DP >= heterozyg1]
    }
    ### Saving and returning of the final object.
    writeVcf(obj=SNVs,
             filename=paste(vcfDir,
                            sub(pattern="vcf",
                                replacement="filtered.vcf",
                                x=vcfFile),
                            sep="/"),
             index=index)
    return(SNVs)
}

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
#' @importFrom data.table fread as.data.table data.table
benchmarkVariants <- function(vcfDir,
                              ref_vcf,
                              query_vcf,
                              featureTable=FALSE,
                              workDir=NULL){
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

