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
#  Last updated: July 28, 2025.

filterVariants <- function(vcfDir=NULL,
                           vcfFile,
                           gr=NULL,
                           depth=1,
                           score=NULL,
                           heterozyg=0,
                           heterozyg1=0,
                           index=FALSE,
                           workDir=NULL){
    ### Loading the required package.
    #   This code was successfully tested with the package
    #   VariantAnnotation v.1.55.0.
    suppressMessages(expr=library(package=VariantAnnotation))
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
    setClassUnion("ExpData", c("matrix", "SummarizedExperiment"))
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
