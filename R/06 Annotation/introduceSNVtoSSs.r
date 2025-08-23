#' Introduce single nucleotide variant(-s) to splice site(-s)
#' @description Introduction of single nucleotide variant(-s) into the splice
#'     site sequence(-s).
#' @param x genomic coordinates of 5' or 3' splice sites. It can be the name of
#'     object of class GRanges (as output of extractSSsCoordsEEJs() or
#'     extractSSsCoordsGTF() function) or tab-delimited TXT file with the
#'     following fields:
#'     i) seqnames - name of chromosome or scaffold with prefix "chr";
#'     ii) start   - start coordinate of splice site;
#'     iii) end    - end coordinate of splice site;
#'     iv) strand  - strand information about splice site;
#'     v-...) ...  - (optionally) any metadata.
#'     If TXT file, it should be in work directory.
#' @param vcfDir character string giving the name (or path to and name) of
#'     directory with VCF file(-s). NULL by default, which means the current
#'     working directory will be used instead of the specified VCF directory.
#' @param vcfFile character string giving the name of input VCF file.
#' @param genome character string giving the name of the BSgenome data package
#'     to be used for reference sequence extraction. Default value is
#'     package BSgenome.Hsapiens.UCSC.hg38.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return object of class list containing extended data frame with splice
#'     sites and DNAStringSet with corrected splice site sequences.
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @examples
#' coords <- "EEJs.Kasumi_1.RR_KD.PMC_BSU.fiveSSs.txt"
#' vcf <- "test_seq4, fisherCalleR, SNVs.filtered.vcf"
#' intrSNVtoSSs <- introduceSNVtoSSs(x=coords,
#'                                   vcfDir="Files_VCF",
#'                                   vcfFile=vcf,
#'                                   genome="BSgenome.Hsapiens.UCSC.hg38",
#'                                   workDir="D:/Vasily Grinev")
#' @export
#  Last updated: August 23, 2025.

introduceSNVtoSSs <- function(x,
                              vcfDir=NULL,
                              vcfFile,
                              genome="BSgenome.Hsapiens.UCSC.hg38",
                              workDir=NULL){
    ### Loading of required package.
    #   This code was successfully tested with package VariantAnnotation v.1.55.0.
    suppressMessages(expr=library(package=VariantAnnotation))
    ### Loading the required auxiliary function.
    source(file=paste(workDir, "extractSeqSSs.r", sep="/"))
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
    #   GRanges.
    setClassUnion("ExpData", c("matrix", "SummarizedExperiment"))
    vcf <- rowRanges(x=readVcf(file=paste(vcfDir, vcfFile, sep="/")))[, "ALT"]
    ALT <- as.character(x=unlist(x=vcf$ALT))
    vcf <- vcf[rep(x=1:length(x=vcf), times=lengths(x=vcf$ALT)), ]
    vcf$ALT <- ALT
    ### Loading of the genomic coordinates of splice sites as an object of
    #   class GRanges.
    if (is.object(x=x)){
        spl.sites <- x
    }else{
        frt <- tools::file_ext(x=x)
        if (!frt %in% "txt"){
            stop("Invalid file format")
        }
        spl.sites <- read.table(file=paste(workDir, x, sep="/"),
                                sep="\t",
                                header=TRUE,
                                quote="\"",
                                as.is=TRUE)
        spl.sites <- makeGRangesFromDataFrame(df=spl.sites,
                                              keep.extra.columns=TRUE)
    }
    ### Find overlaps between splice sites and SNVs.
    hits <- findOverlaps(query=spl.sites,
                         subject=vcf,
                         type="any",
                         ignore.strand=TRUE)
    if (length(x=hits) > 0){
        spl.sites <- spl.sites[queryHits(x=hits)]
        vcf <- vcf[subjectHits(x=hits)]
    ### Extraction of reference sequences of splice sites to be corrected.
        spl.sites$seq_original <- extractSeqSSs(x=spl.sites,
                                                genome=genome,
                                                workDir=workDir)
        spl.sites$seq_original <- as.character(x=spl.sites$seq_original)
        spl.sites$seq_corrected <- spl.sites$seq_original
        rev_compl <- as.character(x=strand(x=spl.sites)) == "-"
        spl.sites$seq_corrected[rev_compl] <-
        reverseComplement(x=DNAStringSet(x=spl.sites$seq_corrected[rev_compl]))
    ### Correction of the splice site sequences.
        if (length(x=spl.sites) != length(x=vcf)){
            stop("The list of splice sites to be corrected and the list of ",
                 "relevant SNVs differ in length.")
        }
        substr(x=spl.sites$seq_corrected,
               start=start(x=vcf) - start(x=spl.sites) + 1,
               stop=start(x=vcf) - start(x=spl.sites) + 1) <- vcf$ALT
        spl.sites$seq_corrected[rev_compl] <-
        reverseComplement(x=DNAStringSet(x=spl.sites$seq_corrected[rev_compl]))
    ### Returning the final object.
        anno <- data.frame(spl.sites)
        anno <- anno[, c(1:7,
                         (ncol(x=anno) - 1):ncol(x=anno),
                         8:(ncol(x=anno) - 2))]
        seqSSs <- DNAStringSet(x=anno$seq_corrected)
        names(x=seqSSs) <- paste(paste(seqnames(x=spl.sites),
                                       paste(start(x=spl.sites),
                                             end(x=spl.sites),
                                             sep="-"),
                                       sep=":"),
                                 strand(x=spl.sites),
                                 sep="_str")
        return(list(annotation=anno, sequences=seqSSs))
    }else{
        warning("WARNING: No matches were found between splice sites and SNVs.")
    }
}
