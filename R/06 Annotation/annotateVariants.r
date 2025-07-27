#' Annotate of single nucleotide variations
#' @description This function performs a multivariate annotation of the single
#'     nucleotide variations given in the standard VCF format.
#' @param vcfDir character string giving the name (or path to and name) of
#'     directory with VCF file(-s). NULL by default, which means the current
#'     working directory will be used instead of the specified VCF directory.
#' @param vcfFile character string giving the name of input VCF file.
#' @param ref_genes character string giving the name of TXT file in
#'     tab-delimited format with basic annotations of reference genes. This
#'     file must contains the following fields:
#'     i) gene_id        - (optionally) gene ID (for example, Ensembl based);
#'     ii) gene_name     - name of gene (for example, HUGO based);
#'     iii) gene_biotype - biotype of gene;
#'     iv) seqnames      - name of chromosome or scaffold with prefix "chr";
#'     v) start          - start coordinate of gene;
#'     vi) end           - end coordinate of gene;
#'     vii) strand       - (optionally) strand information about gene.
#' @param output character string giving the name (without extension) of TXT
#'     file in tab-delimited format for storing annotated single nucleotide
#'     variations. NULL by default, which means that the name will be generated
#'     automatically based on the name of the input VCF file.
#' @param workDir character string giving the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return an object of class data.frame and new TXT file in tab-delimited
#'     format containing annotated single nucleotide variations.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @examples
#' genes <- "Ensembl release 114, GRCh38.p14, annotations of genes.txt"
#' annoSNVs <- annotateVariants(vcfDir="Files_VCF",
#'                              vcfFile="example_seq.vcf",
#'                              ref_genes=,
#'                              output="example_seq, annotated SNVs",
#'                              workDir="D:/Vasily Grinev")
#' @export
#' Last updated: July 27, 2025.

annotateVariants <- function(vcfDir=NULL,
                             vcfFile,
                             ref_genes,
                             output=NULL,
                             workDir=NULL){
    ### Loading the required package.
    #   This code was successfully tested with the package
    #   VariantAnnotation v.1.55.0.
    suppressMessages(expr=library(package=data.table))
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
    vcf <- readVcf(file=paste(vcfDir, vcfFile, sep="/"))
    ### Loading of the reference genes as an object of class GRanges.
    genes <- read.table(file=paste(workDir, ref_genes, sep="/"),
                        sep="\t",
                        header=TRUE,
                        quote="\"",
                        as.is=TRUE)
    genes <- makeGRangesFromDataFrame(df=genes, keep.extra.columns=TRUE)
    ### Annotation of the single nucleotide variations with gene names.
    mcols(x=rowRanges(x=vcf))$gene_name <- ""
    hits <- findOverlaps(query=rowRanges(x=vcf),
                         subject=genes[!is.na(x=genes$gene_name)],
                         type="within",
                         ignore.strand=TRUE)
    hits <- data.table(cbind(hits@from,
                       genes[!is.na(x=genes$gene_name)][hits@to, ]$gene_name))
    hits <- hits[, paste(V2, collapse=","), by=V1]
    colnames(x=hits) <- c("V1", "V2")
    hits$V1 <- as.numeric(x=hits$V1)
    hits$V2 <- unlist(x=lapply(X=hits$V2,
                           FUN=function(y){paste0(sort(x=unique(x=strsplit(x=y,
                                                             split=",")[[1]])),
                                                 collapse=",")}))
    rowRanges(x=vcf)[hits$V1, ]$gene_name <- as.vector(x=hits$V2)
    ### Returning and saving of the final object.
    SNVs <- cbind(data.frame(rowRanges(x=vcf)), info(x=vcf))
    SNVs$ALT <- unlist(x=lapply(X=SNVs$ALT, paste0, collapse=","))
    SNVs$snv_id <- paste(paste(paste(SNVs$seqnames, SNVs$start, sep=":"),
                               SNVs$REF, sep="_"),
                         SNVs$ALT, sep="/")
    SNVs <- SNVs[, c("snv_id", "seqnames", "start", "end", "strand", "REF",
                     "ALT", "DP", "MM", "QUAL", "gene_name")]
    colnames(x=SNVs) <- c(colnames(x=SNVs)[1:5],
                          "ref", "alt", "pos_depth", "alt_depth",
                          "score", "gene_name")
    SNVs$seqnames <- as.character(x=SNVs$seqnames)
    SNVs$strand <- as.character(x=SNVs$strand)
    SNVs$alt_depth <- as.numeric(x=SNVs$alt_depth)
    if (is.null(x=output)){
        output_file <- paste(workDir,
                             sub(pattern="vcf",
                                 replacement=", annotated SNVs.txt",
                                 x=vcfFile),
                             sep="/")
    }else{
        output_file <- paste(workDir, paste(output, "txt", sep="."), sep="/")
    }
    write.table(x=SNVs,
                file=output_file,
                sep="\t",
                quote=FALSE,
                col.names=TRUE,
                row.names=FALSE)
    return(SNVs)
}
