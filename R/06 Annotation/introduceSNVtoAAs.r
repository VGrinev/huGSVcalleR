#' Predict amino acid coding changes for single nucleotide variants
#' @description Prediction of amino acid coding changes for single nucleotide
#'     variants.
#' @param vcfDir character string giving the name (or path to and name) of
#'     directory with VCF file(-s). NULL by default, which means the current
#'     working directory will be used instead of the specified VCF directory.
#' @param vcfFile character string giving the name of input VCF file.
#' @param OrgDb character string giving the name of the OrgDb data package with
#'     genome wide annotations. Default value is package org.Hs.eg.db.
#' @param TxDb character string giving the name of the TxDb data package with
#'     transcriptome wide annotations. Default value is package
#'     TxDb.Hsapiens.UCSC.hg38.knownGene.
#' @param genome character string giving the name of the BSgenome data package
#'     to be used for reference sequence extraction. Default value is
#'     package BSgenome.Hsapiens.UCSC.hg38.
#' @param proteins character string giving the name of TXT file in tab-delimited
#'     format with protein coding annotations of Ensembl transcripts. This file
#'     must contains the following IDs: Ensembl transcript IDs, Ensembl gene
#'     IDs, Ensembl protein _IDs, PDB IDs, UniProtKB/Swiss-Prot IDs and
#'     UniProtKB/TrEMBL IDs.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return object of class data frame with predicted and annotated amino acid
#'     coding changes for single nucleotide variants.
#' @author Vasily V. Grinev.
#' @examples
#' vcfFile="test_seq5, fisherCalleR, SNVs.vcf"
#' proteins="Ensembl_GRCh38.p14_release.114_proteins.txt"
#' intrSNVtoAAs <- introduceSNVtoAAs(vcfDir="Files_VCF",
#'                                   vcfFile=vcfFile,
#'                                   OrgDb="org.Hs.eg.db",
#'                                   TxDb="TxDb.Hsapiens.UCSC.hg38.knownGene",
#'                                   genome="BSgenome.Hsapiens.UCSC.hg38",
#'                                   proteins=proteins,
#'                                   workDir="D:/Vasily Grinev/SNVs")
#' @export
#  Last updated: August 26, 2025.

introduceSNVtoAAs <- function(vcfDir=NULL,
                              vcfFile,
                              OrgDb="org.Hs.eg.db",
                              TxDb="TxDb.Hsapiens.UCSC.hg38.knownGene",
                              genome="BSgenome.Hsapiens.UCSC.hg38",
                              proteins,
                              workDir=NULL){
    ### Loading of required packages.
    #   This code was successfully tested with packages
    #   BSgenome.Hsapiens.UCSC.hg38 v.1.4.5, org.Hs.eg.db v.3.21.0,
    #   TxDb.Hsapiens.UCSC.hg38.knownGene v.3.21.0 & VariantAnnotation v.1.55.0.
    suppressMessages(expr=library(package=genome, character.only=TRUE))
    suppressMessages(expr=library(package=OrgDb, character.only=TRUE))
    suppressMessages(expr=library(package=TxDb, character.only=TRUE))
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
    #   GRanges.
    setClassUnion("ExpData", c("matrix", "SummarizedExperiment"))
    vcf <- readVcf(file=paste(vcfDir, vcfFile, sep="/"))
    ALT <- as.character(x=unlist(x=rowRanges(x=vcf)$ALT))
    vcf <- vcf[rep(x=1:length(x=vcf), times=lengths(x=rowRanges(x=vcf)$ALT)), ]
    alt(x=vcf) <- DNAStringSetList(strsplit(x=ALT, split=""))
    ### Loading gene models as an object of class TxDb.
    txDb <- get(x=TxDb)
    ### Loading genome wide annotations as an object of class OrgDb.
    orgDb <- get(x=OrgDb)
    ### Loading of the protein annotations as an object of class data frame.
    pr <- read.table(file=paste(workDir, proteins, sep="/"),
                     sep="\t",
                     header=TRUE,
                     quote="\"",
                     as.is=TRUE)
    rownames(x=pr) <- pr[, 1]
    ### Prediction of amino acid coding changes for single nucleotide variants.
    AAs <- predictCoding(query=vcf, subject=txDb, seqSource=get(x=genome))
    AAs <- AAs[AAs$CONSEQUENCE == "nonsynonymous", ]
    ### Annotation.
    AAs$snv_id <- names(x=AAs)
    AAs$snv_ref <- as.character(x=AAs$REF)
    AAs$snv_alt <- as.character(x=unlist(x=AAs$ALT))
    AAs$snv_qual <- AAs$QUAL
    AAs$Ensembl_gene_IDs <- select(x=orgDb,
                                   keys=AAs$GENEID,
                                   columns="ENSEMBL",
                                   keytype="ENTREZID")$ENSEMBL
    AAs$Entrez_gene_IDs <- AAs$GENEID
    AAs$HUGO_gene_symbols <- select(x=orgDb,
                                    keys=AAs$GENEID,
                                    columns="SYMBOL",
                                    keytype="ENTREZID")$SYMBOL
    AAs$Ensembl_tx_names <- select(x=txDb,
                                   keys=AAs$TXID,
                                   columns="TXNAME",
                                   keytype="TXID")$TXNAME
    AAs$Ensembl_tx_names <- substr(x=AAs$Ensembl_tx_names, start=1, stop=15)
    AAs$SNV_CDS_loc <- start(x=AAs$CDSLOC)
    AAs$ref_codon <- as.vector(x=as.character(x=AAs$REFCODON))
    AAs$alt_codon <- as.vector(x=as.character(x=AAs$VARCODON))
    AAs$protein_id <- ""
    AAs$pdb_id <- ""
    AAs$uniprot.swissprot_id <- ""
    AAs$uniprot.trembl_id <- ""
    AAs$SNV_prot_loc <- as.vector(x=unlist(x=AAs$PROTEINLOC))
    AAs$ref_aa <- as.vector(x=as.character(x=AAs$REFAA))
    AAs$alt_aa <- as.vector(x=as.character(x=AAs$VARAA))
    AAs <- data.frame(AAs)[, c(25:28, 1:5, 29:42)]
    AAs$seqnames <- as.character(x=AAs$seqnames)
    AAs$strand <- as.character(x=AAs$strand)
    idx <- AAs[, 13] %in% pr[, 1]
    AAs[idx, 17:20] <- pr[AAs[idx, 13], 3:6]
    ### Returning the final object. 
    return(AAs)
}
