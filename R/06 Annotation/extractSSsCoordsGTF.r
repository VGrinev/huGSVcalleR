#' Extract genomic coordinates of splice sites from annotations
#' @description Extraction of the genomic coordinates of splice sites using
#'     genomic/transcriptomic annotations.
#' @param x character string specifying the name of file with input data.
#'     Accepted formats: GFF2.5/GTF or GFF3 (the gz archive of GFF2.5/GTF or
#'     GFF3 file is also accepted). The file should include models of genes in
#'     general/generic feature format with mandatory field "transcript_id".
#' @param gtfDir character string specifying the name of the directory
#'     containing the file with input data. NULL by default that mean the
#'     current working directory.
#' @param workDir character string specifying the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return object of class list containing two GRanges objects with genomic
#'     coordinates of 5' and 3' splice sites.
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @examples
#' res <- extractSSsCoordsGTF(x="Ensembl_GRCh38.p14_release.114.gtf",
#'                            gtfDir="Files_GTF",
#'                            workDir="D:/Vasily Grinev")
#' @export
#  Last updated: August 22, 2025.

extractSSsCoordsGTF <- function(x, gtfDir=NULL, workDir=NULL){
    ### Loading of required package.
    #   This code was successfully tested with package rtracklayer v.1.52.1.
    suppressMessages(expr=library(package=rtracklayer))
    ### Setting the working directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Full path to the file.
    if (is.null(x=gtfDir)){
        path <- workDir
    }else{
        path <- paste(workDir, gtfDir, sep="/")
    }
    ### Retrieving and validation of the file extension.
    frt <- tools::file_ext(x=x)
    if (!frt %in% c("gtf", "gff3", "gz")){
        stop("Invalid file format")
    }
    ### Loading of exon-exon junctions as an object of class GRanges.
    EEJs <- import(con=paste(path, x, sep="/"))
    EEJs <- EEJs[gsub(pattern=".+exon",
                      replacement="exon",
                      x=EEJs$type) == "exon", ]
    if ("transcript_id" %in% colnames(x=elementMetadata(x=EEJs)) == "TRUE"){
        splitField <- "transcript_id"
    }else{
        splitField <- "Name"
    }
    EEJs <- makeGRangesListFromDataFrame(df=as.data.frame(x=EEJs),
                                         split.field=splitField,
                                         keep.extra.columns=FALSE)
    EEJs <- unlist(x=setdiff(x=range(x=EEJs), EEJs))
    start(x=EEJs) <- start(x=EEJs) - 1
    end(x=EEJs) <- end(x=EEJs) + 1
    EEJs$transcript_id <- names(x=EEJs)
    names(x=EEJs) <- NULL
    EEJs <- as.data.frame(EEJs, stringsAsFactors=FALSE)
    EEJs$seqnames <- as.character(x=EEJs$seqnames)
    if (length(x=unique(x=nchar(EEJs$seqnames) == 1)) == 2){
        EEJs$seqnames <- paste("chr", EEJs$seqnames, sep="")
        if ("chrMT" %in% EEJs$seqnames == "TRUE"){
            EEJs[EEJs$seqnames == "chrMT", ]$seqnames <- "chrM"
        }
        EEJs <- EEJs[EEJs$seqnames %in%
                     c(paste("chr", c(1:22, "X", "Y", "M"), sep="")), ]
    }else{
        if ("chrMT" %in% EEJs$seqnames == "TRUE"){
            EEJs[EEJs$seqnames == "chrMT", ]$seqnames <- "chrM"
        }
        EEJs <- EEJs[EEJs$seqnames %in%
                     c(paste("chr", c(1:22, "X", "Y", "M"), sep="")), ]
    }
    EEJs$strand <- as.character(x=EEJs$strand)
    EEJs <- EEJs[!duplicated(x=EEJs), ]
    rownames(x=EEJs) <- NULL
    EEJs <- makeGRangesFromDataFrame(df=EEJs,
                                     keep.extra.columns=TRUE)
    ### Converting of genomic coordinates of exon-exon junctions
    #   in genomic coordinates of splice sites.
    fiveSSs <- EEJs
    start(x=fiveSSs[strand(x=fiveSSs) == "+", ]) <-
                               start(x=fiveSSs[strand(x=fiveSSs) == "+", ]) - 2
    end(x=fiveSSs[strand(x=fiveSSs) == "+", ]) <-
                               start(x=fiveSSs[strand(x=fiveSSs) == "+", ]) + 8
    start(x=fiveSSs[strand(x=fiveSSs) == "-", ]) <-
                               end(x=fiveSSs[strand(x=fiveSSs) == "-", ]) - 6
    end(x=fiveSSs[strand(x=fiveSSs) == "-", ]) <-
                               start(x=fiveSSs[strand(x=fiveSSs) == "-", ]) + 8
    threeSSs <- EEJs
    start(x=threeSSs[strand(x=threeSSs) == "+", ]) <-
                              end(x=threeSSs[strand(x=threeSSs) == "+", ]) - 20
    end(x=threeSSs[strand(x=threeSSs) == "+", ]) <-
                            start(x=threeSSs[strand(x=threeSSs) == "+", ]) + 22
    start(x=threeSSs[strand(x=threeSSs) == "-", ]) <-
                             start(x=threeSSs[strand(x=threeSSs) == "-", ]) - 2
    end(x=threeSSs[strand(x=threeSSs) == "-", ]) <-
                            start(x=threeSSs[strand(x=threeSSs) == "-", ]) + 22
    SSs <- list(fiveSSs=sort(x=fiveSSs), threeSSs=sort(x=threeSSs))
    ### Returning the final object.
    return(SSs)
}
