#' Extract genomic coordinates of splice sites from exon-exon junctions
#' @description Extraction of the genomic coordinates of splice sites using
#'     experimentally detected exon-exon junctions.
#' @param x character string specifying the name of the tab-delimited TXT file
#'     with experimentally detected exon-exon junctions. This file must contain
#'     the following fields:
#'     i) eej_id     - exon-exon junction ID;
#'     ii) gene_id   - gene ID;
#'     iii) seqnames - name of chromosome or scaffold with prefix "chr";
#'     iv) start     - start coordinate of exon-exon junction (genomic
#'                     position of the last nucleotide in the upstream exon);
#'     v) end        - end coordinate of exon-exon junction(genomic
#'                     position of the first nucleotide in the downstream exon);
#'     vi) strand    - strand information about exon-exon junction location;
#'     vii-...) ...  - (optionally) one or more samples with summarized raw
#'                     RNA-Seq reads spanning exon-exon junction.
#' @param eejDir character string specifying the name of the directory
#'     containing the tab-delimited TXT file(-s) with experimentally detected
#'     exon-exon junctions. NULL by default that mean the current working
#'     directory.
#' @param thr an integer argument. It is a threshold for minimal coverage
#'     (sequencing depth) of exon-exon junction (in raw reads). Default value
#'     is 5 but can be specified by the user.
#' @param workDir character string specifying the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return object of class list containing two GRanges objects with genomic
#'     coordinates of 5' and 3' splice sites.
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @examples
#' res <- extractSSsCoordsEEJs(x="EEJs.Kasumi_1.RR_KD.PMC_BSU.txt",
#'                             eejDir="Files_EEJs",
#'                             workDir="D:/Vasily Grinev")
#' @export
#  Last updated: August 22, 2025.

extractSSsCoordsEEJs <- function(x, eejDir=NULL, thr=5, workDir=NULL){
    ### Loading of required package.
    #   This code was successfully tested with package GenomicRanges v.1.56.1.
    suppressMessages(expr=library(package=GenomicRanges))
    ### Setting the working directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Full path to the file.
    if (is.null(x=eejDir)){
        path <- workDir
    }else{
        path <- paste(workDir, eejDir, sep="/")
    }
    ### Retrieving and validation of the file extension.
    frt <- tools::file_ext(x=x)
    if (!frt %in% "txt"){
        stop("Invalid file format")
    }
    ### Loading of exon-exon junctions as an object of class GRanges.
    EEJs <- read.table(file=paste(path, x, sep="/"),
                       sep="\t",
                       header=TRUE,
                       quote="\"",
                       as.is=TRUE)
    if (ncol(x=EEJs) > 7 & !is.null(x=thr)){
        EEJs <- EEJs[rowSums(x=EEJs[, -1:-7] >= thr) >= 1, ]
    }
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
    SSs <- list(fiveSSs=fiveSSs, threeSSs=threeSSs)
    ### Returning the final object.
    return(SSs)
}
