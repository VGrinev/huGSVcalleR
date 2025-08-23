#' Extract of splice site sequences from a reference genome
#' @description Extraction of splice site sequences from the reference genome.
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
#' @param genome character string giving the name of the BSgenome data package
#'     to be used for reference sequence extraction. Default value is
#'     package BSgenome.Hsapiens.UCSC.hg38.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return object of class DNAStringSet containing splice site sequences.
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @examples
#' res <- extractSeqSSs(x="EEJs.Kasumi_1.RR_KD.PMC_BSU.fiveSSs.txt",
#'                      genome="BSgenome.Hsapiens.UCSC.hg38",
#'                      workDir="D:/Vasily Grinev")
#' @export
#  Last updated: August 22, 2025.

extractSeqSSs <- function(x,
                          genome="BSgenome.Hsapiens.UCSC.hg38",
                          workDir=NULL){
    ### Loading the required packages.
    #   This code was successfully tested with the packages BSgenome v.1.72 and
    #   BSgenome.Hsapiens.UCSC.hg38 v.1.4.5.
    suppressMessages(expr=library(genome, character.only=TRUE))
    ### Setting the working directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Loading of the genomic coordinates of splice sites as an object of
    #   class GRanges.
    if (is.object(x=x)){
        coordsSSs <- x
    }else{
        frt <- tools::file_ext(x=x)
        if (!frt %in% "txt"){
            stop("Invalid file format")
        }
        coordsSSs <- read.table(file=paste(workDir, x, sep="/"),
                                sep="\t",
                                header=TRUE,
                                quote="\"",
                                as.is=TRUE)
        coordsSSs <- makeGRangesFromDataFrame(df=coordsSSs,
                                              keep.extra.columns=TRUE)
    }
    ### Extraction of splice site sequences from the reference genome.
    seqSSs <- getSeq(x=get(x=genome), names=coordsSSs)
    names(x=seqSSs) <- paste(paste(seqnames(x=coordsSSs),
                                   paste(start(x=coordsSSs),
                                         end(x=coordsSSs),
                                         sep="-"),
                                   sep=":"),
                             strand(x=coordsSSs),
                             sep="_str")
    ### Returning the final object.
    return(seqSSs)
}
