#' Score of splice sites
#' @description Scoring of the splice sites according to MaxEntScan models
#'     which are available at http://hollywood.mit.edu/burgelab/maxent or
#'     https://github.com/matthdsm/MaxEntScan.
#' @param x splice site sequences in FASTA format. It can be the name of object
#'     of class DNAStringSet (as output, for example, of extractSeqSSs()
#'     function) or standard FA/FASTA file.
#' @param type integer value 5 or 3 for define of donor or acceptor splice
#'     sites, respectively.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return object of class data frame containing three fields: IDs of splice
#'     sites, splice site sequences and MaxEntScan scores.
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @examples
#' res <- scoreSSs(x="EEJs.Kasumi_1.RR_KD.PMC_BSU.fiveSSs.fasta",
#'                 type=5,
#'                 workDir="D:/Vasily Grinev")
#' @export
#  Last updated: August 22, 2025.

scoreSSs <- function(x,
                     type,
                     workDir=NULL){
    ### Loading the required package.
    #   This code was successfully tested with the package Biostrings v.2.72.1.
    suppressMessages(expr=library(package=Biostrings))
    ### Connection to splice site scoring models.
    # models <- system.file("extdata", "me2x5", package="huGSVCalleR")
    # setwd(dir=dirname(path=models))
    models <- "D:/Vasily Grinev/exdata/me2x5"
    setwd(dir=dirname(path=models))
    ### Setting the working directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Loading of the splice site sequences as an object of class DNAStringSet.
    if (is.object(x=x)){
        seqSSs <- x
    }else{
        frt <- tools::file_ext(x=x)
        if (!frt %in% c("fa", "fasta")){
            stop("Invalid file format")
        }
        seqSSs <- readDNAStringSet(filepath=paste(workDir, x, sep="/"),
                                   format=frt)
    }
    ### Quality control of splice site sequences.
    if (!type %in% c(3, 5)){
        stop("Invalid type of splice sites. ",
             "The entered type of splice sites must be numeric value 3 or 5.")
    }
    if (((type == 3) & (unique(x=width(x=seqSSs))) != 23)){
        stop("Incorrect length of splice sites. ",
             "Acceptor splice site sequences must be 23 in length.")
    }
    if (((type == 5) & (unique(x=width(x=seqSSs))) != 9)){
        stop("Incorrect length of splice sites. ",
             "Donor splice site sequences must be 9 in length.")
    }
    if (!all(uniqueLetters(x=seqSSs) %in%
                                   c("a", "c", "g", "t", "G", "C", "T", "A"))){
        warning("Some sequences showed characters apart from A, C, G and T. ",
                "Scores can't be calculated for the incorrect records.")
    }
    ### Creation a new text file with splice site sequences.
    if (file.exists("seqSSs")){
        file.remove("seqSSs")
    }
    cat(as.character(x=seqSSs), file="seqSSs", sep="\n", append=TRUE)
    ### Calculation of the MaxEntScan scores.
    if (type == 3){
        cmd <- paste("score3.pl", "seqSSs")
    }
    if (type == 5){
        cmd <- paste("score5.pl", "seqSSs")
    }
    scores <- system2(command="perl", args=cmd, stdout=TRUE)
    if (length(x=scores) > 0){
        scores <- substr(x=scores,
                         start=regexpr(pattern="\t", text=scores)[[1]] + 1,
                         stop=nchar(x=scores))
    }
    file.remove("seqSSs")
    ### Returning the final object.
    scores <- cbind(names(x=seqSSs),
                    as.vector(x=as.character(x=seqSSs)),
                    scores)
    colnames(x=scores) <- c("ss_id", "ss_seq", "max.ent_score")
    scores <- data.frame(scores)
    scores[, 3] <- as.numeric(scores[, 3])
    scores <- scores[!duplicated(x=scores), ]
    return(scores)
}
