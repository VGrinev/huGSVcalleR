#' Read records from FASTQ-formatted file
#' @description Reading the records from FASTQ-formatted file.
#' @param fastqDir character string specifying the name of the directory
#'     containing the FASTQ file(s). NULL by default that mean the current
#'     working directory.
#' @param fastq character string specifying the name of the FASTQ file
#'     containing short reads. Allowed formats are "fastq" or "fastq.gz".
#' @param n integer value giving the number of reads loaded into the computer's
#'     RAM. The default value is NULL, at which all reads are loaded into RAM.
#' @param workDir character string specifying the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return object of class list containing DNAStringSet object of read
#'     sequences and PhredQuality object of read quality scores.
#' @author Vasily V. Grinev
#' @examples
#' fq <- readFASTQ(fastqDir="Files_FASTQ",
#'                 fastq="example_seq.read1.fastq.gz",
#'                 n=1e5,
#'                 workDir="D:/Vasily Grinev")
#' @export
#' @importFrom Biostrings DNAStringSet PhredQuality
#' @importFrom tools file_ext
readFASTQ <- function(fastqDir=NULL,
                      fastq,
                      n=NULL,
                      workDir=NULL){
    ### Setting the working directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Full path to the FASTQ file(-s).
    if (is.null(x=fastqDir)){
        path <- paste(workDir, fastqDir, sep="")
    }else{
        path <- paste(workDir, fastqDir, sep="/")
    }
    ### Retrieving and validation of the FASTQ file(-s) extension.
    frt <- tools::file_ext(x=fastq)
    if (!frt %in% c("fastq", "gz")){
        stop("Invalid file format")
    }
    ### Loading the FASTQ data.
    if (is.null(x=n)){
        reads <- readLines(con=paste(path, fastq, sep="/"))
        FQs <- list(reads=DNAStringSet(x=reads[seq(from=2,
                                                   to=length(x=reads),
                                                   by=4)]),
                    QScores=PhredQuality(x=reads[seq(from=4,
                                                     to=length(x=reads),
                                                     by=4)]))
        names(x=FQs$reads) <- reads[seq(from=1, to=length(x=reads), by=4)]
        names(x=FQs$QScores) <- reads[seq(from=1, to=length(x=reads), by=4)]
    }else{
        reads <- readLines(con=paste(path, fastq, sep="/"), n=n * 4)
        FQs <- list(reads=DNAStringSet(x=reads[seq(from=2, to=n * 4, by=4)]),
                    QScores=PhredQuality(x=reads[seq(from=4, to=n * 4, by=4)]))
        names(x=FQs$reads) <- reads[seq(from=1, to=n * 4, by=4)]
        names(x=FQs$QScores) <- reads[seq(from=1, to=n * 4, by=4)]
    }
    rm(reads)
    ### Returning the final object.
    return(FQs)
}
