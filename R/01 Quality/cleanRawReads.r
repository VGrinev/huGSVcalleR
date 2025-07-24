#' Calculate quality scores for short reads
#' @description Calculation of the quality scores for short reads.
#' @param fastqDir a character string specifying the name of the directory
#'     containing the FASTQ file(s). NULL by default that mean the current
#'     working directory.
#' @param fastq a character string specifying the name of the FASTQ file
#'     containing short reads. Allowed formats are "fastq" or "fastq.gz".
#' @param n integer value giving the number of reads loaded into the computer's
#'     RAM. The default value is NULL, at which all reads are loaded into RAM.
#' @param adapters a character string specifying the name of the tab-delimited
#'     TXT file containing adapter sequences. The default value is NULL, at
#'     which the search for adaptor sequences in short reads is not carried out.
#'     If so, this file must contain the following three fields:
#'     i) adapter_id         - adapter ID;
#'     ii) adapter_name      - adapter name;
#'     iii) adapter_sequence - adapter sequence.
#' @param contaminants a character string specifying the name of the
#'     tab-delimited TXT file containing foreign sequences. The default value
#'     is NULL, at which the search for foreign sequences in short reads is not
#'     carried out. If so, this file must contain the following three fields:
#'     i) foreign_id         - foreign sequence ID;
#'     ii) foreign_name      - foreign sequence name;
#'     iii) foreign_sequence - foreign sequence itself.
#' @param workDir character string specifying the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return object of class list containing all main quality scores for short
#'     reads.
#' @author Vasily V. Grinev
#' @examples
#' qa <- assessQRawReads(fastqDir="Files_FASTQ",
#'                       fastq="2-Galkina-A_S2_L001_R1_001.fastq.gz",
#'                       n=1e5,
#'                       adapters=NULL,
#'                       contaminants=NULL,
#'                       workDir="D:/Vasily Grinev")
#' @export
#' Last updated: July 24, 2025.

fastqDir="Files_FASTQ"
fastq1="example_seq.read1.fastq"
fastq2="example_seq.read2.fastq"
tr_score=20
tr_start=10
tr_end=30
phred="phred_scores.txt"
k=3
halfwidth=5
successive=TRUE
readl=100
readq=20
N=3
dustScore=875
workDir="D:/Vasily Grinev"

cleanRawReads <- function(fastqDir=NULL,
                          fastq1,
                          fastq2=NULL,
                          tr_score=NULL,
                          tr_start=10,
                          tr_end=30,
                          phred=NULL,
                          k=3,
                          halfwidth=NULL,
                          successive=TRUE,
                          readl=0,
                          readq=0,
                          N=0,
                          dustScore=NULL,
                          workDir=NULL){
    ### Loading the required package.
    #   This code was successfully tested with the package ShortRead v.1.67.0.
    suppressMessages(expr=library(package=ShortRead))
    ### Loading the required auxiliary functions.
    source(file=paste(workDir, "cutLeftAdapter.r", sep="/"))
    source(file=paste(workDir, "cutRigthAdapter.r", sep="/"))
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
    frt1 <- tools::file_ext(x=fastq1)
    if (!frt1 %in% c("fastq", "gz")){
        stop("Invalid format of FASTQ file 1")
    }
    if (!is.null(x=fastq2)){
        frt2 <- tools::file_ext(x=fastq2)
        if (!frt2 %in% c("fastq", "gz")){
            stop("Invalid format of FASTQ file 2")
        }
    }
    ### Loading of the FASTQ file 1.
    FQ1 <- readFastq(dirPath=path, pattern=fastq1, withIds=TRUE)
    ### Raw reads processing.
    ##  Step 1, trimming of the adapter sequences.
    #   This step is under development.
    ##  Step 2, trimming of the low quality ends.
    #   Fixed trimming.
    if (is.null(x=tr_score)){
        FQ1 <- narrow(x=FQ1,
                      start=tr_start + 1,
                      end=max(x=width(x=FQ1)) - tr_end)
    }else{
    #   Dynamic trimming.
        PhSs <- read.table(file=paste(workDir, phred, sep="/"),
                           sep="\t",
                           header=TRUE,
                           quote="",
                           as.is=TRUE,
                           comment.char="z")
        FQ1 <- trimEnds(object=FQ1,
                        a=PhSs[PhSs$score == tr_score, ]$ASCII,
                        left=TRUE, right=FALSE,
                        relation="<=")
        if (!is.null(x=halfwidth)){
            FQ1 <- trimTailw(object=FQ1,
                             k=k,
                             a=PhSs[PhSs$score == tr_score, ]$ASCII,
                             halfwidth=halfwidth)
        }else{
            FQ1 <- trimTails(object=FQ1,
                             k=k,
                             a=PhSs[PhSs$score == tr_score, ]$ASCII,
                             successive=successive)
        }
    }
    ##  Step 3, removing of the too short reads.
    FQ1 <- FQ1[width(x=FQ1) >= readl]
    ##  Step 4, removing of the low quality reads.
    FQ1 <- FQ1[rowMeans(x=as(quality(x=FQ1), "matrix"), na.rm=TRUE) >= readq]
    ##  Step 5, removing of the reads with ambiguous nucleotides.
    FQ1 <- FQ1[polynFilter(threshold=as.integer(x=N), nuc="other")(FQ1)]
    ##  Step 6, removing of the duplicated reads.
    FQ1 <- FQ1[occurrenceFilter(min=1L, max=1L,
                                withSread=TRUE,
                                duplicates="sample")(FQ1)]
    FQ2 <- FQ1[dustyFilter(threshold=dustScore, batchSize=NA)(FQ1)]
    return()
}
