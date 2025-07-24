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
#' Last updated: July 12, 2025.

assessQRawReads <- function(fastqDir=NULL,
                            fastq,
                            n=NULL,
                            adapters=NULL,
                            contaminants=NULL,
                            workDir=NULL){
    ### Loading the required packages.
    #   This code was successfully tested with the packages Biostrings v.2.72.1,
    #   data.table v.1.15.4, matrixStats v.1.3.0 and ShortRead v.1.67.0.
    suppressMessages(expr=library(package=Biostrings))
    suppressMessages(expr=library(package=data.table))
    suppressMessages(expr=library(package=matrixStats))
    suppressMessages(expr=library(package=ShortRead))
    ### Setting the working directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Loading the required auxiliary functions.
    source(file=paste(workDir, "readFASTQ.r", sep="/"))
    source(file=paste(workDir, "findAdapters.r", sep="/"))
    source(file=paste(workDir, "findForeignSeqs.r", sep="/"))
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
    ### Collection of the quality metrics.
    ##  Sample name.
    sample_name <- gsub(pattern=".fastq.+", replacement="", x=fastq)
    ##  File name.
    file_name <- fastq
    ##  Total number of reads.
    TNBs <- countFastq(dirPath=path, pattern=fastq)
    TNRs <- TNBs[[1]]
    ##  Total number of bases.
    TNBs <- TNBs[[2]]
    ##  Distribution the length of reads.
    reads <- readFASTQ(fastqDir=fastqDir, fastq=fastq, n=n, workDir=workDir)
    RLs <- as.numeric(x=summary(object=width(x=reads$reads)))
    ##  Per entire reads quality.
    PERQs <- mean(x=as(object=reads$QScores, Class="IntegerList"))
    names(x=PERQs) <- NULL
    PERQs <- list(means=PERQs,
                  density=do.call(what=cbind, args=density(x=PERQs)[1:2]),
                  metrics=quantile(x=PERQs,
                                   probs=c(0.1, 0.25, 0.5, 0.75, 0.9)))
    ##  Per cycle quality of reads.
    if (is.null(x=n)){
        QA <- qa(dirPath=path, pattern=fastq)
    }else{
        QA <- qa(dirPath=path, pattern=fastq, n=n)
    }
    PCQRs <- do.call(what=rbind,
                     args=lapply(X=split(x=QA[["perCycle"]]$quality,
                                         f=QA[["perCycle"]]$quality$Cycle),
                                 FUN=function(y){rep(x=y[, 3], times=y[, 4])}))
    m <- list()
    for (i in seq(from=1, to=nrow(x=PCQRs), by=5)){
        if (i + 4 < nrow(x=PCQRs)){
            m[[i]] <- c(paste(i, i:i + 4, sep="-"),
                        quantile(x=PCQRs[i:i + 4, ],
                                 probs=c(0.1, 0.25, 0.5, 0.75, 0.9)))
        }else{
            m[[i]] <- c(paste(i, nrow(x=PCQRs), sep="-"),
                        quantile(x=PCQRs[i:nrow(x=PCQRs), ],
                                 probs=c(0.1, 0.25, 0.5, 0.75, 0.9)))

        }
    }
    PCQRs <- data.frame(do.call(args=m, what=rbind))
    colnames(x=PCQRs) <- c("cycles", "minimum", "Q1", "median", "Q3", "maximum")
    PCQRs[, -1] <- apply(X=PCQRs[, -1], MARGIN=2, FUN=as.numeric)
    ##  Bases composition of the entire set of reads.
    BCs <- letterFrequency(x=reads$reads,
                           letters=c("A", "C", "G", "T", "N"))
    BCs <- BCs/rowSums(x=BCs)
    metrics=rbind(quantile(x=BCs[, 1], probs=c(0.1, 0.25, 0.5, 0.75, 0.9)),
                  quantile(x=BCs[, 2], probs=c(0.1, 0.25, 0.5, 0.75, 0.9)),
                  quantile(x=BCs[, 3], probs=c(0.1, 0.25, 0.5, 0.75, 0.9)),
                  quantile(x=BCs[, 4], probs=c(0.1, 0.25, 0.5, 0.75, 0.9)),
                  quantile(x=BCs[, 5], probs=c(0.1, 0.25, 0.5, 0.75, 0.9)))
    metrics <- data.frame(cbind(colnames(x=BCs), metrics))
    colnames(x=metrics) <- c("nucleotide",
                             "minimum", "Q1", "median", "Q3", "maximum")
    metrics[, -1] <- apply(X=metrics[, -1], MARGIN=2, FUN=as.numeric)
    BCs <- list(frequency=BCs, metrics=metrics)
    ##  Per cycle content of bases.
    PCCBs <- QA[["perCycle"]]$baseCall[, -4]
    PCCBs$Base <- as.character(x=PCCBs$Base)
    PCCBs <- do.call(what=rbind,
                     args=lapply(X=split(x=PCCBs,
                                         f=PCCBs$Cycle),
                                 FUN=function(y){rep(x=y[, 2], times=y[, 3])}))
    m <- data.frame(cycles="", A=0, C=0, G=0, T=0, N=0)
    for (i in seq(from=1, to=nrow(x=PCCBs), by=5)){
        if (i + 4 < nrow(x=PCCBs)){
            bases <- c(paste(i, i:i + 4, sep="-"),
                       table(PCCBs[i:i + 4, ]))
            m[i, 1] <- bases[1]
            m[i, -1] <- bases[match(colnames(x=m), names(x=bases))][-1]
        }else{
            bases <- c(paste(i, nrow(x=PCCBs), sep="-"),
                       table(PCCBs[i:nrow(x=PCCBs), ]))
            m[i, 1] <- bases[1]
            m[i, -1] <- bases[match(colnames(x=m), names(x=bases))][-1]
        }
    }
    PCCBs <- m[!is.na(x=m$cycles), ]
    PCCBs[is.na(x=PCCBs$N), ]$N <- 0
    PCCBs[, -1] <- apply(X=PCCBs[, -1], MARGIN=2, FUN=as.numeric)
    rownames(x=PCCBs) <- NULL
    PCCBs[, -1] <- PCCBs[, -1]/rowSums(x=PCCBs[, -1])
    ##  GC composition of the entire set of reads.
    GCs <- letterFrequency(x=reads$reads, letters=c("G", "C"))
    GCs <- rowSums(x=GCs)/width(x=reads$reads)
    GCs <- list(frequency=GCs,
                density=do.call(what=cbind, args=density(x=GCs)[1:2]),
                metrics=quantile(x=GCs, probs=c(0.1, 0.25, 0.5, 0.75, 0.9)))
    ##  Summarization the low-complexity of reads (using dusty scores).
    DSs <- dustyScore(x=reads$reads)
    DSs <- list(scores=DSs,
                density=do.call(what=cbind, args=density(x=log10(x=DSs))[1:2]),
                metrics=quantile(x=DSs, probs=c(0.1, 0.25, 0.5, 0.75, 0.9)))
    ##  Frequency of reads.
    FRs <- QA[["sequenceDistribution"]][, 1:2]
    colnames(x=FRs) <- c("occurrence", "reads")
    FRs$frequency <- FRs$reads/sum(x=FRs$reads)
    ##  Overrepresented reads.
    ORs <- ShortRead:::.freqSequences(QA, "read")[, 1:2]
    ##  Identification of irrelevant sequences.
    #   Identification of adapters.
    ASs <- NA
    if (!is.null(x=adapters)){
        ASs <- findAdapters(fastqDir=fastqDir,
                            fastq=fastq,
                            n=n,
                            adapters=adapters,
                            workDir=workDir)
    }
    #   Identification of foreign sequences.
    FSs <- NA
    if (!is.null(x=contaminants)){
        FSs <- findForeignSeqs(fastqDir=fastqDir,
                               fastq=fastq,
                               n=n,
                               contaminants=contaminants,
                               workDir=workDir)
    }
    ### Quality metrics.
    return(list(SN=sample_name,
                FN=file_name,
                TNRs=TNRs,
                TNBs=TNBs,
                RLs=RLs,
                PERQs=PERQs,
                PCQRs=PCQRs,
                BCs=BCs,
                PCCBs=PCCBs,
                GCs=GCs,
                DSs=DSs,
                FRs=FRs,
                ORs=ORs,
                ASs=ASs,
                FSs=FSs)
           )
}

