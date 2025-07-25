#' Calculate pileup statistics for a BAM file
#' @description This is a high-level function that performs a multi-threaded
#'     calculation of the nucleotide-wise read coverage based on the BAM file
#'     of interest.
#' @param bamDir character string giving the name (or path to and name) of
#'     directory with BAM file(s). NULL by default, which means the current
#'     working directory will be used instead of the specified BAM directory.
#' @param bamFile character string giving the name of input BAM file.
#' @param gr character string indicating the method for specifying the genomic
#'     location(-s) for which pileup statistics will be collected. Default
#'     value is "SeqInfo". In this case, the locations will be extracted
#'     directly from the BAM file and pileup statistics will be collected for
#'     the entire genome. Alternatively, the user can specify the genomic
#'     locations of interest by applying to the argument a name of the
#'     tab-delimited TXT file containing genomic coordinates. This file mus
#'     contains the following four fields:
#'     i) seqnames - the name of chromosome or scaffold with prefix "chr";
#'     ii) start   - start coordinate of genomic interval of interest;
#'     iii) end    - end coordinate of genomic interval of interest;
#'     iv) strand  - (optionally) strand information about genomic interval
#'                   of interest.
#'     TXT file with genomic location(-s) must be in working directory.
#' @param canChr logical argument specifying that the records from non-canonical
#'     chromosomes/scaffolds should be removed. Default value is TRUE.
#' @param depth integer value specifying the maximum number of overlapping
#'     alignments considered for each position in the pileup.
#'     Default value is 1e6.
#' @param baseq integer value specifying the minimal base quality for each
#'     nucleotide in an alignment to be included in pileup. Default value is 20.
#' @param mapq integer value specifying the minimal mapping quality for an
#'     alignment to be included in pileup. Default value is 20.
#' @param coverage integer value specifying the minimal coverage threshold for
#'     filtering of successful genomic positions.
#' @param genome character string giving the name of the BSgenome data package
#'     to be used for reference sequence extraction. Default value is
#'     package BSgenome.Hsapiens.UCSC.hg38.
#' @param tmpDir character string giving the name (or path to and name) of
#'     directory for storing of temporary file(-s) with pileup statistics.
#'     NULL by default, which means the current working directory will be used
#'     instead of the specified directory.
#' @param postfix character string that will be added to the CSV output file
#'     name as postfix. NULL by default. In this case, the postfix
#'     "PileupResults" will be used.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @return a CSV file containing pileup statistics will be written to the
#'     directory where input BAM file is stored.
#' @details 
#'     - The function processes the BAM file in parallel chunks, each
#'     corresponding to a genomic range.
#'     - It computes the number of occurrences of each nucleotide (A, C, G, T)
#'     at each position within the specified range.
#'     - The function includes a coverage threshold to filter positions with
#'     low coverage.
#'     - It uses reference genome to retrieve the reference nucleotide for each
#'     position in the pileup data.
#'     - The final pileup statistics is saved in the specified CSV file, which
#'     contains information about the position, reference nucleotide, and
#'     counts for A, C, G and T.
#' @author Vasily V. Grinev, Dzianis D. Sarnatski, Mikalai M. Yatskou
#' @examples
#' res <- collectBasesPileup(bamDir="Files_BAM",
#'                           bamFile="SRR21721515.bam",
#'                           gr="genomic_intervals.txt",
#'                           canChr=TRUE,
#'                           depth=1e6,
#'                           baseq=20,
#'                           mapq=20,
#'                           coverage=10,
#'                           genome="BSgenome.Hsapiens.UCSC.hg38",
#'                           tmpDir=NULL,
#'                           postfix="PileupResults",
#'                           workDir="D:/Vasily Grinev")
#' @export
#' Last updated: July 24, 2025.

collectBasesPileup <- function(bamDir=NULL,
                               bamFile,
                               gr="SeqInfo",
                               canChr=TRUE,
                               depth=1e6, baseq=20, mapq=20,
                               coverage=NULL,
                               genome="BSgenome.Hsapiens.UCSC.hg38",
                               tmpDir=NULL,
                               postfix=NULL,
                               workDir=NULL){
    ### Loading the required package.
    #   This code was successfully tested with the packages BSgenome v.1.72,
    #   BSgenome.Hsapiens.UCSC.hg38 v.1.4.5, data.table v.1.15.4,
    #   doParallel v.1.0.17, Rsamtools v.2.20.0 and tidyverse v.2.0.0.
    suppressMessages(expr=library(data.table))
    suppressMessages(expr=library(doParallel))
    suppressMessages(expr=library(genome, character.only=TRUE))
    suppressMessages(expr=library(Rsamtools))
    suppressMessages(expr=library(tidyverse))
    ### Setting the working directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Setting the temporary directory.
    if (is.null(x=tmpDir)){
        tmpDir <- workDir
    }
    ### Full path to the BAM file.
    if (is.null(x=bamDir)){
        path <- workDir
    }else{
        path <- paste(workDir, bamDir, sep="/")
    }
    ### Setting the full path to the BAM file and reference to this file.
    bam <- paste(path, bamFile, sep="/")
    bam_link <- BamFile(file=bam)
    ### Setting the genomic location(-s) of interest.
    if (gr == "SeqInfo"){
        GRs <- as.data.frame(cbind(seqinfo(x=bam_link)@seqnames, 1,
                                   seqinfo(x=bam_link)@seqlengths))
        colnames(x=GRs) <- c("seqnames", "start", "end")
        GRs$start <- as.numeric(x=GRs$start)
        GRs$end <- as.numeric(x=GRs$end)
        if (isTRUE(x=canChr)){
            if (grepl(pattern="chr",
                      x=seqinfo(x=bam_link)@seqnames[1]) == "TRUE"){
                GRs <- GRs[nchar(x=GRs$seqnames) <= 5, ]
            }else{
                GRs <- GRs[nchar(x=GRs$seqnames) == 1, ]
            }
        }
        rownames(x=GRs) <- NULL
        GRs <- makeGRangesFromDataFrame(df=GRs)
    }else{
        GRs <- read.table(file=paste(workDir, gr, sep="/"),
                          sep="\t",
                          header=TRUE,
                          quote="\"",
                          as.is=TRUE)
        if (isTRUE(x=canChr)){
            if (grepl(pattern="chr",
                      x=seqinfo(x=bam_link)@seqnames[1]) == "TRUE"){
                GRs <- GRs[nchar(x=GRs$seqnames) <= 5, ]
            }else{
                GRs <- GRs[nchar(x=GRs$seqnames) == 1, ]
            }
        }
        if (grepl(pattern="chr", x=seqinfo(x=bam_link)@seqnames[1]) == "TRUE"){
            GRs <- makeGRangesFromDataFrame(df=GRs)
        }else{
            GRs$seqnames <- gsub(pattern="chr", replacement="", x=GRs$seqnames)
            GRs <- makeGRangesFromDataFrame(df=GRs)
        }
    }
    ### Setting the multithreads process.
    numCores <- detectCores() - 1
    cl <- makeCluster(numCores)
    registerDoParallel(cl)
    ### Collection of pileup statistics.
    ##  Processing of genomic location(-s) as chunk(-s) based on memory limits.
    foreach(i=seq_along(GRs), .combine=rbind, .packages=c("tidyverse",
                                                          "data.table",
                                                          "GenomicRanges",
                                                          "Rsamtools",
                                                          genome)) %dopar% {
    #   Get pileup statistics for the current chunk.
    chunk <- GRs[i]
    pileup_data <- pileup(file=bam_link,
                          scanBamParam=ScanBamParam(which=chunk),
                          pileupParam=PileupParam(max_depth=depth,
                                                  min_base_quality=baseq,
                                                  min_mapq=mapq,
                                                  distinguish_strands=FALSE,
                                                  include_deletions=FALSE,
                                                  include_insertions=FALSE))
    #   Transformation of the primary pileup object to a six-fields table
    #   object of class data.frame.
    pileup_data$position <- paste(pileup_data$seqnames,
                                  pileup_data$pos,
                                  sep=":")
    #   Cumulative statistics.
    temp_stat <- pileup_data %>%
                    mutate(reference=0, A=0, C=0, G=0, T=0) %>%
                       mutate(A=case_when(nucleotide == "A" ~ count, TRUE ~ A),
                              C=case_when(nucleotide == "C" ~ count, TRUE ~ C),
                              G=case_when(nucleotide == "G" ~ count, TRUE ~ G),
                              T=case_when(nucleotide == "T" ~ count, TRUE ~ T))
    temp_stat <- as.data.table(temp_stat)[, .(seqnames=seqnames,
                                              pos=pos,
                                              reference=sum(x=reference),
                                              A=sum(x=A),
                                              C=sum(x=C),
                                              G=sum(x=G),
                                              T=sum(x=T)),
                                          by=position][, 
                             .(position, seqnames, pos, reference, A, C, G, T)]
    #   Filtration against too low coverage.
        if (!is.null(x=coverage)){
            temp_stat <- temp_stat %>%
                            filter(rowSums(x=across(c(A, C, G, T))) >= coverage)
        }
    #   Reference nucleotides.
    ref <- temp_stat %>%
              select(seqnames, pos) %>%
              mutate(start=pos) %>%
              mutate(end=pos) %>%
              select(-pos)
    ref <- makeGRangesFromDataFrame(df=ref)
    temp_stat <- temp_stat %>%
                    mutate(reference=getSeq(x=get(x=genome),
                                            names=ref,
                                            as.character=TRUE)) %>%
                    select(-pos, -seqnames)
    #   Writing temporary file with pileup statistics.
        if (nrow(x=temp_stat) > 0){
            fwrite(x=temp_stat,
                   file=paste0(tmpDir,
                               "/",
                               "temp_stat_",
                               gsub(pattern=":",
                                    replacement="_",
                                    x=temp_stat[1, "position"]),
                               ".csv"))
        }
    }
    stopCluster(cl)
    ##  Consolidation of pileup statistics.
    #   Listing of temporary file(-s).
    setwd(dir=tmpDir)
    files <- list.files(path=tmpDir, pattern="temp_stat_.*\\.csv$")
    if (length(x=files) == 0){
        files <- Sys.glob(paths="temp_stat*")
    }
    #   Extracting of chromosome and position data.
    extract_chr_pos <- function(filename){
        matches <- str_match(string=filename,
                             pattern="temp_stat_chr(\\d+|[XY])_(\\d+)\\.csv")
        chr <- matches[, 2]
        pos <- as.numeric(matches[, 3])
        chr <- ifelse(chr == "X", 23, ifelse(chr == "Y", 24, as.numeric(x=chr)))
        list(chr=chr, pos=pos)
    }
    file_info <- map_dfr(files, function(file){
        info <- extract_chr_pos(file)
        tibble(file=file, chr=info$chr, pos=info$pos)
    })
    #   Sorting files by chromosomes and positions.
    files <- file_info %>%
                 arrange(chr, pos) %>%
    pull(file)
    #   Name of the final output file.
    if (!is.null(x=postfix)){
        csv_file <- paste(gsub(pattern=".bam", replacement="", x=bam),
                          postfix,
                          "csv",
                          sep=".")
    }else{
        csv_file <- paste(gsub(pattern=".bam", replacement="", x=bam),
                          "PileupResults",
                          "csv",
                          sep=".")
    }
    #   Removing of old output file (if so).
    invisible(if (file.exists(csv_file)){
                  file.remove(csv_file)
              }
    )
    #   Collecting data into a single file.
    for (file in files){
        temp_df <- fread(file=file)
        if (!file.exists(csv_file)){
            fwrite(x=temp_df, file=csv_file)
        }else{
            fwrite(x=temp_df, file=csv_file, append=TRUE)
        }
        file.remove(file)
    }
    temp_df <- fread(file=csv_file)
    temp_df <- temp_df[!duplicated(x=temp_df), ]
    fwrite(x=temp_df, file=csv_file)
}
