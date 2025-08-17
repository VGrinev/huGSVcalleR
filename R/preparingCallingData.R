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
#'     tab-delimited TXT file containing genomic coordinates. This file must
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
#' @importFrom data.table data.table fwrite as.data.table
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom Rsamtools BamFile ScanBamParam seqinfo pileup PileupParam
#' @importFrom GenomicRanges makeGRangesFromDataFrame GRanges
#' @importFrom BSgenome getSeq
#' @importFrom dplyr mutate case_when filter arrange pull select across
#' @importFrom stringr str_match
#' @importFrom purrr map_dfr
#' @importFrom tidyr %>%
#' @importFrom utils read.table
#' @importFrom methods is
#' @importFrom foreach foreach %dopar%
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


#' Filter a BAM file
#' @description Multiparameter BAM file filtering.
#' @param bamDir character string giving the name (or path to and name) of
#'     directory with BAM file(s). NULL by default, which means the current
#'     working directory will be used instead of the specified BAM directory.
#' @param bamFile character string giving the name of input BAM file.
#' @param gr character string specifying the name of the tab-delimited TXT file
#'     containing coordinates of genomic location(-s) of interest. The default
#'     value is NULL. If so, this file must contains the following four fields:
#'     i) seqnames - the name of chromosome or scaffold with prefix "chr";
#'     ii) start   - start coordinate of genomic interval of interest;
#'     iii) end    - end coordinate of genomic interval of interest;
#'     iv) strand  - (optionally) strand information about genomic interval
#'                   of interest.
#'     TXT file with genomic location(-s) must be in working directory.
#' @param flag logical vector used to filter reads based on their flag entry.
#'     This is most easily created with the helper function scanBamFlag() from
#'     package Rsamtools. The function scanBamFlag() recognizes the following
#'     flags: isPaired, isProperPair, isUnmappedQuery, hasUnmappedMate,
#'     isMinusStrand, isMateMinusStrand, isFirstMateRead, isSecondMateRead,
#'     isNotPrimaryRead, isSecondaryAlignment, isNotPassingQualityControls,
#'     isDuplicate and isSupplementaryAlignment. The chosen flag must have
#'     value TRUE to be involved in filtration process and only such flag
#'     should be submitted to the function. Filtering by flag is not performed
#'     by default.
#' @param tlen integer value limiting the maximum template length. NULL by
#'     default that mean no any limitations.
#' @param mapq integer value limiting the minimum mapping quality. NULL by
#'     default that mean no any limitations.
#' @param index logical, FALSE by default. If TRUE, it allows to sort and index
#'     a new filtered BAM file.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @return filtered BAM file (optionally, sorted and indexed if index=TRUE)
#'     with postfix "filtered".
#' @author Vasily V. Grinev
#' @examples
#' filterBamFile(bamDir="Files_BAM",
#'               bamFile="SRR21721515.bam",
#'               gr="genomic_intervals.txt",
#'               flag=scanBamFlag(isPaired=TRUE, isProperPair=TRUE),
#'               tlen=300,
#'               mapq=20,
#'               index=TRUE,
#'               workDir="D:/Vasily Grinev")
#' @export
#' @importFrom Rsamtools scanBamFlag ScanBamParam filterBam sortBam indexBam
#' @importFrom GenomicRanges makeGRangesListFromDataFrame
#' @importFrom IRanges IRangesList
#' @importFrom parallel detectCores
#' @importFrom utils read.table
filterBamFile <- function(bamDir=NULL,
                          bamFile,
                          gr=NULL,
                          flag=scanBamFlag(),
                          tlen=NULL,
                          mapq=NULL,
                          index=FALSE,
                          workDir=NULL){
  ### Setting the working directory.
  if (is.null(x=workDir)){
    workDir <- getwd()
  }
  ### Full path to the BAM file.
  if (is.null(x=bamDir)){
    path <- workDir
  }else{
    path <- paste(workDir, bamDir, sep="/")
  }
  ### Development an object of class IRangesList to filter the BAM file of
  #   interest in given genomic interval(-s) (if so).
  if (!is.null(x=gr)){
    GRs <- read.table(file=paste(workDir, gr, sep="/"),
                      sep="\t",
                      header=TRUE,
                      quote="\"",
                      as.is=TRUE)
    GRs$split <- GRs$seqnames
    GRs <- makeGRangesListFromDataFrame(df=GRs, split.field="split")
    GRs <- IRangesList(GRs)
  }
  ### Filtering the BAM file of interest.
  param <- ScanBamParam(flag=flag)
  if (!is.null(x=gr)){
    param@which <- GRs
  }
  if (!is.null(x=tlen)){
    param@what <- "isize"
  }
  if (!is.null(x=mapq)){
    param@mapqFilter <- as.integer(x=mapq)
  }
  if (!is.null(x=tlen)){
    filter <- FilterRules(list(function(x){x$isize <= tlen}))
  }else{
    filter <- FilterRules()
  }
  dest <- gsub(pattern="bam",
               replacement="filtered.bam",
               x=paste(path, bamFile, sep="/"))
  f <- filterBam(file=paste(path, bamFile, sep="/"),
                 index=paste(path, bamFile, sep="/"),
                 filter=filter,
                 param=param,
                 destination=dest,
                 indexDestination=FALSE)
  if (isTRUE(x=index)){
    s <- sortBam(file=f,
                 destination=gsub(pattern=".bam", replacement="", x=f),
                 nThreads=detectCores(logical=TRUE) - 1)
    indexBam(files=s)
  }
}


#' Group short reads in BAM file
#' @description This is a wrapper function that uses the low-level function
#'     AddOrReplaceReadGroup() from package GATK to add read group tags to BAM
#'     file of interest. This stage could be not nessesery, if BAM file already
#'     contains read group tags.
#' @param bamDir character string giving the name (or path to and name) of
#'     directory with BAM file(s). NULL by default, which means the current
#'     working directory will be used instead of the specified BAM directory.
#' @param bamFile character string giving the name of input BAM file.
#' @param postfix character string that will be added to the output file name
#'     as postfix. "RG" by default.
#' @param SE if TRUE (by default), standart error determined by function
#'     AddOrReplaceReadGroup() will be writen in file STDERR.RG.TXT.
#' @param RGID RGID value for appropriate BAM field.
#' @param RGLB RGLB value for appropriate BAM field.
#' @param RGPL RGPL value for appropriate BAM field.
#' @param RGPU RGPU value for appropriate BAM field.
#' @param RGSM RGSM value for appropriate BAM field.
#' @param gatk_path character string giving the path to internal or external
#'     (system) GATK executive script. This path can be obtained by running
#'     Find.Gatk() or settled manually. Default value is "gatk", i.e. function
#'     will uses GATK systemly from $PATH.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @return a new BAM file with additional information about read groups and
#'     (optionally) auxillary file with standart error log.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @examples
#' ...
#' @export
groupReadsGATK <- function(bamDir=NULL,
                           bamFile,
                           postfix="RG",
                           SE=TRUE,
                           RGID="Bam.Recal.RG",
                           RGLB="default.LIBRARY.ID",
                           RGPL="default.ILLUMINA",
                           RGPU="default.PLATFORM",
                           RGSM="default.SAMPLE",
                           gatk_path="gatk",
                           workDir=NULL){
  ### Setting the working directory.
  if (is.null(x=workDir)){
    workDir <- getwd()
  }
  ### Full path to the BAM file.
  if (is.null(x=bamDir)){
    path <- workDir
  }else{
    path <- paste(workDir, bamDir, sep="/")
  }
  ### Naming of output file.
  output <- paste(gsub(pattern=".bam",
                       replacement="",
                       x=bamFile),
                  postfix,
                  "bam",
                  sep=".")
  ### Settings for saving of standart error.
  se <- ""
  if (SE == TRUE){
    se <- paste(gsub(pattern=".bam", replacement="", x=bamFile),
                "STDERR.RG.TXT",
                sep=".")
  }
  ### Running of command.
  system2(command=gatk_path,
          args=c("AddOrReplaceReadGroups",
                 paste("-I", paste(path, bamFile, sep="/"), sep=" "),
                 paste("-O", paste(path, output, sep="/"), sep=" "),
                 paste("--RGID", RGID, "--RGLB", RGLB, "--RGPL", RGPL,
                       "--RGPU", RGPU, "--RGSM", RGSM, sep=" ")),
          if (se == ""){
            stderr=se
          }else{
            stderr=paste(path, se, sep="/")
          })
}


#' Index a reference genome
#' @description This is a wrapper function that uses the low-level function
#'     CreateSequenceDictionary() from package GATK to built .dict and/or .fai
#'     indices for reference genome of interest.
#' @param fastaDir character string giving the name (or path to and name) of
#'     folder for storing of reference genome.
#' @param faFile charater string giving the name of file containing reference
#'     genome in format *.FASTA or *.FA.
#' @param fai logical. If TRUE (by default), .fai index of FASTA file will be
#'     created in working diectory.
#' @param dict logical. If TRUE (by default), .dict index of FASTA file will be
#'     created in working diectory.
#' @param gatk_path character string giving the path to internal or external
#'     (system) GATK executive script. This path can be obtained by running
#'     Find.Gatk() or settled manualy. Default value is "gatk", i.e. function
#'     will uses GATK systemly from $PATH.
#' @param workDir character string giving the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return index files in working directory.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @examples
#' ...
#' @export
#' @importFrom Rsamtools indexFa
indexGenomeGATK <- function(fastaDir=NULL,
                            faFile,
                            fai=TRUE,
                            dict=TRUE,
                            gatk_path="gatk",
                            workDir=NULL){
  ### Setting the working directory.
  if (is.null(x=workDir)){
    workDir <- getwd()
  }
  ### Full path to the FASTA/FA file.
  if (is.null(x=fastaDir)){
    path <- workDir
  }else{
    path <- paste(workDir, fastaDir, sep="/")
  }
  ### Building the .fai index.
  if (fai == TRUE){
    indexFa(file=paste(path, faFile, sep="/"))
  }
  ### Building the .dict index.
  if (dict == TRUE){
    system2(command=gatk_path,
            args=c("CreateSequenceDictionary",
                   "-R",
                   paste(path, faFile, sep="/")))
  }
}


#' Index a VCF file with known variants
#' @description This is a wrapper function that uses the low-level function
#'     IndexFeatureFile() from package GATK to built a .tbi index for VCF file
#'     with known variants.
#' @param vcfDir character string giving the name (or path to and name) of
#'     directory with VCF file(s). NULL by default, which means that the
#'     current working directory will be used instead of the specified VCF
#'     directory.
#' @param vcfFile character string giving the name of input VCF file.
#' @param output_name character string giving the name of index. NULL by
#'     default, which means that the name of input VCF file will be uased.
#'     Important note: to be used correctly by the GATK, the VCF file and its
#'     index must have the same name.
#' @param gatk_path character string giving the path to internal or external
#'     (system) GATK executive script. This path can be obtained by running
#'     Find.Gatk() or settled manually. Default value is "gatk", i.e. function
#'     will uses GATK systemly from $PATH.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @return .tbi file in working directory.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @examples
#' ...
#' @export
indexVcfGATK <- function(vcfDir=NULL,
                         vcfFile=NULL,
                         output_name=NULL,
                         gatk_path="gatk",
                         workDir=NULL){
  ### Setting the working directory.
  if (is.null(x=workDir)){
    workDir <- getwd()
  }
  ### Full path to the VCF file.
  if (is.null(x=vcfDir)){
    path <- workDir
  }else{
    path <- paste(workDir, vcfDir, sep="/")
  }
  ### Running of command.
  if (is.null(x=output_name)){
    system2(command=gatk_path,
            args=c("IndexFeatureFile -I",
                   paste(path, vcfFile, sep="/")))
  }else{
    system2(command=gatk_path,
            args=c("IndexFeatureFile -I",
                   paste(path, vcfFile, sep="/"),
                   "-O", output_name))
  }
}


#' Detect and mark duplicated reads
#' @description This is a wrapper function that uses the low-level function
#'     MarkDuplicates() from package GATK to find out and mark duplicated reads
#'     in BAM file(-s) of interest.
#' @param bamDir character string giving the name (or path to and name) of
#'     directory with BAM file(s). NULL by default, which means the current
#'     working directory will be used instead of the specified BAM directory.
#' @param bamFile character string giving the name of input BAM file.
#' @param postfix character string that will be added to the output file name
#'     as postfix. "markedDuplicates" by default.
#' @param removeDupls a way to deal with detected duplicated reads. With "no"
#'     all duplicated reads will be kept in output BAM file, with "yes" all
#'     duplicated reads will be removed, and with "seq" only optical and other
#'     sequencing duplicates will be removed.
#' @param SE if TRUE (by default), standart error determined by function
#'     MarkDuplicates() will be writen in file STDERR.MARKDUP.TXT.
#' @param gatk_path character string giving the path to internal or external
#'     (system) GATK executive script. This path can be obtained by running
#'     Find.Gatk() or settled manually. Default value is "gatk", i.e. function
#'     will uses GATK systemly from $PATH.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @return a new BAM file with processed duplicated reads, a file with metrics
#'     of input BAM file and (optionally) a file with standart error log.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @examples
#' ...
#' @export
markDuplicatesGATK <- function(bamDir=NULL,
                               bamFile,
                               postfix="markedDuplicates",
                               removeDupls="no",
                               SE=TRUE,
                               gatk_path="gatk",
                               workDir=NULL){
  ### Setting the working directory.
  if (is.null(x=workDir)){
    workDir <- getwd()
  }
  ### Full path to the BAM file.
  if (is.null(x=bamDir)){
    path <- workDir
  }else{
    path <- paste(workDir, bamDir, sep="/")
  }
  ### Settings of run.
  input <- paste("-I",
                 paste(path, bamFile, sep="/"),
                 sep=" ")
  output <- paste("-O",
                  paste(path,
                        paste(gsub(pattern=".bam",
                                   replacement="",
                                   x=bamFile),
                              postfix,
                              "bam",
                              sep="."),
                        sep="/"),
                  sep=" ")
  metrics <- paste("-M",
                   paste(path,
                         paste(gsub(pattern=".bam",
                                    replacement="",
                                    x=bamFile),
                               postfix,
                               "metrics.txt",
                               sep="."),
                         sep="/"),
                   sep=" ")
  ### Settings for removing duplicated reads.
  if (removeDupls == "no"){
    remove.d <- paste("--REMOVE_DUPLICATES", "FALSE")
  }else{
    if (removeDupls == "yes"){
      remove.d <- paste("--REMOVE_DUPLICATES", "TRUE")
    }else{
      if (removeDupls == "seq"){
        remove.d <- paste("--REMOVE_DUPLICATES", "FALSE",
                          "--REMOVE_SEQUENCING_DUPLICATES", "TRUE")
      }
    }
  }
  ### Settings for saving of standart error.
  se <- ""
  if (SE == TRUE){
    se <- paste(gsub(pattern=".bam", replacement="", x=bamFile),
                "STDERR.MARKDUP.TXT",
                sep=".")
  }
  ### Running of command.
  system2(command=gatk_path,
          args=c("MarkDuplicates", input, output, metrics, remove.d),
          if (se == ""){
            stderr=se
          }else{
            stderr=paste(path, se, sep="/")
          })
}


#' Recalibrate a BAM file
#' @description This is a wrapper function that uses the low-level functions
#'     BaseRecalibrator() and ApplyBQSR() from package GATK to recalibrate
#'     sequencing quality scores of bases in BAM file of interest.
#' @param bamDir character string giving the name (or path to and name) of
#'     directory with BAM file(s). NULL by default, which means the current
#'     working directory will be used instead of the specified BAM directory.
#' @param bamFile character string giving the name of input BAM file.
#' @param vcfDir character string giving the name (or path to and name) of
#'     directory with VCF file(s). NULL by default, which means that the
#'     current working directory will be used instead of the specified VCF
#'     directory.
#' @param vcfFile character string giving the name of input VCF file.
#' @param fastaDir character string giving the name (or path to and name) of
#'     folder for storing of reference genome.
#' @param faFile charater string giving the name of file containing reference
#'     genome in format *.FASTA or *.FA.
#' @param createBQSR logical, TRUE by default. It allows to create a new BQSR
#'     table with provided BAM file, known variation sites and reference genome.
#' @param recalibrate logical, TRUE by default. It allows to run a recalibration
#'     process. For this purpose, the BQSR table provided trough nameBQSR or
#'     created by function itself if createBQSR=TRUE, will be used.
#' @param nameBQSR character string giving the name of pre-computed BQSR table.
#'     NULL by default.
#' @param postfix character string that will be added to the output file name
#'     as postfix. "Recalibrated" by default.
#' @param SE if TRUE (by default), standart error determined by function
#'     MarkDuplicates() will be writen in file STDERR.RECALIBRATION.TXT.
#' @param gatk_path character string giving the path to internal or external
#'     (system) GATK executive script. This path can be obtained by running
#'     Find.Gatk() or settled manually. Default value is "gatk", i.e. function
#'     will uses GATK systemly from $PATH.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @return if createBQSR=TRUE, a new BQSR table in working directory;
#'     if recalibrate=TRUE, a new BAM file with recalibrated sequencing quality
#'     scores of bases in working directory; if SE=TRUE, a file with standart
#'     error log in working directory.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @examples
#' ...
#' @export
recalibrateBamGATK <- function(bamDir=NULL,
                               bamFile,
                               vcfDir=NULL,
                               vcfFile,
                               fastaDir=NULL,
                               faFile,
                               createBQSR=TRUE,
                               recalibrate=TRUE,
                               nameBQSR,
                               postfix="Recalibrated",
                               SE=TRUE,
                               gatk_path="gatk",
                               workDir=NULL){
  ### Setting the working directory.
  if (is.null(x=workDir)){
    workDir <- getwd()
  }
  ### Full path to the BAM file.
  if (is.null(x=bamDir)){
    path <- workDir
  }else{
    path <- paste(workDir, bamDir, sep="/")
  }
  ### Full path to the VCF file.
  if (is.null(x=vcfDir)){
    path2 <- workDir
  }else{
    path2 <- paste(workDir, vcfDir, sep="/")
  }
  ### Full path to the FASTA/FA file.
  if (is.null(x=fastaDir)){
    path3 <- workDir
  }else{
    path3 <- paste(workDir, fastaDir, sep="/")
  }
  ### Generation of BQSR table.
  if (isTRUE(x=createBQSR)){
    bqsr_name <- paste(gsub(pattern=".bam", replacement="", x=bamFile),
                       "bqsr.table",
                       sep=".")
    system2(command=gatk_path,
            args=c("BaseRecalibrator",
                   paste("-I",
                         paste(path, bamFile, sep="/"),
                         sep = " "),
                   paste("-O",
                         paste(path, bqsr_name, sep="/"),
                         sep = " "),
                   paste("--known-sites",
                         paste(path2, vcfFile, sep="/"),
                         sep = " "),
                   paste("-R",
                         paste(path3, faFile, sep="/"),
                         sep = " ")))
  }
  ### Recalibration process.
  if (isTRUE(x=recalibrate)){
    if (isTRUE(x=createBQSR)){
      bqsr_name <- paste(gsub(pattern=".bam", replacement="", x=bamFile),
                         "bqsr.table",
                         sep=".")
    }else{
      bqsr_name <- nameBQSR
    }
    ### Settings for saving of standart error.
    se <- ""
    if (SE == TRUE){
      se <- paste(gsub(pattern=".bam", replacement="", x=bamFile),
                  "STDERR.RECALIBRATION.TXT",
                  sep=".")
    }
    ### Running of command.
    system2(command=gatk_path,
            args=c("ApplyBQSR",
                   paste("-bqsr", paste(path, bqsr_name, sep="/"), sep=" "),
                   paste("-I", paste(path, bamFile, sep="/"), sep=" "),
                   paste("-O", paste(path,
                                     paste(gsub(pattern=".bam",
                                                replacement="",
                                                x=bamFile),
                                           postfix,
                                           "bam",
                                           sep="."),
                                     sep="/"), sep=" ")))
    if (se == ""){
      stderr=se
    }else{
      stderr=paste(path, se, sep="/")
    }
  }
}


