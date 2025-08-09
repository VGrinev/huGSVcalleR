#' Call single nucleotide variations using different statistical tests
#' @description This function is a wrapper that allows to call and control
#'     various low-level functions for identification of single nucleotide
#'     variations from NGS data.
#' @param fr_data either a name of CSV file containing pileup statistics (as
#'     returned by collectBasesPileup() function) or an object of class
#'     data.frame with frequancy table.
#' @param criteria character vector specifying which test(-s) to be performed.
#'     The following values are possible: "binomial", "counting", "entropy"
#'     (by default), "fisher", "fisherSmyth", "gatk" and "poisson".
#' @param vcfDir character string giving the name (or path to and name) of
#'     directory for storing of generated VCF file(-s). NULL by default, which
#'     means the current working directory will be used instead of the
#'     specified BAM directory.
#' @param outputVCF character string giving the name of output VCF file.
#' @param bamDir character string giving the name (or path to and name) of
#'     directory with BAM file(-s). NULL by default, which means the current
#'     working directory will be used instead of the specified BAM directory.
#' @param bamFile character string giving the name of input BAM file.
#' @param fastaDir character string giving the name (or path to and name) of
#'     folder for storing of reference genome (with .fai indices if GATK).
#' @param faFile charater string giving the name of file containing reference
#'     genome in format *.FASTA or *.FA.
#' @param baseq integer value, minimal "base quality" for each nucleotide in an
#'     alignment. Default value is 20 that means read nucleotides with quality
#'     scores less than 20 will be excluded from analysis.
#' @param mindepth integer value giving the minimal coverage of location
#'     containing the variation. Default value is 1.
#' @param maxdepth integer value giving the maximal coverage of location
#'     containing the variation. Default value is 1e6. This option is useful
#'     for removing PCR artifacts.
#' @param qvalue numeric value giving the q-value cutoff for variations calling
#'     at sequencing depth 50X. Default value is 12. The q-value is calcuated
#'     as -log10(p), where p is a p-value yielded from the Fisher’s exact test.
#'     The function exactSNP() automatically adjusts the q-value cutoff for
#'     each chromosomal location according to its sequencing depth, based on
#'     this cutoff.
#' @param trim integer value giving the number of nucleotides trimmed off from
#'     each end of the read. Default value is 0.
#' @param threads numeric value giving the number of threads/CPUs used. Default
#'     value is 1.
#' @param postfix character string that will be added to the output VCF file
#'     name as postfix. NA by default at which the postfix "GATK based SNVs"
#'     will be used. Only with GATK.
#' @param intervals character string specifying the name of the tab-delimited
#'     TXT file containing coordinates of genomic location(-s) of interest. The
#'     default value is NULL. If so, this file must contains the following four
#'     fields:
#'     i) seqnames - the name of chromosome or scaffold with prefix "chr";
#'     ii) start   - start coordinate of genomic interval of interest;
#'     iii) end    - end coordinate of genomic interval of interest;
#'     iv) strand  - (optionally) strand information about genomic interval
#'                   of interest.
#'     TXT file with genomic location(-s) must be in working directory.
#' @param pcr_model GATK argument --pcr-indel-model. Possible values: NONE,
#'     AGRESSIVE, HOSTILE, CONSERVATIVE. See HaplotypeCaller
#'     at https://gatk.broadinstitute.org for further details.
#' @param soft_clip GATK argument --dont-use-soft-clipped. TRUE by default.
#' @param ERC GATK ERC mode. "NONE" by default. Possible values: NONE, GVCF or
#'     BP_RESOLUTION. See HaplotypeCaller at https://gatk.broadinstitute.org
#'     for further details.
#' @param cores number of CPU cores to use for Hidden Markov Models. If "MAX"
#'     (by default), all available cores will be detected and used.
#' @param SE if TRUE (by default), standart error determined by GATK will be
#'     writen in file STDERR.HC.TXT.
#' @param gatk_path character string giving the path to internal or external
#'     (system) GATK executive script. This path can be obtained by running
#'     Find.Gatk() or settled manually. Default value is "gatk", i.e. function
#'     will uses GATK systemly from $PATH.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @return VCF file(-s) with detected single nucleotide variations.
#' @author Dzianis D. Sarnatski, Mikalai M. Yatskou, Vasily V. Grinev
#' @examples
#' ...
#' @export
#' @importFrom data.table data.table
callSNVs <- function(fr_data,
                     criteria="entropy",
                     vcfDir=NULL,
                     outputVCF,
                     #common arguments to be passed to "fisherSmyth" and "gatk"
                     bamDir=NULL,
                     bamFile=NA,
                     fastaDir=NULL,
                     faFile=NA,
                     #specific arguments to be passed to "fisherSmyth"
                     baseq=20,
                     mindepth=1,
                     maxdepth=1e6,
                     qvalue=12,
                     trim=0,
                     threads=1,
                     #specific arguments to be passed to "gatk"
                     postfix=NA,
                     intervals=NULL,
                     pcr_model=NA,
                     soft_clip=TRUE,
                     ERC="NONE",
                     cores="MAX",
                     SE=TRUE,
                     gatk_path=NULL,
                     workDir=NULL){
  ### Validation of the input criteria.
  if (length(x=criteria) == 0){
    stop("No test(-s) will be performed. Please specify valid criteria")
  }
  invalid_criteria <- setdiff(x=criteria, y=SUPPORTED_CLASSICAL_TESTS)
  if (length(x=invalid_criteria) > 0){
    stop(paste("Invalid criteria specified:",
               paste(invalid_criteria, collapse=", ")))
  }
  ### Setting the working directory.
  if (is.null(x=workDir)){
    workDir <- getwd()
  }
  ### Classical tests.
  if (criteria %in% c("binomial", "counting", "entropy", "fisher", "poisson")){
    ##  Retrieving of the frequancy table.
    if (is.object(x=fr_data)){
      fr_table=fr_data
    }else{
      #   Full path to the CSV file.
      path <- paste(workDir, fr_data, sep="/")
      #   Retrieving and validation of the file extension.
      frt <- tools::file_ext(x=path)
      if (!frt %in% "csv"){
        stop("Invalid format of CSV file")
      }
      #   Does the file actually exist in the specified location?
      if (!file.exists(path)){
        stop(paste("File not found:", path))
      }
      #   Loading of the CSV file as an object of class data.frame.
      fr_table <- read.csv(file=path, header=TRUE)
    }
    ##  Initialization.
    results <- list()
    ##  Perform tests based on the specified criteria.
    tests <- list(counting=snp_counting_cpp_parallel,
                  binomial=snp_binomial_cpp_parallel,
                  entropy=snp_entropy_cpp_parallel,
                  fisher=snp_fisher_cpp_parallel,
                  poisson=snp_poisson_cpp_parallel)
    results <- lapply(X=criteria,
                      FUN=function(test){tests[[test]](fr_table)})
    names(x=results) <- criteria
    ##  Transformation and writing of the results to VCF file.
    input <- data.table(cbind(fr_table, results))
    generateVCF(input=input,
                version="4.0",
                vcfDir=vcfDir,
                outputVCF=outputVCF,
                workDir=workDir)
  }
  ### Gordon K. Smyth's based approach.
  if (criteria == "fisherSmyth"){
    fisherSmythSNVs <- fisherCalleR(bamDir=bamDir,
                                    bamFile=bamFile,
                                    fastaDir=fastaDir,
                                    faFile=faFile,
                                    baseq=baseq,
                                    mindepth=mindepth,
                                    maxdepth=maxdepth,
                                    qvalue=qvalue,
                                    trim=trim,
                                    threads=threads,
                                    workDir=workDir)
  }
  ### GATK approach.
  if (criteria == "gatk"){
    gatkSNVs <- gatkCalleR(bamDir=bamDir,
                           bamFile=bamFile,
                           postfix=postfix,
                           fastaDir=fastaDir,
                           faFile=faFile,
                           intervals=NULL,
                           pcr_model,
                           soft_clip=TRUE,
                           ERC="NONE",
                           cores="MAX",
                           SE=TRUE,
                           gatk_path=NULL,
                           workDir=workDir)
  }
}


#' Call of single nucleotide variations using Fisher’s exact test
#' @description This is a wrapper function that uses the low-level function
#'     exactSNP() from package Rsubread for accurate and efficient SNVs calling.
#' @param bamDir character string giving the name (or path to and name) of
#'     directory with BAM file(-s). NULL by default, which means the current
#'     working directory will be used instead of the specified BAM directory.
#' @param bamFile character string giving the name of input BAM file.
#' @param fastaDir character string giving the name (or path to and name) of
#'     folder for storing of reference genome.
#' @param faFile charater string giving the name of file containing reference
#'     genome in format *.FASTA or *.FA.
#' @param baseq integer value, minimal "base quality" for each nucleotide in an
#'     alignment. Default value is 20 that means read nucleotides with quality
#'     scores less than 20 will be excluded from analysis.
#' @param mindepth integer value giving the minimal coverage of location
#'     containing the variation. Default value is 1.
#' @param maxdepth integer value giving the maximal coverage of location
#'     containing the variation. Default value is 1e6. This option is useful
#'     for removing PCR artifacts.
#' @param qvalue numeric value giving the q-value cutoff for variations calling
#'     at sequencing depth 50X. Default value is 12. The q-value is calcuated
#'     as -log10(p), where p is a p-value yielded from the Fisher’s exact test.
#'     The function exactSNP() automatically adjusts the q-value cutoff for
#'     each chromosomal location according to its sequencing depth, based on
#'     this cutoff.
#' @param trim integer value giving the number of nucleotides trimmed off from
#'     each end of the read. Default value is 0.
#' @param threads numeric value giving the number of threads/CPUs used. Default
#'     value is 1.
#' @param workDir character string giving the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return file of VCF format is written to the current working directory.
#' @author Vasily V. Grinev.
#' @examples
#' res <- fisherCalleR(bamDir="Files_BAM",
#'                     bamFile="SRR21721515.bam",
#'                     fastaDir="Files_FASTA",
#'                     faFile="hg38.fa.gz",
#'                     baseq=20,
#'                     mindepth=1,
#'                     maxdepth=1e6,
#'                     qvalue=12,
#'                     trim=0,
#'                     threads=3,
#'                     workDir="D:/Vasily Grinev")
#' @export
#' @importFrom Rsubread exactSNP
fisherCalleR <- function(bamDir=NULL,
                         bamFile,
                         fastaDir=NULL,
                         faFile,
                         baseq=20,
                         mindepth=1, maxdepth=1e6,
                         qvalue=12,
                         trim=0,
                         threads=1,
                         workDir=NULL){
  ### Setting the working directory.
  if (is.null(x=workDir)){
    workDir <- getwd()
  }
  ### Full path to the BAM file.
  if (is.null(x=bamDir)){
    path1 <- workDir
  }else{
    path1 <- paste(workDir, bamDir, sep="/")
  }
  ### Full path to the reference genome.
  if (is.null(x=fastaDir)){
    path2 <- workDir
  }else{
    path2 <- paste(workDir, fastaDir, sep="/")
  }
  ### Setting the BAM file.
  bam <- paste(path1, bamFile, sep="/")
  ### Setting the FASTA file.
  fa <- paste(path2, faFile, sep="/")
  ### Setting the output VCF file.
  vcf <- gsub(pattern=".bam",
              replacement=", fisherCalleR, SNVs.vcf",
              x=bam)
  ### Detection of SNVs.
  SNVs <- exactSNP(readFile=bam,
                   isBAM=TRUE,
                   refGenomeFile=fa,
                   minBaseQuality=baseq,
                   minReads=mindepth,
                   maxReads=maxdepth,
                   qvalueCutoff=qvalue,
                   nTrimmedBases=trim,
                   nthreads=threads,
                   outputFile=vcf)
  ### Returning a final object.
  return(SNVs)
}


#' Call of single nucleotide variations using haplotype caller GATK
#' @description This is a wrapper function to run GATK within the R environment
#'     for accurate and efficient SNVs calling.
#' @param bamDir character string giving the name (or path to and name) of
#'     directory with BAM file(s). NULL by default, which means the current
#'     working directory will be used instead of the specified BAM directory.
#' @param bamFile character string giving the name of input BAM file.
#' @param postfix character string that will be added to the output VCF file
#'     name as postfix. NA by default at which the postfix "GATK based SNVs"
#'     will be used.
#' @param fastaDir character string giving the name (or path to and name) of
#'     folder for storing of fai-indexed reference genome.
#' @param faFile charater string giving the name of file containing reference
#'     genome in format *.FASTA or *.FA.
#' @param intervals character string specifying the name of the tab-delimited
#'     TXT file containing coordinates of genomic location(-s) of interest. The
#'     default value is NULL. If so, this file must contains the following four
#'     fields:
#'     i) seqnames - the name of chromosome or scaffold with prefix "chr";
#'     ii) start   - start coordinate of genomic interval of interest;
#'     iii) end    - end coordinate of genomic interval of interest;
#'     iv) strand  - (optionally) strand information about genomic interval
#'                   of interest.
#'     TXT file with genomic location(-s) must be in working directory.
#' @param pcr_model GATK argument --pcr-indel-model. Possible values: NONE,
#'     AGRESSIVE, HOSTILE, CONSERVATIVE. See HaplotypeCaller
#'     at https://gatk.broadinstitute.org for further details.
#' @param soft_clip GATK argument --dont-use-soft-clipped. TRUE by default.
#' @param ERC GATK ERC mode. "NONE" by default. Possible values: NONE, GVCF or
#'     BP_RESOLUTION. See HaplotypeCaller at https://gatk.broadinstitute.org
#'     for further details.
#' @param cores number of CPU cores to use for Hidden Markov Models. If "MAX"
#'     (by default), all available cores will be detected and used.
#' @param SE if TRUE (by default), standart error determined by GATK will be
#'     writen in file STDERR.HC.TXT.
#' @param gatk_path character string giving the path to internal or external
#'     (system) GATK executive script. This path can be obtained by running
#'     Find.Gatk() or settled manually. Default value is "gatk", i.e. function
#'     will uses GATK systemly from $PATH.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @details GATK performs SNVs and short indels calling by local de novo
#'     assembly of mapped reads
#' @return a new VCF file with detected SNVs.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @examples
#' ...
#' @export
#' @importFrom parallel detectCores
gatkCalleR <- function(bamDir=NULL,
                       bamFile,
                       postfix=NA,
                       fastaDir=NULL,
                       faFile,
                       intervals=NULL,
                       pcr_model,
                       soft_clip=TRUE,
                       ERC="NONE",
                       cores="MAX",
                       SE=TRUE,
                       gatk_path=NULL,
                       workDir=NULL){
  ### Validation of GATK path.
  if(is.null(x=gatk_path)){
    return("Cannot run GATK without path. Please provide full path to GATK
                folder and execution script! Format of
                path: /home/Bioinformatics/Soft/gatk-4.2.6.1/gatk. Or use a
                just *gatk*, if it can be called directly from bash console")
  }
  ### Setting the working directory.
  if (is.null(x=workDir)){
    workDir <- getwd()
  }
  ### Full path to the BAM file.
  if (is.null(x=bamDir)){
    path1 <- workDir
  }else{
    path1 <- paste(workDir, bamDir, sep="/")
  }
  ### Full path to the reference genome.
  if (is.null(x=fastaDir)){
    path2 <- workDir
  }else{
    path2 <- paste(workDir, fastaDir, sep="/")
  }
  ### Settings of run.
  #   Settings for input BAM file.
  input <- paste("-I",
                 paste(path1, bamFile, sep="/"),
                 sep=" ")
  #   Settings for output VCF file.
  if (!is.na(x=postfix)){
    output <- paste("-O",
                    paste(path1,
                          paste(gsub(pattern=".bam",
                                     replacement="",
                                     x=bamFile),
                                postfix,
                                "vcf",
                                sep="."),
                          sep="/"),
                    sep=" ")
  }else{
    output <- paste("-O",
                    paste(path1,
                          paste(gsub(pattern=".bam",
                                     replacement="",
                                     x=bamFile),
                                "GATK based SNVs",
                                "vcf",
                                sep="."),
                          sep="/"),
                    sep=" ")
  }
  #   Settings for fai-indexed reference genome.
  fa <- paste("-R", paste(path2, faFile, sep="/"), sep=" ")
  #   Settings for intervals.
  if (!is.null(x=intervals)){
    GRs <- read.table(file=paste(workDir, intervals, sep="/"),
                      sep="\t",
                      header=TRUE,
                      quote="\"",
                      as.is=TRUE)
    GRs <- paste(paste(GRs[, 1], GRs[, 2], sep=":"), GRs[, 3], sep="-")
    GRs <- paste("--intervals", GRs)
    fa <- paste(fa, GRs)
  }
  #   Settings for PCR indel model.
  pcr <- paste("--pcr-indel-model", pcr_model)
  #   Settings for using for discovery the soft clipped regions
  softclip <- paste("--dont-use-soft-clipped-bases", soft_clip)
  # Settings for ERC output.
  erc <- paste("-ERC", ERC)
  #   Setting the multithreads process.
  if (cores == "MAX"){
    cores <- detectCores()
  }
  cores <- paste("--native-pair-hmm-threads", cores)
  ### Settings for saving of standart error.
  se <- ""
  if (SE == TRUE){
    se <- paste(gsub(pattern=".bam", replacement="", x=bamFile),
                "STDERR.HC.TXT",
                sep=".")
  }
  ### Running of command.
  system2(command=gatk_path,
          args=c("HaplotypeCaller", input, output, fa, pcr, softclip, erc, cores),
          if (se == ""){
            stderr=se
          }else{
            stderr=paste(path1, se, sep="/")
          })
}


#' Generate VCF file based on call table
#' @description This function converts call table to standard VCF file.
#' @param input character string giving the name of tab-delimited TXT file or
#'     R object of class data.frame containing call data. File or object must
#'     contains the following fields:
#'     i) position   - location of single nucleotide variant in format
#'                     chrZ:coordinate;
#'     ii) reference - reference variant of nucleotide;
#'     iii) A        - number of redas supporting A in a given position;
#'     iv) C         - number of redas supporting C in a given position;
#'     v) G          - number of redas supporting G in a given position;
#'     vi) T         - number of redas supporting T in a given position;
#'     vii-xi) ...   - one to five columns with statistical test results.
#'     If file, it must be located in working directory.
#' @param version character string specifying the version of VCF format.
#'     Default value is "4.0". It is planned that other versions of the format
#'     will be added in the future.
#' @param vcfDir character string giving the name (or path to and name) of
#'     directory with VCF file(-s). NULL by default, which means the current
#'     working directory will be used instead of the specified VCF directory.
#' @param outputVCF character string giving the name of output VCF file.
#' @param workDir character string giving the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return file(-s) of VCF format is/are written to the current working
#'     directory or specified VCF directory.
#' @author Vasily V. Grinev.
#' @examples
#' generateVCF(input="example_seq.call_table.txt",
#'             version="4.0",
#'             vcfDir="Files_VCF",
#'             outputVCF="example_seq.call_table.identified_SNVs",
#'             workDir="D:/Vasily Grinev")
#' @export
#' @importFrom data.table fread
generateVCF <- function(input,
                        version="4.0",
                        vcfDir=NULL,
                        outputVCF,
                        workDir=NULL){
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
  ### Retrieving of the input data.
  if (is.object(x=input)){
    vcf_pre=input
  }else{
    vcf_pre <- fread(file=paste(workDir, input, sep="/"))
  }
  tests <- colnames(x=vcf_pre)[-1:-6]
  ### Transformation of data into a frame of the appropriate structure.
  vcf_pre$DP <- rowSums(x=vcf_pre[, 3:6])
  A <- vcf_pre[vcf_pre$A > 0 & vcf_pre$reference != "A", c(-4:-6)]
  A$ALT <- "A"
  A$MM <- A$A
  A <- A[, -3]
  C <- vcf_pre[vcf_pre$C > 0 & vcf_pre$reference != "C", c(-3, -5:-6)]
  C$ALT <- "C"
  C$MM <- C$C
  C <- C[, -3]
  G <- vcf_pre[vcf_pre$G > 0 & vcf_pre$reference != "G", c(-3:-4, -6)]
  G$ALT <- "G"
  G$MM <- G$G
  G <- G[, -3]
  T <- vcf_pre[vcf_pre$T > 0 & vcf_pre$reference != "T", c(-3:-5)]
  T$ALT <- "T"
  T$MM <- T$T
  T <- T[, -3]
  vcf_pre <- rbind(A, C, G, T)
  if (version == "4.0"){
    VCF <- unique(vcf_pre[, 1:(ncol(x=vcf_pre) - 2)])
    VCF <- VCF[order(x=VCF$position), ]
    ALT <- vcf_pre[, paste(ALT, collapse=","), by=position]
    ALT <- ALT[order(x=ALT$position), ]
    VCF$ALT <- ALT$V1
    MMsum <- vcf_pre[, sum(MM), by=position]
    MMsum <- MMsum[order(x=MMsum$position), ]
    VCF$MMsum <- MMsum$V1
    MM <- vcf_pre[, paste(MM, collapse=","), by=position]
    MM <- MM[order(x=MM$position), ]
    VCF$MM <- MM$V1
    VCF$CHROM <- gsub(pattern=":.+", replacement="", x=VCF$position)
    VCF$POS <- as.numeric(x=gsub(pattern=".+:",
                                 replacement="",
                                 x=VCF$position))
    VCF$ID <- paste(paste(VCF$position,
                          VCF$reference,
                          sep="_"),
                    VCF$ALT,
                    sep="/")
    VCF$QUAL <- "."
    VCF$FILTER <- "."
    VCF$INFO <- paste(paste0("DP=", VCF$DP),
                      paste0("MMsum=", VCF$MMsum),
                      paste0("MM=", VCF$MM),
                      sep=";")
    VCF$REF <- VCF$reference
    VCF <- data.frame(VCF)
    VCF <- VCF[, c("CHROM", "POS", "ID",
                   "REF", "ALT",
                   "QUAL",
                   "FILTER",
                   "INFO",
                   tests)]
    ### Writing of VCF file(-s).
    meta <- c("##fileformat=VCFv4.0",
              "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">",
              "##INFO=<ID=MMsum,Number=1,Type=Integer,Description=\"Number of supporting reads for alternative allele(-s)\">",
              "##INFO=<ID=MM,Number=1,Type=String,Description=\"Number of supporting reads for each alternative allele\">",
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    for (i in 1:length(x=tests)){
      VCF$QUAL <- round(x=VCF[, tests[i]], digits=4)
      vcf_name <- paste(paste(vcfDir,
                              paste(outputVCF,
                                    tests[i],
                                    sep="."),
                              sep="/"),
                        "vcf",
                        sep=".")
      empty_file <- file(description=vcf_name, open="wt")
      writeLines(text=meta, con=empty_file)
      write.table(x=rbind(VCF[1:6, 1:8], VCF[, 1:8]),
                  file=vcf_name,
                  sep="\t",
                  quote=FALSE,
                  col.names=FALSE,
                  row.names=FALSE)
      close(con=empty_file)
      cleaning <- readLines(con=vcf_name)
      writeLines(text=cleaning[-6], con=vcf_name)
    }
  }
}


