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
#  Last updated: July 28, 2025.

generateVCF <- function(input,
                        version="4.0",
                        vcfDir=NULL,
                        outputVCF,
                        workDir=NULL){
    ### Loading the required package.
    #   This code was successfully tested with the package data.table v.1.15.4.
    suppressMessages(expr=library(package=data.table))
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
