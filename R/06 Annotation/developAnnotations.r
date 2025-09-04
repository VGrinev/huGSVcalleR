#' Develop Ensembl-like annotations of experimental transcriptome
#' @description This function converts a standard GTF/GFF file with experimental
#'     transcriptome to Ensembl-like annotations.
#' @param gtfDir character string giving the name (or path to and name) of
#'     directory with GTF file(-s). NULL by default, which means the current
#'     working directory will be used instead of the specified GTF directory.
#' @param gtfFile character string giving the name of input GTF/GFF file with
#'     experimental transcriptome. It is typically StringTie output. Valid
#'     format is "gtf" or "gff".
#' @param orf character string giving the name of tab-delimited TXT file with
#'     coordinates of open reading frames in experimental transcripts. It is
#'     typically ORFhunteR output. This file should include three fields:
#'     i) transcript_id - transcript ID;
#'     ii) orf_start    - start coordinate of the open reading frame in a
#'                        transcript;
#'     iii) orf_end     - end coordinate of the open reading frame in a
#'                        transcript.
#' @param gtf character string giving the name of output GTF/GFF file with
#'     Ensembl-like annotated experimental transcriptome.
#' @param src character string describing the origin of the GTF/GFF file.
#' @param organism the proper scientific name of the organism for which
#'     transcriptomic data were generated.
#' @param sqlite character string giving the name of output SQLite database
#'     with Ensembl-like annotated experimental transcriptome.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return GTF/GFF file and local SQLite database with Ensembl-like annotations
#'     of the experimental transcriptome.
#' @author Vasily V. Grinev.
#' @examples
#' anno <- developAnnotations(gtfDir="Files_GTF",
#'                            gtfFile="test.data_6.gtf",
#'                            orf="test.data_6.txt",
#'                            gtf="test.data_6.anno.gtf",
#'                            src="StringTie",
#'                            organism="Homo sapiens",
#'                            sqlite="test.data_6.anno.sqlite",
#'                            workDir="D:/Vasily Grinev")
#' @export
#  Last updated: September 4, 2025.

developAnnotations <- function(gtfDir=NULL,
                               gtfFile,
                               orf,
                               gtf,
                               src,
                               organism,
                               sqlite,
                               workDir=NULL){
    ### Loading of required packages.
    #   This code was successfully tested with packages GenomicFeatures
    #   v.1.40.1 and rtracklayer v.1.48.0.
    suppressMessages(expr=library(package=GenomicFeatures))
    suppressMessages(expr=library(package=rtracklayer))
    ### Setting the working directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Setting the GTF directory.
    if (is.null(x=gtfDir)){
        gtfDir <- workDir
    }else{
        gtfDir <- paste(workDir, gtfDir, sep="/")
    }
    ### Loading of the experimental annotations as an object of class GRanges.
    GTFs <- import(con=paste(gtfDir, gtfFile, sep="/"))
    ### Loading of the ORF coordinates as an object of class data frame.
    ORFs <- read.table(file=paste(workDir, orf, sep="/"),
                       sep="\t",
                       header=TRUE,
                       quote="\"",
                       as.is=TRUE)
    ### Development of the Ensembl-like annotations.
    GTFs$gene_biotype <- "noncoding"
    GTFs$transcript_biotype <- "noncoding"
    GTFs[GTFs$transcript_id %in% ORFs$transcript_id, ]$transcript_biotype <-
                                                               "protein_coding"
    GTFs[GTFs$gene_id %in%
         GTFs[GTFs$transcript_biotype == "protein_coding", ]$gene_id,
                                             ]$gene_biotype <- "protein_coding"
    genes <- split(x=GTFs, f=GTFs$gene_id)
    genes <- lapply(X=genes,
                    FUN=function(y){start(x=y) <- min(x=start(x=y))
                                    end(x=y) <- max(x=end(x=y))
                                    y <- y[1]
                                    y$type <- "gene"
                                    mcols(x=y) <-
                                           mcols(x=y)[!names(x=mcols(x=y)) %in%
                                    c("transcript_id", "exon_number",
                                      "transcript_biotype")]
                                    return(y)})
    coding <- GTFs[GTFs$transcript_biotype == "protein_coding", ]
    coding <- split(x=coding, f=coding$transcript_id)
    for (i in 1:length(x=coding)){
        tr <- coding[[i]]
        orf_start <- ORFs[ORFs$transcript_id == tr$transcript_id[1], ]$orf_start
        orf_end <- ORFs[ORFs$transcript_id == tr$transcript_id[1], ]$orf_end
        if (as.character(x=strand(x=tr))[1] == "+"){
            if (orf_start > width(x=tr)[1]){
                m1 <- min(x=which(x=cumsum(x=width(x=tr)) - orf_start > 0))
                m2 <- max(x=cumsum(x=width(x=tr)[1:(m1 - 1)]))
                start_codon <- start(x=tr)[m1] + orf_start - m2 - 1
                start_codon <- c(start_codon, start_codon + 2)
                m1 <- min(x=which(x=cumsum(x=width(x=tr)) - orf_end > 0))
                m2 <- max(x=cumsum(x=width(x=tr)[1:(m1 - 1)]))
                stop_codon <- start(x=tr)[m1] + orf_end - m2 - 3
                stop_codon <- c(stop_codon, stop_codon + 2)
            }else{
                start_codon <- start(x=tr)[1] + orf_start - 1
                start_codon <- c(start_codon, start_codon + 2)
                if (length(x=tr) == 1){
                    stop_codon <- start(x=tr)[1] + orf_end - 1
                    stop_codon <- c(stop_codon - 2, stop_codon)
                }else{
                    m1 <- min(x=which(x=cumsum(x=width(x=tr)) - orf_end > 0))
                    if (m1 == 1){
                        stop_codon <- start(x=tr)[1] + orf_end - 1
                        stop_codon <- c(stop_codon - 2, stop_codon)
                    }else{
                        m2 <- max(x=cumsum(x=width(x=tr)[1:(m1 - 1)]))
                        stop_codon <- start(x=tr)[m1] + orf_end - m2 - 3
                        stop_codon <- c(stop_codon, stop_codon + 2)
                    }
                }
            }
            tr <- sort(x=tr)
            CDS <- c(min(x=start_codon), min(x=stop_codon) - 1)
            CDS <- data.frame(seqnames=as.character(x=seqnames(x=tr))[1],
                              start=CDS[1],
                              end=CDS[2],
                              strand=as.character(x=strand(x=tr))[1])
            CDS <- makeGRangesFromDataFrame(df=CDS)
            hits <- findOverlaps(query=CDS, subject=tr, type="any")
            CDS_exons <- tr[subjectHits(x=hits)]
            CDS_exons$type <- "CDS"
            start(x=CDS_exons)[1] <- start(x=CDS)
            end(x=CDS_exons)[length(x=CDS_exons)] <- end(x=CDS)
            start <- data.frame(seqnames=as.character(x=seqnames(x=tr))[1],
                                start=start_codon[1],
                                end=start_codon[2],
                                strand=as.character(x=strand(x=tr))[1])
            start <- makeGRangesFromDataFrame(df=start)
            hits <- findOverlaps(query=start, subject=tr, type="any")
            start_codon <- tr[subjectHits(x=hits)]
            start_codon$type <- "start_codon"
            start(x=start_codon)[1] <- start(x=start)
            end(x=start_codon)[1] <- end(x=start)
            stop <- data.frame(seqnames=as.character(x=seqnames(x=tr))[1],
                               start=stop_codon[1],
                               end=stop_codon[2],
                               strand=as.character(x=strand(x=tr))[1])
            stop <- makeGRangesFromDataFrame(df=stop)
            hits <- findOverlaps(query=stop, subject=tr, type="any")
            stop_codon <- tr[subjectHits(x=hits)]
            stop_codon$type <- "stop_codon"
            start(x=stop_codon)[1] <- start(x=stop)
            end(x=stop_codon)[1] <- end(x=stop)
            UTRs <- disjoin(c(tr, CDS))
            hits <- findOverlaps(query=UTRs, subject=CDS, type="any")
            UTRs <- UTRs[-queryHits(x=hits)]
            five_prime_utr <- UTRs[!UTRs > CDS]
            if (length(x=five_prime_utr) > 0){
                hits <- findOverlaps(query=five_prime_utr,
                                     subject=tr,
                                     type="any")
                mcols(x=five_prime_utr) <- mcols(x=tr[subjectHits(x=hits)])
                five_prime_utr$type <- "five_prime_utr"
            }
            three_prime_utr <- UTRs[UTRs > CDS]
            if (length(x=three_prime_utr) > 0){
                hits <- findOverlaps(query=three_prime_utr,
                                     subject=tr,
                                     type="any")
                mcols(x=three_prime_utr) <- mcols(x=tr[subjectHits(x=hits)])
                three_prime_utr$type <- "three_prime_utr"
                start(x=three_prime_utr) <- start(x=three_prime_utr) + 3
            }
            transcript <- tr[1]
            mcols(x=transcript) <-
             mcols(x=transcript)[names(x=mcols(x=transcript)) != "exon_number"]
            transcript$type <- "transcript"
            start(x=transcript) <- min(x=start(x=tr))
            end(x=transcript) <- max(x=end(x=tr))
            gene <- genes[names(x=genes) == unique(x=tr$gene_id)][[1]]
            tr <- c(gene,
                    transcript,
                    sort(c(tr, five_prime_utr, start_codon, CDS_exons,
                           stop_codon, three_prime_utr)))
        }
        if (as.character(x=strand(x=tr))[1] == "-"){
            tr <- sort(x=tr, decreasing=TRUE)
            if (orf_start > width(x=tr)[1]){
                m1 <- min(x=which(x=cumsum(x=width(x=tr)) - orf_start > 0))
                m2 <- max(x=cumsum(x=width(x=tr)[1:(m1 - 1)]))
                start_codon <- end(x=tr)[m1] - (orf_start - m2) + 1
                start_codon <- c(start_codon - 2, start_codon)
                m1 <- min(x=which(x=cumsum(x=width(x=tr)) - orf_end > 0))
                m2 <- max(x=cumsum(x=width(x=tr)[1:(m1 - 1)]))
                stop_codon <- end(x=tr)[m1] - (orf_end - m2) + 1
                stop_codon <- c(stop_codon, stop_codon + 2)
            }else{
                start_codon <- end(x=tr)[1] - orf_start + 1
                start_codon <- c(start_codon - 2, start_codon)
                if (length(x=tr) == 1){
                    stop_codon <- end(x=tr)[1] - orf_end + 1
                    stop_codon <- c(stop_codon, stop_codon + 2)
                }else{
                    if (tail(x=cumsum(x=width(x=tr)), n=1) == orf_end){
                        stop_codon <- start(x=tail(x=tr, n=1))
                        stop_codon <- c(stop_codon, stop_codon + 2)
                    }else{
                        m1 <- min(x=which(x=cumsum(x=width(x=tr)) - orf_end > 0))
                        if (m1 == 1){
                            stop_codon <- end(x=tr)[1] - orf_end + 1
                            stop_codon <- c(stop_codon, stop_codon + 2)
                        }else{
                            m2 <- max(x=cumsum(x=width(x=tr)[1:(m1 - 1)]))
                            stop_codon <- end(x=tr)[m1] - (orf_end - m2) + 1
                            stop_codon <- c(stop_codon, stop_codon + 2)
                        }
                    }
                }
            }
            tr <- sort(x=tr)
            tr$exon_number <- rev(x=tr$exon_number)
            CDS <- c(max(x=stop_codon) + 1, max(x=start_codon))
            CDS <- data.frame(seqnames=as.character(x=seqnames(x=tr))[1],
                              start=CDS[1],
                              end=CDS[2],
                              strand=as.character(x=strand(x=tr))[1])
            CDS <- makeGRangesFromDataFrame(df=CDS)
            hits <- findOverlaps(query=CDS, subject=tr, type="any")
            CDS_exons <- tr[subjectHits(x=hits)]
            CDS_exons$type <- "CDS"
            start(x=CDS_exons)[1] <- start(x=CDS)
            end(x=CDS_exons)[length(x=CDS_exons)] <- end(x=CDS)
            start <- data.frame(seqnames=as.character(x=seqnames(x=tr))[1],
                                start=start_codon[1],
                                end=start_codon[2],
                                strand=as.character(x=strand(x=tr))[1])
            start <- makeGRangesFromDataFrame(df=start)
            hits <- findOverlaps(query=start, subject=tr, type="any")
            start_codon <- tr[subjectHits(x=hits)]
            start_codon$type <- "start_codon"
            start(x=start_codon)[1] <- start(x=start)
            end(x=start_codon)[1] <- end(x=start)
            stop <- data.frame(seqnames=as.character(x=seqnames(x=tr))[1],
                               start=stop_codon[1],
                               end=stop_codon[2],
                               strand=as.character(x=strand(x=tr))[1])
            stop <- makeGRangesFromDataFrame(df=stop)
            hits <- findOverlaps(query=stop, subject=tr, type="any")
            stop_codon <- tr[subjectHits(x=hits)]
            stop_codon$type <- "stop_codon"
            start(x=stop_codon)[1] <- start(x=stop)
            end(x=stop_codon)[1] <- end(x=stop)
            UTRs <- disjoin(c(tr, CDS))
            hits <- findOverlaps(query=UTRs, subject=CDS, type="any")
            UTRs <- UTRs[-queryHits(x=hits)]
            five_prime_utr <- UTRs[UTRs > CDS]
            if (length(x=five_prime_utr) > 0){
                hits <- findOverlaps(query=five_prime_utr,
                                     subject=tr,
                                     type="any")
                mcols(x=five_prime_utr) <- mcols(x=tr[subjectHits(x=hits)])
                five_prime_utr$type <- "five_prime_utr"
            }
            three_prime_utr <- UTRs[!UTRs > CDS]
            if (length(x=three_prime_utr) > 0){
                hits <- findOverlaps(query=three_prime_utr,
                                     subject=tr,
                                     type="any")
                mcols(x=three_prime_utr) <- mcols(x=tr[subjectHits(x=hits)])
                three_prime_utr$type <- "three_prime_utr"
                end(x=three_prime_utr) <- end(x=three_prime_utr) - 3
            }
            transcript <- tr[1]
            mcols(x=transcript) <-
             mcols(x=transcript)[names(x=mcols(x=transcript)) != "exon_number"]
            transcript$type <- "transcript"
            start(x=transcript) <- min(x=start(x=tr))
            end(x=transcript) <- max(x=end(x=tr))
            gene <- genes[names(x=genes) == unique(x=tr$gene_id)][[1]]
            tr <- c(gene,
                    transcript,
                    sort(c(tr, five_prime_utr, start_codon, CDS_exons,
                           stop_codon, three_prime_utr), decreasing=TRUE))
        }
        coding[[i]] <- tr
    }
    coding <- unlist(x=coding)
    noncoding <- GTFs[GTFs$transcript_biotype == "noncoding", ]
    noncoding <- split(x=noncoding, f=noncoding$transcript_id)
    for (i in 1:length(x=noncoding)){
        tr <- noncoding[[i]]
        if (as.character(x=strand(x=tr))[1] == "+"){
            tr <- sort(x=tr)
            transcript <- tr[1]
            mcols(x=transcript) <-
             mcols(x=transcript)[names(x=mcols(x=transcript)) != "exon_number"]
            transcript$type <- "transcript"
            start(x=transcript) <- min(x=start(x=tr))
            end(x=transcript) <- max(x=end(x=tr))
            gene <- genes[names(x=genes) == unique(x=tr$gene_id)][[1]]
            tr <- c(gene, transcript, tr)
        }
        if (as.character(x=strand(x=tr))[1] == "-"){
            tr <- sort(x=tr)
            tr$exon_number <- rev(x=tr$exon_number)
            transcript <- tr[1]
            mcols(x=transcript) <-
             mcols(x=transcript)[names(x=mcols(x=transcript)) != "exon_number"]
            transcript$type <- "transcript"
            start(x=transcript) <- min(x=start(x=tr))
            end(x=transcript) <- max(x=end(x=tr))
            gene <- genes[names(x=genes) == unique(x=tr$gene_id)][[1]]
            tr <- c(gene, transcript, tr)
        }
        noncoding[[i]] <- tr
    }
    noncoding <- unlist(x=noncoding)
    GTFs <- c(coding, noncoding)
    GTFs <- split(x=GTFs, f=GTFs$gene_id)
    for (i in 1:length(x=GTFs)){
        GTFs[[i]] <- c(GTFs[[i]][1],
                       GTFs[[i]][GTFs[[i]]$type != "gene", ])
        names(x=GTFs[[i]]) <- NULL
    }
    GTFs <- unlist(x=GTFs)
    ### Export of the developed annotations as file in format GTF/GFF.
    export(object=GTFs, con=paste(gtfDir, gtf, sep="/"))
    ### Export of the developed annotations as SQLite database.
    TxDb <- makeTxDbFromGFF(file=paste(gtfDir, gtf, sep="/"),
                            format="auto",
                            dataSource=src,
                            organism=organism)
    saveDb(x=TxDb, file=paste(gtfDir, sqlite, sep="/"))
}
