#' Annotate of single nucleotide variations
#' @description This function performs a multivariate annotation of the single
#'     nucleotide variations given in the standard VCF format.
#' @param vcfDir character string giving the name (or path to and name) of
#'     directory with VCF file(-s). NULL by default, which means the current
#'     working directory will be used instead of the specified VCF directory.
#' @param vcfFile character string giving the name of input VCF file. File must
#'     not contain NON_REF fields (for HaplotypeCaller it means run with -ERC NONE).
#' @param ref_genes character string giving the name of TXT file in
#'     tab-delimited format with basic annotations of reference genes. This
#'     file must contains the following fields:
#'     i) gene_id        - (optionally) gene ID (for example, Ensembl based);
#'     ii) gene_name     - name of gene (for example, HUGO based);
#'     iii) gene_biotype - biotype of gene;
#'     iv) seqnames      - name of chromosome or scaffold with prefix "chr";
#'     v) start          - start coordinate of gene;
#'     vi) end           - end coordinate of gene;
#'     vii) strand       - (optionally) strand information about gene.
#' @param output character string giving the name (without extension) of TXT
#'     file in tab-delimited format for storing annotated single nucleotide
#'     variations. The postfix "basic" will be added to this name for a file
#'     with a basic annotation and "extended" for a file with annotation
#'     according to the dbSNP. NULL by default, which means that the name
#'     will be generated automatically based on the name of the input
#'     VCF file.
#' @param annotateByDB logical. If TRUE, extended annotation by dbSNP
#'     will be performed and saved in second slot of output object. Also,
#'     tab-delimited txt will be created. FALSE by default.
#'     WARNING: The current version of this feature requires a significant
#'     amount of RAM to load full dbSNP VCF. It is not recommended to run it
#'     on a computer with less than 20 GB RAM.
#' @param dbDir character string giving the name (or path to and name) of
#'     directory with VCF file of dbSNP release. NULL by default, which means the current
#'     working directory will be used instead of the specified directory. Ignored when
#'     annotateByDB is FALSE.
#' @param dbFile character string giving the name of dbSNP VCF file for
#'     extended annotation. Ignored when annotateByDB is FALSE.
#' @param workDir character string giving the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return an object of class list with two slots and two new TXT file in
#' tab-delimited format with basic and extended annotation of provided single nucleotide variations.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @examples
#' genes <- "Ensembl release 114, GRCh38.p14, annotations of genes.txt"
#' annoSNVs <- annotateVariants(vcfDir="Files_VCF",
#'                              vcfFile="example_seq.vcf",
#'                              ref_genes,
#'                              output="example_seq, annotated SNVs",
#'                              workDir="D:/Vasily Grinev")
#' @export
#' @importFrom VariantAnnotation readVcf info
#' @importFrom SparseArray rowRanges
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom methods setClassUnion
annotateVariants <- function(vcfDir=NULL,
                             vcfFile,
                             ref_genes,
                             output=NULL,
                             annotateByDB = FALSE,
                             dbDir = NULL,
                             dbFile = NULL,
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
  ### Loading of the single nucleotide variations as an object of class
  #   CollapsedVCF.
  setClassUnion("ExpData", c("matrix", "SummarizedExperiment"))
  vcf <- readVcf(file=paste(vcfDir, vcfFile, sep="/"))
  ### Loading of the reference genes as an object of class GRanges.
  genes <- read.table(file=paste(workDir, ref_genes, sep="/"),
                      sep="\t",
                      header=TRUE,
                      quote="\"",
                      as.is=TRUE)
  genes <- makeGRangesFromDataFrame(df=genes, keep.extra.columns=TRUE)
  ### Annotation of the single nucleotide variations with gene names.
  mcols(x=rowRanges(x=vcf))$gene_name <- ""
  hits <- findOverlaps(query=rowRanges(x=vcf),
                       subject=genes[!is.na(x=genes$gene_name)],
                       type="within",
                       ignore.strand=TRUE)
  hits <- data.table(cbind(hits@from,
                           genes[!is.na(x=genes$gene_name)][hits@to, ]$gene_name))
  hits <- hits[, paste(V2, collapse=","), by=V1]
  colnames(x=hits) <- c("V1", "V2")
  hits$V1 <- as.numeric(x=hits$V1)
  hits$V2 <- unlist(x=lapply(X=hits$V2,
                             FUN=function(y){paste0(sort(x=unique(x=strsplit(x=y,
                                                                             split=",")[[1]])),
                                                    collapse=",")}))
  rowRanges(x=vcf)[hits$V1, ]$gene_name <- as.vector(x=hits$V2)
  ### Building an object with basic annotation
  SNVs <- cbind(data.frame(rowRanges(x=vcf)), info(x=vcf))
  if (annotateByDB == TRUE){
    alt.exp <- SNVs$ALT
  }
  SNVs$ALT <- unlist(x=lapply(X=SNVs$ALT, paste0, collapse=","))
  SNVs$snv_id <- paste(paste(paste(SNVs$seqnames, SNVs$start, sep=":"),
                             SNVs$REF, sep="_"),
                       SNVs$ALT, sep="/")
  SNVs <- SNVs[, c("snv_id", "seqnames", "start", "end", "strand", "REF",
                   "ALT",
                   "DP", "MM",
                   "QUAL", "gene_name")]
  colnames(x=SNVs) <- c(colnames(x=SNVs)[1:5],
                        "ref", "alt",
                        "pos_depth", "alt_depth",
                        "score", "gene_name")
  SNVs$seqnames <- as.character(x=SNVs$seqnames)
  SNVs$strand <- as.character(x=SNVs$strand)
  SNVs$alt_depth <- as.numeric(x=SNVs$alt_depth)
  out.anno <- list(no_matches = NULL, matches = NULL)
  out.anno[[1]] <- SNVs
  names(out.anno) <- c("No dbSNP matches", "dbSNP matches")
  ### Extended annotation by dbSNP.
  ### Cleaning RAM from previous step VCF
  if (annotateByDB == TRUE){
    rm(vcf)
    gc()
  }
  ### Loading of dbSNP annotation and finding overlaps with SNVs
  ### from previous steps by coordinats
  ref.snp <- readVcf(file = paste(vcfDir, dbFile, sep = "/"),
                     param = ScanVcfParam(info = c("RS", "CAF", "dbSNPBuildID", "SAO",
                                                   "SLO", "PM", "VP", "OTH")))
  hits.ref <- findOverlaps(subject = GRanges(seqnames = SNVs$seqnames,
                                             ranges = IRanges(start = SNVs$start,
                                                              end = SNVs$end)),
                           query = rowRanges(ref.snp), type = "equal")
  hits.ref <- cbind(hits.ref@from, as.character(rowRanges(ref.snp)@seqnames)[hits.ref@to],
                    start(rowRanges(ref.snp)@ranges[hits.ref@to]),
                    end(rowRanges(ref.snp)@ranges[hits.ref@to]),
                    fixed(ref.snp)[hits.ref@to, c("REF", "ALT")],
                    info(ref.snp)[hits.ref@to, ])
  colnames(hits.ref)[1:4] <- c("Overlap.SNVs", "seqnames", "start", "end")
  ### Comparison ALT field with reference SNVs
  alt.hits <- mapply(function(exp, ref) {
    any(exp %in% ref)
  }, alt.exp[1:10000],
  hits.ref$ALT[1:10000])
  hits.ref <- hits.ref[alt.hits,]
  ### Building frame with extended annotation
  hits.ref$ALT <- unlist(x=lapply(X=hits.ref$ALT, paste0, collapse=","))
  hits.ref$CAF <- unlist(x=lapply(X=hits.ref$CAF, paste0, collapse=","))
  colnames(hits.ref)[6] <- "alt_dbSNP"
  annot.SNVs <- cbind(SNVs[hits.ref$Overlap.SNVs,],
                      hits.ref[,6:ncol(hits.ref)])
  ### Removing rows with dbSNP matches from first table
  out.anno[[1]] <- out.anno[[1]][!out.anno[[1]]$snv_id %in% annot.SNVs$snv_id,]
  out.anno[[2]] <- annot.SNVs
  ###Writing output files
  if (is.null(x=output)){
    output_file <- paste(workDir,
                         sub(pattern=".vcf",
                             replacement=", based annotation.txt",
                             x=vcfFile),
                         sep="/")
  }else{
    output_file <- paste(workDir, paste(output, "based annotation", "txt", sep="."), sep="/")
  }
  if (isTRUE(annotateByDB)){
    if (is.null(x=output)){
      output_file_extended <- paste(workDir,
                                    sub(pattern=".vcf",
                                        replacement=", extended annotation.txt",
                                        x=vcfFile),
                                    sep="/")
    }else{
      output_file_extended <- paste(workDir, paste(output, "extended annotation", "txt", sep="."), sep="/")
    }
  }


  write.table(x=out.anno[[1]],
              file=output_file,
              sep="\t",
              quote=FALSE,
              col.names=TRUE,
              row.names=FALSE)
  if (isTRUE(annotateByDB)){
    write.table(x=out.anno[[2]],
                file=output_file_extended,
                sep="\t",
                quote=FALSE,
                col.names=TRUE,
                row.names=FALSE)
  }
  return(out.anno)
}

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
#' @importFrom rtracklayer import export
#' @importFrom txdbmaker makeTxDbFromGFF
#' @importFrom AnnotationDbi saveDb
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom IRanges disjoin
#' @importFrom S4Vectors queryHits subjectHits mcols mcols<- split
#' @importFrom methods as
developAnnotations <- function(gtfDir=NULL,
                               gtfFile,
                               orf,
                               gtf,
                               src,
                               organism,
                               sqlite,
                               workDir=NULL){
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
#' @importFrom BSgenome getSeq
#' @importFrom GenomicRanges makeGRangesFromDataFrame seqnames start end strand
#' @importFrom methods is
extractSeqSSs <- function(x,
                          genome="BSgenome.Hsapiens.UCSC.hg38",
                          workDir=NULL){
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

#' Extract genomic coordinates of splice sites from exon-exon junctions
#' @description Extraction of the genomic coordinates of splice sites using
#'     experimentally detected exon-exon junctions.
#' @param x character string specifying the name of the tab-delimited TXT file
#'     with experimentally detected exon-exon junctions. This file must contain
#'     the following fields:
#'     i) eej_id     - exon-exon junction ID;
#'     ii) gene_id   - gene ID;
#'     iii) seqnames - name of chromosome or scaffold with prefix "chr";
#'     iv) start     - start coordinate of exon-exon junction (genomic
#'                     position of the last nucleotide in the upstream exon);
#'     v) end        - end coordinate of exon-exon junction(genomic
#'                     position of the first nucleotide in the downstream exon);
#'     vi) strand    - strand information about exon-exon junction location;
#'     vii-...) ...  - (optionally) one or more samples with summarized raw
#'                     RNA-Seq reads spanning exon-exon junction.
#' @param eejDir character string specifying the name of the directory
#'     containing the tab-delimited TXT file(-s) with experimentally detected
#'     exon-exon junctions. NULL by default that mean the current working
#'     directory.
#' @param thr an integer argument. It is a threshold for minimal coverage
#'     (sequencing depth) of exon-exon junction (in raw reads). Default value
#'     is 5 but can be specified by the user.
#' @param workDir character string specifying the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return object of class list containing two GRanges objects with genomic
#'     coordinates of 5' and 3' splice sites.
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @examples
#' res <- extractSSsCoordsEEJs(x="EEJs.Kasumi_1.RR_KD.PMC_BSU.txt",
#'                             eejDir="Files_EEJs",
#'                             workDir="D:/Vasily Grinev")
#' @export
#' @importFrom GenomicRanges makeGRangesFromDataFrame start end strand start<- end<-
#' @importFrom methods is
extractSSsCoordsEEJs <- function(x, eejDir=NULL, thr=5, workDir=NULL){
  ### Setting the working directory.
  if (is.null(x=workDir)){
    workDir <- getwd()
  }
  ### Full path to the file.
  if (is.null(x=eejDir)){
    path <- workDir
  }else{
    path <- paste(workDir, eejDir, sep="/")
  }
  ### Retrieving and validation of the file extension.
  frt <- tools::file_ext(x=x)
  if (!frt %in% "txt"){
    stop("Invalid file format")
  }
  ### Loading of exon-exon junctions as an object of class GRanges.
  EEJs <- read.table(file=paste(path, x, sep="/"),
                     sep="\t",
                     header=TRUE,
                     quote="\"",
                     as.is=TRUE)
  if (ncol(x=EEJs) > 7 & !is.null(x=thr)){
    EEJs <- EEJs[rowSums(x=EEJs[, -1:-7] >= thr) >= 1, ]
  }
  rownames(x=EEJs) <- NULL
  EEJs <- makeGRangesFromDataFrame(df=EEJs,
                                   keep.extra.columns=TRUE)
  ### Converting of genomic coordinates of exon-exon junctions
  #   in genomic coordinates of splice sites.
  fiveSSs <- EEJs
  start(x=fiveSSs[strand(x=fiveSSs) == "+", ]) <-
    start(x=fiveSSs[strand(x=fiveSSs) == "+", ]) - 2
  end(x=fiveSSs[strand(x=fiveSSs) == "+", ]) <-
    start(x=fiveSSs[strand(x=fiveSSs) == "+", ]) + 8
  start(x=fiveSSs[strand(x=fiveSSs) == "-", ]) <-
    end(x=fiveSSs[strand(x=fiveSSs) == "-", ]) - 6
  end(x=fiveSSs[strand(x=fiveSSs) == "-", ]) <-
    start(x=fiveSSs[strand(x=fiveSSs) == "-", ]) + 8
  threeSSs <- EEJs
  start(x=threeSSs[strand(x=threeSSs) == "+", ]) <-
    end(x=threeSSs[strand(x=threeSSs) == "+", ]) - 20
  end(x=threeSSs[strand(x=threeSSs) == "+", ]) <-
    start(x=threeSSs[strand(x=threeSSs) == "+", ]) + 22
  start(x=threeSSs[strand(x=threeSSs) == "-", ]) <-
    start(x=threeSSs[strand(x=threeSSs) == "-", ]) - 2
  end(x=threeSSs[strand(x=threeSSs) == "-", ]) <-
    start(x=threeSSs[strand(x=threeSSs) == "-", ]) + 22
  SSs <- list(fiveSSs=fiveSSs, threeSSs=threeSSs)
  ### Returning the final object.
  return(SSs)
}


#' Extract genomic coordinates of splice sites from annotations
#' @description Extraction of the genomic coordinates of splice sites using
#'     genomic/transcriptomic annotations.
#' @param x character string specifying the name of file with input data.
#'     Accepted formats: GFF2.5/GTF or GFF3 (the gz archive of GFF2.5/GTF or
#'     GFF3 file is also accepted). The file should include models of genes in
#'     general/generic feature format with mandatory field "transcript_id".
#' @param gtfDir character string specifying the name of the directory
#'     containing the file with input data. NULL by default that mean the
#'     current working directory.
#' @param workDir character string specifying the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return object of class list containing two GRanges objects with genomic
#'     coordinates of 5' and 3' splice sites.
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @examples
#' res <- extractSSsCoordsGTF(x="Ensembl_GRCh38.p14_release.114.gtf",
#'                            gtfDir="Files_GTF",
#'                            workDir="D:/Vasily Grinev")
#' @export
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges makeGRangesListFromDataFrame makeGRangesFromDataFrame
#' @importFrom GenomicRanges start end strand start<- end<- elementMetadata
#' @importFrom IRanges setdiff
#' @importFrom S4Vectors colnames
#' @importFrom methods as
extractSSsCoordsGTF <- function(x, gtfDir=NULL, workDir=NULL){
  ### Setting the working directory.
  if (is.null(x=workDir)){
    workDir <- getwd()
  }
  ### Full path to the file.
  if (is.null(x=gtfDir)){
    path <- workDir
  }else{
    path <- paste(workDir, gtfDir, sep="/")
  }
  ### Retrieving and validation of the file extension.
  frt <- tools::file_ext(x=x)
  if (!frt %in% c("gtf", "gff3", "gz")){
    stop("Invalid file format")
  }
  ### Loading of exon-exon junctions as an object of class GRanges.
  EEJs <- import(con=paste(path, x, sep="/"))
  EEJs <- EEJs[gsub(pattern=".+exon",
                    replacement="exon",
                    x=EEJs$type) == "exon", ]
  if ("transcript_id" %in% colnames(x=elementMetadata(x=EEJs)) == "TRUE"){
    splitField <- "transcript_id"
  }else{
    splitField <- "Name"
  }
  EEJs <- makeGRangesListFromDataFrame(df=as.data.frame(x=EEJs),
                                       split.field=splitField,
                                       keep.extra.columns=FALSE)
  EEJs <- unlist(x=setdiff(x=range(x=EEJs), EEJs))
  start(x=EEJs) <- start(x=EEJs) - 1
  end(x=EEJs) <- end(x=EEJs) + 1
  EEJs$transcript_id <- names(x=EEJs)
  names(x=EEJs) <- NULL
  EEJs <- as.data.frame(EEJs, stringsAsFactors=FALSE)
  EEJs$seqnames <- as.character(x=EEJs$seqnames)
  if (length(x=unique(x=nchar(EEJs$seqnames) == 1)) == 2){
    EEJs$seqnames <- paste("chr", EEJs$seqnames, sep="")
    if ("chrMT" %in% EEJs$seqnames == "TRUE"){
      EEJs[EEJs$seqnames == "chrMT", ]$seqnames <- "chrM"
    }
    EEJs <- EEJs[EEJs$seqnames %in%
                   c(paste("chr", c(1:22, "X", "Y", "M"), sep="")), ]
  }else{
    if ("chrMT" %in% EEJs$seqnames == "TRUE"){
      EEJs[EEJs$seqnames == "chrMT", ]$seqnames <- "chrM"
    }
    EEJs <- EEJs[EEJs$seqnames %in%
                   c(paste("chr", c(1:22, "X", "Y", "M"), sep="")), ]
  }
  EEJs$strand <- as.character(x=EEJs$strand)
  EEJs <- EEJs[!duplicated(x=EEJs), ]
  rownames(x=EEJs) <- NULL
  EEJs <- makeGRangesFromDataFrame(df=EEJs,
                                   keep.extra.columns=TRUE)
  ### Converting of genomic coordinates of exon-exon junctions
  #   in genomic coordinates of splice sites.
  fiveSSs <- EEJs
  start(x=fiveSSs[strand(x=fiveSSs) == "+", ]) <-
    start(x=fiveSSs[strand(x=fiveSSs) == "+", ]) - 2
  end(x=fiveSSs[strand(x=fiveSSs) == "+", ]) <-
    start(x=fiveSSs[strand(x=fiveSSs) == "+", ]) + 8
  start(x=fiveSSs[strand(x=fiveSSs) == "-", ]) <-
    end(x=fiveSSs[strand(x=fiveSSs) == "-", ]) - 6
  end(x=fiveSSs[strand(x=fiveSSs) == "-", ]) <-
    start(x=fiveSSs[strand(x=fiveSSs) == "-", ]) + 8
  threeSSs <- EEJs
  start(x=threeSSs[strand(x=threeSSs) == "+", ]) <-
    end(x=threeSSs[strand(x=threeSSs) == "+", ]) - 20
  end(x=threeSSs[strand(x=threeSSs) == "+", ]) <-
    start(x=threeSSs[strand(x=threeSSs) == "+", ]) + 22
  start(x=threeSSs[strand(x=threeSSs) == "-", ]) <-
    start(x=threeSSs[strand(x=threeSSs) == "-", ]) - 2
  end(x=threeSSs[strand(x=threeSSs) == "-", ]) <-
    start(x=threeSSs[strand(x=threeSSs) == "-", ]) + 22
  SSs <- list(fiveSSs=sort(x=fiveSSs), threeSSs=sort(x=threeSSs))
  ### Returning the final object.
  return(SSs)
}


#' Predict amino acid coding changes for single nucleotide variants
#' @description Prediction of amino acid coding changes for single nucleotide
#'     variants.
#' @param vcfDir character string giving the name (or path to and name) of
#'     directory with VCF file(-s). NULL by default, which means the current
#'     working directory will be used instead of the specified VCF directory.
#' @param vcfFile character string giving the name of input VCF file.
#' @param OrgDb character string giving the name of the OrgDb data package with
#'     genome wide annotations. Default value is package org.Hs.eg.db.
#' @param TxDb character string giving the name of the TxDb data package with
#'     transcriptome wide annotations. Default value is package
#'     TxDb.Hsapiens.UCSC.hg38.knownGene.
#' @param genome character string giving the name of the BSgenome data package
#'     to be used for reference sequence extraction. Default value is
#'     package BSgenome.Hsapiens.UCSC.hg38.
#' @param proteins character string giving the name of TXT file in tab-delimited
#'     format with protein coding annotations of Ensembl transcripts. This file
#'     must contains the following IDs: Ensembl transcript IDs, Ensembl gene
#'     IDs, Ensembl protein _IDs, PDB IDs, UniProtKB/Swiss-Prot IDs and
#'     UniProtKB/TrEMBL IDs.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return object of class data frame with predicted and annotated amino acid
#'     coding changes for single nucleotide variants.
#' @author Vasily V. Grinev.
#' @examples
#' vcfFile="test_seq5, fisherCalleR, SNVs.vcf"
#' proteins="Ensembl_GRCh38.p14_release.114_proteins.txt"
#' intrSNVtoAAs <- introduceSNVtoAAs(vcfDir="Files_VCF",
#'                                   vcfFile=vcfFile,
#'                                   OrgDb="org.Hs.eg.db",
#'                                   TxDb="TxDb.Hsapiens.UCSC.hg38.knownGene",
#'                                   genome="BSgenome.Hsapiens.UCSC.hg38",
#'                                   proteins=proteins,
#'                                   workDir="D:/Vasily Grinev/SNVs")
#' @export
#' @importFrom VariantAnnotation readVcf alt alt<- predictCoding
#' @importFrom SparseArray rowRanges
#' @importFrom Biostrings DNAStringSetList
#' @importFrom dplyr select
#' @importFrom GenomicRanges start
#' @importFrom methods setClassUnion
introduceSNVtoAAs <- function(vcfDir=NULL,
                              vcfFile,
                              OrgDb="org.Hs.eg.db",
                              TxDb="TxDb.Hsapiens.UCSC.hg38.knownGene",
                              genome="BSgenome.Hsapiens.UCSC.hg38",
                              proteins,
                              workDir=NULL){
  ### Loading of required packages.
  suppressMessages(expr=library(package=genome, character.only=TRUE))
  suppressMessages(expr=library(package=OrgDb, character.only=TRUE))
  suppressMessages(expr=library(package=TxDb, character.only=TRUE))
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
  ### Loading of the single nucleotide variations as an object of class
  #   GRanges.
  setClassUnion("ExpData", c("matrix", "SummarizedExperiment"))
  vcf <- readVcf(file=paste(vcfDir, vcfFile, sep="/"))
  ALT <- as.character(x=unlist(x=rowRanges(x=vcf)$ALT))
  vcf <- vcf[rep(x=1:length(x=vcf), times=lengths(x=rowRanges(x=vcf)$ALT)), ]
  alt(x=vcf) <- DNAStringSetList(strsplit(x=ALT, split=""))
  ### Loading gene models as an object of class TxDb.
  txDb <- get(x=TxDb)
  ### Loading genome wide annotations as an object of class OrgDb.
  orgDb <- get(x=OrgDb)
  ### Loading of the protein annotations as an object of class data frame.
  pr <- read.table(file=paste(workDir, proteins, sep="/"),
                   sep="\t",
                   header=TRUE,
                   quote="\"",
                   as.is=TRUE)
  rownames(x=pr) <- pr[, 1]
  ### Prediction of amino acid coding changes for single nucleotide variants.
  AAs <- predictCoding(query=vcf, subject=txDb, seqSource=get(x=genome))
  AAs <- AAs[AAs$CONSEQUENCE == "nonsynonymous", ]
  ### Annotation.
  AAs$snv_id <- names(x=AAs)
  AAs$snv_ref <- as.character(x=AAs$REF)
  AAs$snv_alt <- as.character(x=unlist(x=AAs$ALT))
  AAs$snv_qual <- AAs$QUAL
  AAs$Ensembl_gene_IDs <- select(x=orgDb,
                                 keys=AAs$GENEID,
                                 columns="ENSEMBL",
                                 keytype="ENTREZID")$ENSEMBL
  AAs$Entrez_gene_IDs <- AAs$GENEID
  AAs$HUGO_gene_symbols <- select(x=orgDb,
                                  keys=AAs$GENEID,
                                  columns="SYMBOL",
                                  keytype="ENTREZID")$SYMBOL
  AAs$Ensembl_tx_names <- select(x=txDb,
                                 keys=AAs$TXID,
                                 columns="TXNAME",
                                 keytype="TXID")$TXNAME
  AAs$Ensembl_tx_names <- substr(x=AAs$Ensembl_tx_names, start=1, stop=15)
  AAs$SNV_CDS_loc <- start(x=AAs$CDSLOC)
  AAs$ref_codon <- as.vector(x=as.character(x=AAs$REFCODON))
  AAs$alt_codon <- as.vector(x=as.character(x=AAs$VARCODON))
  AAs$protein_id <- ""
  AAs$pdb_id <- ""
  AAs$uniprot.swissprot_id <- ""
  AAs$uniprot.trembl_id <- ""
  AAs$SNV_prot_loc <- as.vector(x=unlist(x=AAs$PROTEINLOC))
  AAs$ref_aa <- as.vector(x=as.character(x=AAs$REFAA))
  AAs$alt_aa <- as.vector(x=as.character(x=AAs$VARAA))
  AAs <- data.frame(AAs)[, c(25:28, 1:5, 29:42)]
  AAs$seqnames <- as.character(x=AAs$seqnames)
  AAs$strand <- as.character(x=AAs$strand)
  idx <- AAs[, 13] %in% pr[, 1]
  AAs[idx, 17:20] <- pr[AAs[idx, 13], 3:6]
  ### Returning the final object.
  return(AAs)
}


#' Introduce single nucleotide variant(-s) to splice site(-s)
#' @description Introduction of single nucleotide variant(-s) into the splice
#'     site sequence(-s).
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
#' @param vcfDir character string giving the name (or path to and name) of
#'     directory with VCF file(-s). NULL by default, which means the current
#'     working directory will be used instead of the specified VCF directory.
#' @param vcfFile character string giving the name of input VCF file.
#' @param genome character string giving the name of the BSgenome data package
#'     to be used for reference sequence extraction. Default value is
#'     package BSgenome.Hsapiens.UCSC.hg38.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return object of class list containing extended data frame with splice
#'     sites and DNAStringSet with corrected splice site sequences.
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @examples
#' coords <- "EEJs.Kasumi_1.RR_KD.PMC_BSU.fiveSSs.txt"
#' vcf <- "test_seq4, fisherCalleR, SNVs.filtered.vcf"
#' intrSNVtoSSs <- introduceSNVtoSSs(x=coords,
#'                                   vcfDir="Files_VCF",
#'                                   vcfFile=vcf,
#'                                   genome="BSgenome.Hsapiens.UCSC.hg38",
#'                                   workDir="D:/Vasily Grinev")
#' @export
#' @importFrom VariantAnnotation readVcf
#' @importFrom SparseArray rowRanges
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom GenomicRanges start end strand seqnames
#' @importFrom Biostrings DNAStringSet reverseComplement
#' @importFrom methods setClassUnion
introduceSNVtoSSs <- function(x,
                              vcfDir=NULL,
                              vcfFile,
                              genome="BSgenome.Hsapiens.UCSC.hg38",
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
  ### Loading of the single nucleotide variations as an object of class
  #   GRanges.
  setClassUnion("ExpData", c("matrix", "SummarizedExperiment"))
  vcf <- rowRanges(x=readVcf(file=paste(vcfDir, vcfFile, sep="/")))[, "ALT"]
  ALT <- as.character(x=unlist(x=vcf$ALT))
  vcf <- vcf[rep(x=1:length(x=vcf), times=lengths(x=vcf$ALT)), ]
  vcf$ALT <- ALT
  ### Loading of the genomic coordinates of splice sites as an object of
  #   class GRanges.
  if (is.object(x=x)){
    spl.sites <- x
  }else{
    frt <- tools::file_ext(x=x)
    if (!frt %in% "txt"){
      stop("Invalid file format")
    }
    spl.sites <- read.table(file=paste(workDir, x, sep="/"),
                            sep="\t",
                            header=TRUE,
                            quote="\"",
                            as.is=TRUE)
    spl.sites <- makeGRangesFromDataFrame(df=spl.sites,
                                          keep.extra.columns=TRUE)
  }
  ### Find overlaps between splice sites and SNVs.
  hits <- findOverlaps(query=spl.sites,
                       subject=vcf,
                       type="any",
                       ignore.strand=TRUE)
  if (length(x=hits) > 0){
    spl.sites <- spl.sites[queryHits(x=hits)]
    vcf <- vcf[subjectHits(x=hits)]
    ### Extraction of reference sequences of splice sites to be corrected.
    spl.sites$seq_original <- extractSeqSSs(x=spl.sites,
                                            genome=genome,
                                            workDir=workDir)
    spl.sites$seq_original <- as.character(x=spl.sites$seq_original)
    spl.sites$seq_corrected <- spl.sites$seq_original
    rev_compl <- as.character(x=strand(x=spl.sites)) == "-"
    spl.sites$seq_corrected[rev_compl] <-
      reverseComplement(x=DNAStringSet(x=spl.sites$seq_corrected[rev_compl]))
    ### Correction of the splice site sequences.
    if (length(x=spl.sites) != length(x=vcf)){
      stop("The list of splice sites to be corrected and the list of ",
           "relevant SNVs differ in length.")
    }
    substr(x=spl.sites$seq_corrected,
           start=start(x=vcf) - start(x=spl.sites) + 1,
           stop=start(x=vcf) - start(x=spl.sites) + 1) <- vcf$ALT
    spl.sites$seq_corrected[rev_compl] <-
      reverseComplement(x=DNAStringSet(x=spl.sites$seq_corrected[rev_compl]))
    ### Returning the final object.
    anno <- data.frame(spl.sites)
    anno <- anno[, c(1:7,
                     (ncol(x=anno) - 1):ncol(x=anno),
                     8:(ncol(x=anno) - 2))]
    seqSSs <- DNAStringSet(x=anno$seq_corrected)
    names(x=seqSSs) <- paste(paste(seqnames(x=spl.sites),
                                   paste(start(x=spl.sites),
                                         end(x=spl.sites),
                                         sep="-"),
                                   sep=":"),
                             strand(x=spl.sites),
                             sep="_str")
    return(list(annotation=anno, sequences=seqSSs))
  }else{
    warning("WARNING: No matches were found between splice sites and SNVs.")
  }
}


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
#' @importFrom Biostrings readDNAStringSet width uniqueLetters
#' @importFrom methods is
scoreSSs <- function(x,
                     type,
                     workDir=NULL){
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ### Connection to splice site scoring models.
  models <- system.file("extdata", "me2x5", package="huGSVcalleR")
  setwd(dir=dirname(path=models))
  # models <- "D:/Vasily Grinev/exdata/me2x5"
  # setwd(dir=dirname(path=models))
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

