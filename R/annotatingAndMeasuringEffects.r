#' Annotation of single nucleotide variations by dbSNP
#' @description This function compares basic annotation of variants,
#'     produced by annotateBasic() function, with annotation from dbSNVs
#'     current release. Current version of function was successfully
#'     tested on GCF_000001405.40 release of dbSNP. List of extracted
#'     dbSNP fields is fixed for now and cannot be changed.
#' @param annoDir character string giving the name (or path to and name)
#'     of directory with TXT file in tab-delimited format with annotated
#'     single nucleotide variations. annoDir is NULL by default, which
#'     means the current working directory will be used instead of the
#'     specified annotation directory.
#' @param annoFile character string giving the name of input TXT
#'     tab-delimited file. By default, using of annotateBasic() output
#'     file is assumed. File must contains next columns:
#'     i) seqnames - name of chromosome or scaffold;
#'     ii) start - start coordinate of variation;
#'     iii) end - end coordinate of variation;
#'     iv) strand - strand of DNA;
#'     iv) ref - reference allele;
#'     v) alt - alt allele, called in this position.
#'     All extra columns wont be used in comparison and will be added in
#'     output file without changes.
#' @param output character string giving the name (without extension) of
#'     TXT file in tab-delimited format for storing annotation. The
#'     postfix "dbSNP annotation" will be added to this name. NULL by
#'     default, which means that the name will be generated
#'     automatically based on the name of the input VCF file.
#' @param dbDir character string giving the path to directory with VCF
#'     file of dbSNP release. NULL by default, which means the current
#'     working directory will be used instead of the specified directory.
#'     Can be same with annoDir.
#' @param dbFile character string giving the name of dbSNP VCF file for
#'     extended annotation. "GCF_000001405.40.gz" by default.
#' @param genome character string giving the version of genome. "hg38" by
#'     default. Used for annotaded SNV subsetting from dbSNP.
#' @param chrToNC logical. If TRUE, sequence names in basic annotation
#'     will be converted from UCSC-style format (chr) to NCBI-style
#'     (NC_***). Requires an file with possible sequence names, provided
#'     in chrID. At the end of annotation, UCSC format of names will be
#'     automatically restored. FALSE by default.
#' @param chrID character string giving the name tab-delimited TXT file
#'     with names of chrmosome sequences in NCBI-RefSeq format (first
#'     column) and UCSC format (second column). File must be placed into
#'     dbDir.
#' @param workDir character string giving the path to and name of working
#'     directory. NULL by default that means the current working
#'     directory.
#' @return an tab-delimited TXT file with basic SNVs annotation and
#'     additional fields for variations that were found into dbSNP (for
#'     others columns will contain NA).
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @examples
#' ex <- system.file("extdata", package = "huGSVcalleR")
#' workDir <- tempdir()
#' file.copy(
#'     file.path(ex, "test_annotation.vcf"),
#'     file.path(workDir, "test.vcf"),
#'     overwrite = TRUE
#' )
#' file.copy(file.path(ex, "ref_genes.txt"), workDir, overwrite = TRUE)
#' annotateBasic(
#'     vcfDir = NULL,
#'     vcfFile = "test.vcf",
#'     ref_genes = "ref_genes.txt",
#'     output = "example",
#'     workDir = workDir
#' )
#' anno <- read.delim(
#'     file.path(workDir, "example.basic annotation.txt"),
#'     sep = "\t"
#' )
#' snv <- anno[1, ]
#' fields <- c(
#'     RS = "Integer", dbSNPBuildID = "Integer", SAO = "Integer",
#'     PM = "Integer", FREQ = "String", COMMON = "Integer",
#'     CLNHGVS = "String", CLNVI = "String", CLNORIGIN = "String",
#'     CLNSIG = "String", CLNDISDB = "String", CLNDN = "String",
#'     CLNREVSTAT = "String", CLNACC = "String"
#' )
#' header <- c(
#'     "##fileformat=VCFv4.2",
#'     "##contig=<ID=chr1,length=1000000>",
#'     paste0(
#'         "##INFO=<ID=", names(fields),
#'         ",Number=1,Type=", fields,
#'         ",Description=\"x\">"
#'     ),
#'     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
#' )
#' info <- paste(
#'     "RS=1", "dbSNPBuildID=150", "SAO=0", "PM=0", "FREQ=0.01",
#'     "COMMON=1", "CLNHGVS=.", "CLNVI=.", "CLNORIGIN=.", "CLNSIG=.",
#'     "CLNDISDB=.", "CLNDN=.", "CLNREVSTAT=.", "CLNACC=.",
#'     sep = ";"
#' )
#' record <- paste(
#'     snv$seqnames, snv$start, ".", snv$ref, snv$alt,
#'     ".", "PASS", info,
#'     sep = "\t"
#' )
#' writeLines(c(header, record), file.path(workDir, "dbsnp.vcf"))
#' gz <- Rsamtools::bgzip(
#'     file.path(workDir, "dbsnp.vcf"),
#'     dest = file.path(workDir, "dbsnp.vcf.gz"),
#'     overwrite = TRUE
#' )
#' Rsamtools::indexTabix(gz, format = "vcf")
#' annotateDbSNP(
#'     annoDir = NULL,
#'     annoFile = "example.basic annotation.txt",
#'     dbDir = NULL,
#'     dbFile = "dbsnp.vcf.gz",
#'     chrToNC = FALSE,
#'     chrID = "unused.txt",
#'     workDir = workDir
#' )
#' @export
#' @importFrom data.table data.table
#' @importFrom VariantAnnotation readVcf info ScanVcfParam fixed
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#'     GRanges
#' @importFrom IRanges IRanges
annotateDbSNP <- function(
        annoDir = NULL,
        annoFile,
        output = NULL,
        dbDir = NULL,
        dbFile = "GCF_000001405.40.gz",
        genome = "hg38",
        chrToNC = FALSE,
        chrID,
        workDir = NULL) {
    ### Resolving the annotation and database directories.
    paths <- .resolveAnnotationPaths(annoDir, dbDir, workDir)
    ### Loading of basic SNVs annotation.
    SNVs <- read.table(
        file = paste(paths$annoDir, annoFile, sep = "/"),
        header = TRUE,
        sep = "\t"
    )
    ### Conversion chr -> NCBI chromosome ID (if requested).
    chr.backup <- NULL
    if (isTRUE(chrToNC)) {
        chr.backup <- SNVs$seqnames
        SNVs <- .applyChrToNC(SNVs, paths$dbDir, chrID)
    }
    ### Query ranges for dbSNP subsetting and overlap search.
    qry <- .buildSnvQueryRanges(SNVs)
    ### Loading of dbSNP annotation and finding overlaps.
    fields <- .dbSnpInfoFields()
    ref.snp <- .loadDbSnpVcf(paths$dbDir, dbFile, genome, qry, fields)
    ### Pair-by-pair comparison of ALT alleles.
    hits.ref <- .findDbSnpOverlaps(qry, ref.snp)
    alt.hits <- .compareAltAlleles(SNVs, ref.snp, hits.ref)
    ### Extracting dbSNP fields and building the final data.frame.
    hits.to <- .extractDbSNPFields(ref.snp, hits.ref, alt.hits, fields)
    anno <- .buildDbSNPAnnoFrame(SNVs, hits.to, hits.ref, alt.hits)
    ### Restoring chromosome names.
    if (isTRUE(chrToNC)) {
        anno$seqnames <- chr.backup
    }
    ### Writing the output file.
    .writeDbSnpAnnotation(anno, workDir, annoFile, output)
}

#' Resolve annotation and database directories
#' @description Internal helper for annotateDbSNP(). Returns the absolute
#'     paths to the annotation and database directories.
#' @param annoDir character string or NULL.
#' @param dbDir character string or NULL.
#' @param workDir character string or NULL.
#' @return list with annoDir and dbDir.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.resolveAnnotationPaths <- function(annoDir, dbDir, workDir) {
    if (is.null(x = workDir)) {
        workDir <- getwd()
    }
    list(
        workDir = workDir,
        annoDir = if (is.null(x = annoDir)) {
            workDir
        } else {
            paste(workDir, annoDir, sep = "/")
        },
        dbDir = if (is.null(x = dbDir)) {
            workDir
        } else {
            paste(workDir, dbDir, sep = "/")
        }
    )
}

#' Apply UCSC-to-NCBI chromosome name conversion
#' @description Internal helper for annotateDbSNP(). Reads the
#'     correspondence table and converts sequence names.
#' @param SNVs data.frame with basic SNVs annotation.
#' @param dbDir character string with the database directory.
#' @param chrID character string with the correspondence table name.
#' @return data.frame with converted sequence names.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.applyChrToNC <- function(SNVs, dbDir, chrID) {
    NC <- read.table(
        file = paste(dbDir, chrID, sep = "/"),
        header = TRUE,
        sep = "\t"
    )
    .convertChrToNCBI(SNVs, NC)
}

#' Build GRanges with SNV query ranges
#' @description Internal helper for annotateDbSNP(). Builds the query
#'     ranges reused for VCF subsetting and overlap search.
#' @param SNVs data.frame with basic SNVs annotation.
#' @return GRanges object with query ranges.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.buildSnvQueryRanges <- function(SNVs) {
    GRanges(
        seqnames = SNVs$seqnames,
        ranges = IRanges(start = SNVs$start, end = SNVs$end),
        strand = SNVs$strand
    )
}

#' Return the dbSNP INFO fields used by annotateDbSNP()
#' @description Internal helper for annotateDbSNP(). Returns the fixed
#'     list of extracted dbSNP fields.
#' @return character vector with field names.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.dbSnpInfoFields <- function() {
    c(
        "RS", "dbSNPBuildID", "SAO",
        "PM", "FREQ", "COMMON", "CLNHGVS", "CLNVI",
        "CLNORIGIN", "CLNSIG", "CLNDISDB", "CLNDN",
        "CLNREVSTAT", "CLNACC"
    )
}

#' Load the dbSNP VCF subset
#' @description Internal helper for annotateDbSNP(). Reads the dbSNP VCF
#'     restricted to the query ranges and INFO fields.
#' @param dbDir character string with the database directory.
#' @param dbFile character string with the dbSNP VCF file name.
#' @param genome character string with the genome version.
#' @param qry GRanges object with query ranges.
#' @param fields character vector with INFO fields.
#' @return CollapsedVCF object with dbSNP annotation.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.loadDbSnpVcf <- function(dbDir, dbFile, genome, qry, fields) {
    vcf <- readVcf(
        file = paste(dbDir, dbFile, sep = "/"),
        genome = genome,
        param = ScanVcfParam(info = fields, which = qry)
    )
    gc()
    vcf
}

#' Find overlaps between SNVs and dbSNP
#' @description Internal helper for annotateDbSNP(). Wraps findOverlaps()
#'     with type = "equal".
#' @param qry GRanges object with SNV query ranges.
#' @param ref.snp CollapsedVCF object with dbSNP annotation.
#' @return Hits object.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.findDbSnpOverlaps <- function(qry, ref.snp) {
    findOverlaps(
        query = qry,
        subject = rowRanges(ref.snp),
        type = "equal"
    )
}

#' Compare ALT alleles between SNVs and dbSNP
#' @description Internal helper for annotateDbSNP(). Returns a logical
#'     vector indicating which overlaps have matching ALT alleles.
#' @param SNVs data.frame with basic SNVs annotation.
#' @param ref.snp CollapsedVCF object with dbSNP annotation.
#' @param hits.ref Hits object from findOverlaps().
#' @return logical vector.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.compareAltAlleles <- function(SNVs, ref.snp, hits.ref) {
    alt.compare <- list(
        experimental = strsplit(SNVs$alt[hits.ref@from], ","),
        reference = fixed(ref.snp)[hits.ref@to, c("ALT")]
    )
    vapply(
        seq_along(alt.compare[[1]]),
        function(i) {
            any(
                as.character(alt.compare[[1]][[i]]) %in%
                    as.character(alt.compare[[2]][[i]])
            )
        },
        logical(1)
    )
}

#' Write dbSNP annotation table
#' @description Internal helper for annotateDbSNP(). Resolves the output
#'     file name and writes the annotation table.
#' @param anno data.frame with merged annotation.
#' @param workDir character string with the working directory.
#' @param annoFile character string with the input annotation file name.
#' @param output character string or NULL.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.writeDbSnpAnnotation <- function(anno, workDir, annoFile, output) {
    outputFile <- .resolveDbSNPOutputFile(workDir, annoFile, output)
    write.table(
        x = anno,
        file = outputFile,
        sep = "\t",
        quote = FALSE,
        col.names = TRUE,
        row.names = FALSE
    )
    invisible(NULL)
}

#' Convert chromosome names from UCSC to NCBI style
#' @description Internal helper for annotateDbSNP(). Converts UCSC-style
#'     sequence names (chr...) to NCBI-RefSeq style (NC_...) using the
#'     correspondence table provided by the user.
#' @param SNVs data.frame with basic SNVs annotation.
#' @param NC data.frame with two columns: NCBI-RefSeq accession (first)
#'     and UCSC name (second).
#' @return data.frame with converted sequence names.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.convertChrToNCBI <- function(SNVs, NC) {
    for (i in seq_len(nrow(NC))) {
        SNVs$seqnames[SNVs$seqnames == NC[i, 2]] <- NC[i, 1]
    }
    SNVs
}

#' Extract dbSNP fields for matched positions
#' @description Internal helper for annotateDbSNP(). Extracts the REF/ALT
#'     alleles and the selected INFO fields from the dbSNP VCF for the
#'     positions that matched the experimental SNVs.
#' @param ref.snp CollapsedVCF object with dbSNP annotation.
#' @param hits.ref Hits object from findOverlaps().
#' @param alt.hits logical vector indicating which overlaps have matching
#'     ALT alleles.
#' @param fields character vector with dbSNP INFO fields to extract.
#' @return data.frame with dbSNP annotation for matched positions.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.extractDbSNPFields <- function(ref.snp, hits.ref, alt.hits, fields) {
    refalt <- data.table(
        unlist(x = lapply(
            X = rowRanges(ref.snp)$REF[hits.ref@to[which(alt.hits)]],
            paste0,
            collapse = ","
        )),
        unlist(x = lapply(
            X = rowRanges(ref.snp)$ALT[hits.ref@to[which(alt.hits)]],
            paste0,
            collapse = ","
        ))
    )
    colnames(refalt) <- c("dbSNP.ref", "dbSNP.alt")
    info.db <- info(ref.snp)[hits.ref@to[which(alt.hits)], fields]
    info.db$FREQ <- unlist(lapply(info.db$FREQ, paste, collapse = " "))
    for (f in c(
        "CLNHGVS", "CLNVI", "CLNORIGIN", "CLNSIG",
        "CLNDISDB", "CLNDN", "CLNREVSTAT", "CLNACC"
    )) {
        info.db[[f]] <- unlist(x = lapply(
            X = info.db[[f]],
            paste0,
            collapse = ","
        ))
    }
    as.data.frame(cbind(refalt, as.data.frame(info.db)))
}

#' Build final dbSNP annotation data.frame
#' @description Internal helper for annotateDbSNP(). Merges the basic SNVs
#'     annotation with the dbSNP fields, filling NA for positions without
#'     dbSNP matches.
#' @param SNVs data.frame with basic SNVs annotation.
#' @param hits.to data.frame with dbSNP annotation for matched positions.
#' @param hits.ref Hits object from findOverlaps().
#' @param alt.hits logical vector indicating which overlaps have matching
#'     ALT alleles.
#' @return data.frame with merged annotation.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.buildDbSNPAnnoFrame <- function(SNVs, hits.to, hits.ref, alt.hits) {
    anno <- as.data.frame(matrix(
        data = NA,
        nrow = nrow(SNVs),
        ncol = ncol(hits.to)
    ))
    colnames(anno) <- colnames(hits.to)
    for (i in seq_len(ncol(anno))) {
        class(anno[, i]) <- class(hits.to[, i])
    }
    anno[hits.ref@from[alt.hits], ] <- hits.to
    cbind(SNVs, anno)
}

#' Resolve output file name for annotateDbSNP()
#' @description Internal helper for annotateDbSNP(). Builds the full path
#'     to the output TXT file based on the user-provided output name or
#'     the input annotation file name.
#' @param workDir character string giving the path to working directory.
#' @param annoFile character string giving the name of input annotation
#'     file.
#' @param output character string giving the user-provided output name
#'     (without extension), or NULL.
#' @return character string with the full path to the output file.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.resolveDbSNPOutputFile <- function(workDir, annoFile, output) {
    if (is.null(x = output)) {
        paste(
            workDir,
            sub(
                pattern = ".txt",
                replacement = ", dbSNP annotation.txt",
                x = annoFile
            ),
            sep = "/"
        )
    } else {
        paste(
            workDir,
            paste(output, "dbSNP annotation", "txt", sep = "."),
            sep = "/"
        )
    }
}


#' Annotate of single nucleotide variations
#' @description This function performs a basic annotation of the single
#'     nucleotide variations given in the standard VCF format, and
#'     defines its allegiace to gene.
#' @param vcfDir character string giving the name (or path to and name) of
#'     directory with VCF file(-s). NULL by default, which means the
#'     current working directory will be used instead of the specified
#'     VCF directory.
#' @param vcfFile character string giving the name of input VCF file.
#'     File must not contain NON_REF fields (for HaplotypeCaller it means
#'     run with -ERC NONE).
#' @param ref_genes character string giving the name of TXT file in
#'     tab-delimited format with basic annotations of reference genes.
#'     This file must contains the following fields:
#'     i) gene_id        - (optionally) gene ID (for example, Ensembl
#'                         based);
#'     ii) gene_name     - name of gene (for example, HUGO based);
#'     iii) gene_biotype - biotype of gene;
#'     iv) seqnames      - name of chromosome or scaffold with prefix
#'                         "chr";
#'     v) start          - start coordinate of gene;
#'     vi) end           - end coordinate of gene;
#'     vii) strand       - (optionally) strand information about gene.
#' @param output character string giving the name (without extension) of
#'     TXT file in tab-delimited format for storing annotated single
#'     nucleotide variations. The words "basic" will be added to this
#'     name. NULL by default, which means that the name will be generated
#'     automatically based on the name of the input VCF file.
#' @param workDir character string giving the path to and name of working
#'     directory. NULL by default that means the current working
#'     directory.
#' @return an object of class list with two slots and two new TXT file in
#'     tab-delimited format with basic and extended annotation of
#'     provided single nucleotide variations.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @examples
#' ex <- system.file("extdata", package = "huGSVcalleR")
#' workDir <- tempdir()
#' file.copy(
#'     file.path(ex, "test_annotation.vcf"),
#'     file.path(workDir, "test.vcf"),
#'     overwrite = TRUE
#' )
#' file.copy(file.path(ex, "ref_genes.txt"), workDir, overwrite = TRUE)
#' annotateBasic(
#'     vcfDir = NULL,
#'     vcfFile = "test.vcf",
#'     ref_genes = "ref_genes.txt",
#'     output = "example",
#'     workDir = workDir
#' )
#' @export
#' @importFrom data.table data.table
#' @importFrom VariantAnnotation readVcf info
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom SummarizedExperiment "rowRanges<-"
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
annotateBasic <- function(
        vcfDir = NULL,
        vcfFile,
        ref_genes,
        output = NULL,
        workDir = NULL) {
    ### Setting the working directory.
    if (is.null(x = workDir)) {
        workDir <- getwd()
    }
    ### Setting the VCF directory.
    if (is.null(x = vcfDir)) {
        vcfDir <- workDir
    } else {
        vcfDir <- paste(workDir, vcfDir, sep = "/")
    }
    ### Loading of the single nucleotide variations.
    vcf <- readVcf(file = paste(vcfDir, vcfFile, sep = "/"))
    ### Loading of the reference genes.
    genes <- .loadRefGenes(workDir, ref_genes)
    ### Annotation of the single nucleotide variations with gene names.
    vcf <- .annotateWithGenes(vcf, genes)
    ### Building an object with basic annotation.
    SNVs <- .buildBasicAnnoFrame(vcf)
    ### Writing output file.
    output_file <- .resolveBasicOutputFile(workDir, vcfFile, output)
    write.table(
        x = SNVs,
        file = output_file,
        sep = "\t",
        quote = FALSE,
        col.names = TRUE,
        row.names = FALSE
    )
}

#' Load reference gene annotations
#' @description Internal helper for annotateBasic(). Reads a tab-delimited
#'     TXT file with reference gene annotations and converts it to a
#'     GRanges object.
#' @param workDir character string giving the path to working directory.
#' @param ref_genes character string giving the name of the input file.
#' @return GRanges object with reference gene annotations.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.loadRefGenes <- function(workDir, ref_genes) {
    genes <- read.table(
        file = paste(workDir, ref_genes, sep = "/"),
        sep = "\t",
        header = TRUE,
        quote = "\"",
        as.is = TRUE
    )
    makeGRangesFromDataFrame(df = genes, keep.extra.columns = TRUE)
}

#' Annotate SNVs with gene names
#' @description Internal helper for annotateBasic(). Adds a gene_name
#'     metadata column to the VCF with comma-separated gene names for
#'     each SNV.
#' @param vcf CollapsedVCF object with SNVs.
#' @param genes GRanges object with reference gene annotations.
#' @return CollapsedVCF object with gene_name metadata column.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.annotateWithGenes <- function(vcf, genes) {
    mcols(x = rowRanges(x = vcf))$gene_name <- ""
    hits <- findOverlaps(
        query = rowRanges(x = vcf),
        subject = genes[!is.na(x = genes$gene_name)],
        type = "within",
        ignore.strand = TRUE
    )
    hits <- data.table(cbind(
        hits@from,
        genes[!is.na(x = genes$gene_name)][hits@to, ]$gene_name
    ))
    hits <- hits[, paste(V2, collapse = ","), by = V1]
    colnames(x = hits) <- c("V1", "V2")
    hits$V1 <- as.numeric(x = hits$V1)
    hits$V2 <- vapply(
        hits$V2,
        function(y) {
            paste0(
                sort(x = unique(x = strsplit(x = y, split = ",")[[1]])),
                collapse = ","
            )
        },
        character(1)
    )
    rowRanges(x = vcf)[hits$V1, ]$gene_name <- as.vector(x = hits$V2)
    vcf
}

#' Build basic annotation data.frame
#' @description Internal helper for annotateBasic(). Converts the
#'     annotated CollapsedVCF into a data.frame with the standard columns
#'     required by annotateDbSNP().
#' @param vcf CollapsedVCF object with annotated SNVs.
#' @return data.frame with basic annotation.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.buildBasicAnnoFrame <- function(vcf) {
    SNVs <- cbind(data.frame(rowRanges(x = vcf)), info(x = vcf))
    SNVs$ALT <- unlist(x = lapply(X = SNVs$ALT, paste0, collapse = ","))
    SNVs$snv_id <- paste(
        paste(paste(SNVs$seqnames, SNVs$start, sep = ":"),
            SNVs$REF,
            sep = "_"
        ),
        SNVs$ALT,
        sep = "/"
    )
    SNVs <- SNVs[, c(
        "snv_id", "seqnames", "start", "end", "strand", "REF",
        "ALT", "DP", "QUAL", "gene_name"
    )]
    colnames(x = SNVs) <- c(
        colnames(x = SNVs)[seq_len(5)],
        "ref", "alt", "pos_depth", "score", "gene_name"
    )
    SNVs$seqnames <- as.character(x = SNVs$seqnames)
    SNVs$strand <- as.character(x = SNVs$strand)
    SNVs
}

#' Resolve output file name for annotateBasic()
#' @description Internal helper for annotateBasic(). Builds the full path
#'     to the output TXT file.
#' @param workDir character string giving the path to working directory.
#' @param vcfFile character string giving the name of input VCF file.
#' @param output character string giving the user-provided output name
#'     (without extension), or NULL.
#' @return character string with the full path to the output file.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.resolveBasicOutputFile <- function(workDir, vcfFile, output) {
    if (is.null(x = output)) {
        paste(
            workDir,
            sub(
                pattern = ".vcf.bgz",
                replacement = ", basic annotation.txt",
                x = vcfFile
            ),
            sep = "/"
        )
    } else {
        paste(
            workDir,
            paste(output, "basic annotation", "txt", sep = "."),
            sep = "/"
        )
    }
}


#' Develop Ensembl-like annotations of experimental transcriptome
#' @description This function converts a standard GTF/GFF file with
#'     experimental transcriptome to Ensembl-like annotations.
#' @param gtfDir character string giving the name (or path to and name) of
#'     directory with GTF file(-s). NULL by default, which means the
#'     current working directory will be used instead of the specified
#'     GTF directory.
#' @param gtfFile character string giving the name of input GTF/GFF file
#'     with experimental transcriptome. It is typically StringTie output.
#'     Valid format is "gtf" or "gff".
#' @param orf character string giving the name of tab-delimited TXT file
#'     with coordinates of open reading frames in experimental
#'     transcripts. It is typically ORFhunteR output. This file should
#'     include three fields:
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
#' @param sqlite character string giving the name of output SQLite
#'     database with Ensembl-like annotated experimental transcriptome.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return GTF/GFF file and local SQLite database with Ensembl-like
#'     annotations of the experimental transcriptome.
#' @author Vasily V. Grinev.
#' @examples
#' ex <- system.file("extdata", package = "huGSVcalleR")
#' workDir <- tempdir()
#' file.copy(
#'     file.path(ex, "test_annotation.gtf"),
#'     workDir,
#'     overwrite = TRUE
#' )
#' write.table(
#'     data.frame(
#'         transcript_id = "TX1",
#'         orf_start = 1L,
#'         orf_end = 30L
#'     ),
#'     file = file.path(workDir, "orfs.txt"),
#'     sep = "\t",
#'     quote = FALSE,
#'     row.names = FALSE
#' )
#' if (requireNamespace("GenomeInfoDbData", quietly = TRUE)) {
#' developAnnotations(
#'     gtfDir = NULL,
#'     gtfFile = "test_annotation.gtf",
#'     orf = "orfs.txt",
#'     gtf = "example.anno.gtf",
#'     src = "example",
#'     organism = "Homo sapiens",
#'     sqlite = "example.anno.sqlite",
#'     workDir = workDir
#' )
#' }
#' @export
#' @importFrom rtracklayer import export
#' @importFrom txdbmaker makeTxDbFromGFF
#' @importFrom AnnotationDbi saveDb
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom IRanges disjoin
#' @importFrom S4Vectors queryHits subjectHits mcols mcols<- split
#' @importFrom methods as
developAnnotations <- function(gtfDir = NULL,
    gtfFile,
    orf,
    gtf,
    src,
    organism,
    sqlite,
    workDir = NULL) {
    ### Setting the working directory.
    if (is.null(x = workDir)) {
        workDir <- getwd()
    }
    ### Setting the GTF directory.
    if (is.null(x = gtfDir)) {
        gtfDir <- workDir
    } else {
        gtfDir <- paste(workDir, gtfDir, sep = "/")
    }
    ### Loading of the experimental annotations and ORFs.
    GTFs <- import(con = paste(gtfDir, gtfFile, sep = "/"))
    ORFs <- read.table(
        file = paste(workDir, orf, sep = "/"),
        sep = "\t",
        header = TRUE,
        quote = "\"",
        as.is = TRUE
    )
    ### Development of the Ensembl-like annotations.
    GTFs <- .assignBiotypes(GTFs, ORFs)
    genes <- .buildGeneRanges(GTFs)
    coding <- GTFs[GTFs$transcript_biotype == "protein_coding", ]
    coding <- .processCodingTranscripts(coding, genes, ORFs)
    noncoding <- GTFs[GTFs$transcript_biotype == "noncoding", ]
    noncoding <- .processNoncodingTranscripts(noncoding, genes)
    GTFs <- .combineTranscripts(coding, noncoding)
    ### Export of the developed annotations as file in format GTF/GFF.
    export(object = GTFs, con = paste(gtfDir, gtf, sep = "/"))
    ### Export of the developed annotations as SQLite database.
    TxDb <- makeTxDbFromGFF(
        file = paste(gtfDir, gtf, sep = "/"),
        format = "auto",
        dataSource = src,
        organism = organism
    )
    saveDb(x = TxDb, file = paste(gtfDir, sqlite, sep = "/"))
}

#' Assign gene and transcript biotypes
#' @description Internal helper for developAnnotations(). Adds gene_biotype
#'     and transcript_biotype columns to GTFs based on ORF presence.
#' @param GTFs GRanges object with experimental annotations.
#' @param ORFs data.frame with ORF coordinates.
#' @return GRanges object with biotype columns.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.assignBiotypes <- function(GTFs, ORFs) {
    GTFs$gene_biotype <- "noncoding"
    GTFs$transcript_biotype <- "noncoding"
    GTFs[GTFs$transcript_id %in% ORFs$transcript_id, ]$transcript_biotype <-
        "protein_coding"
    GTFs[GTFs$gene_id %in%
        GTFs[GTFs$transcript_biotype ==
            "protein_coding", ]$gene_id, ]$gene_biotype <- "protein_coding"
    GTFs
}

#' Build gene-level ranges
#' @description Internal helper for developAnnotations(). Collapses
#'     transcript-level annotations into gene-level ranges.
#' @param GTFs GRanges object with experimental annotations.
#' @return list of GRanges objects, one per gene.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.buildGeneRanges <- function(GTFs) {
    genes <- split(x = GTFs, f = GTFs$gene_id)
    lapply(
        X = genes,
        FUN = function(y) {
            start(x = y) <- min(x = start(x = y))
            end(x = y) <- max(x = end(x = y))
            y <- y[1]
            y$type <- "gene"
            mcols(x = y) <-
                mcols(x = y)[!names(x = mcols(x = y)) %in%
                    c(
                        "transcript_id", "exon_number",
                        "transcript_biotype"
                    )]
            return(y)
        }
    )
}

#' Process coding transcripts
#' @description Internal helper for developAnnotations(). Adds CDS, start
#'     and stop codons, UTRs to coding transcripts.
#' @param coding GRanges object with coding transcripts.
#' @param genes list of gene-level GRanges.
#' @param ORFs data.frame with ORF coordinates.
#' @return list of processed coding transcripts.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.processCodingTranscripts <- function(coding, genes, ORFs) {
    coding <- split(x = coding, f = coding$transcript_id)
    for (i in seq_along(coding)) {
        tr <- coding[[i]]
        orf_start <- ORFs[
            ORFs$transcript_id == tr$transcript_id[1],
        ]$orf_start
        orf_end <- ORFs[
            ORFs$transcript_id == tr$transcript_id[1],
        ]$orf_end
        if (as.character(x = strand(x = tr))[1] == "+") {
            tr <- .processPlusStrandTranscript(tr, orf_start, orf_end, genes)
        }
        if (as.character(x = strand(x = tr))[1] == "-") {
            tr <- .processMinusStrandTranscript(tr, orf_start, orf_end, genes)
        }
        coding[[i]] <- tr
    }
    unlist(x = coding)
}

#' Process plus-strand coding transcript
#' @description Internal helper for developAnnotations(). Adds CDS, start
#'     and stop codons, UTRs to a plus-strand coding transcript.
#' @param tr GRanges object with a single transcript.
#' @param orf_start integer value, start of ORF in transcript coordinates.
#' @param orf_end integer value, end of ORF in transcript coordinates.
#' @param genes list of gene-level GRanges.
#' @return GRanges object with transcript features.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.processPlusStrandTranscript <- function(tr, orf_start, orf_end, genes) {
    codons <- .buildCodonsPlus(tr, orf_start, orf_end)
    tr <- sort(x = tr)
    CDS_exons <- .buildCdsExons(tr, codons$start_codon, codons$stop_codon)
    start_codon <- .buildCodonRange(
        tr, codons$start_codon, "start_codon"
    )
    stop_codon <- .buildCodonRange(
        tr, codons$stop_codon, "stop_codon"
    )
    utrs <- .buildUtrs(tr, CDS_exons, "+")
    transcript <- .buildTranscriptRange(tr)
    gene <- genes[names(x = genes) == unique(x = tr$gene_id)][[1]]
    c(
        gene,
        transcript,
        sort(c(
            tr, utrs$five_prime_utr, start_codon, CDS_exons,
            stop_codon, utrs$three_prime_utr
        ))
    )
}

#' Process minus-strand coding transcript
#' @description Internal helper for developAnnotations(). Adds CDS, start
#'     and stop codons, UTRs to a minus-strand coding transcript.
#' @param tr GRanges object with a single transcript.
#' @param orf_start integer value, start of ORF in transcript coordinates.
#' @param orf_end integer value, end of ORF in transcript coordinates.
#' @param genes list of gene-level GRanges.
#' @return GRanges object with transcript features.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.processMinusStrandTranscript <- function(tr, orf_start, orf_end, genes) {
    codons <- .buildCodonsMinus(tr, orf_start, orf_end)
    tr <- sort(x = tr)
    tr$exon_number <- rev(x = tr$exon_number)
    CDS_exons <- .buildCdsExons(tr, codons$start_codon, codons$stop_codon)
    start_codon <- .buildCodonRange(
        tr, codons$start_codon, "start_codon"
    )
    stop_codon <- .buildCodonRange(
        tr, codons$stop_codon, "stop_codon"
    )
    utrs <- .buildUtrs(tr, CDS_exons, "-")
    transcript <- .buildTranscriptRange(tr)
    gene <- genes[names(x = genes) == unique(x = tr$gene_id)][[1]]
    c(
        gene,
        transcript,
        sort(c(
            tr, utrs$five_prime_utr, start_codon, CDS_exons,
            stop_codon, utrs$three_prime_utr
        ), decreasing = TRUE)
    )
}

#' Process noncoding transcripts
#' @description Internal helper for developAnnotations(). Adds transcript
#'     and gene ranges to noncoding transcripts.
#' @param noncoding GRanges object with noncoding transcripts.
#' @param genes list of gene-level GRanges.
#' @return list of processed noncoding transcripts.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.processNoncodingTranscripts <- function(noncoding, genes) {
    noncoding <- split(x = noncoding, f = noncoding$transcript_id)
    for (i in seq_along(noncoding)) {
        tr <- noncoding[[i]]
        if (as.character(x = strand(x = tr))[1] == "-") {
            tr <- sort(x = tr)
            tr$exon_number <- rev(x = tr$exon_number)
        } else {
            tr <- sort(x = tr)
        }
        transcript <- .buildTranscriptRange(tr)
        gene <- genes[names(x = genes) == unique(x = tr$gene_id)][[1]]
        noncoding[[i]] <- c(gene, transcript, tr)
    }
    unlist(x = noncoding)
}

#' Combine transcript-level ranges into final GTFs
#' @description Internal helper for developAnnotations(). Recombines
#'     processed coding and noncoding transcripts into a single GRanges
#'     object, with gene ranges first.
#' @param coding list of processed coding transcripts.
#' @param noncoding list of processed noncoding transcripts.
#' @return GRanges object with final annotations.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.combineTranscripts <- function(coding, noncoding) {
    GTFs <- c(coding, noncoding)
    GTFs <- split(x = GTFs, f = GTFs$gene_id)
    for (i in seq_along(GTFs)) {
        GTFs[[i]] <- c(
            GTFs[[i]][1],
            GTFs[[i]][GTFs[[i]]$type != "gene", ]
        )
        names(x = GTFs[[i]]) <- NULL
    }
    unlist(x = GTFs)
}

#' Build start/stop codons for plus-strand transcript
#' @description Internal helper for developAnnotations(). Computes the
#'     genomic coordinates of the start and stop codons for a plus-strand
#'     transcript.
#' @param tr GRanges object with a single transcript.
#' @param orf_start integer value, start of ORF in transcript coordinates.
#' @param orf_end integer value, end of ORF in transcript coordinates.
#' @return list with start_codon and stop_codon integer vectors.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.buildCodonsPlus <- function(tr, orf_start, orf_end) {
    if (orf_start > width(x = tr)[1]) {
        m1 <- min(x = which(x = cumsum(x = width(x = tr)) - orf_start > 0))
        m2 <- max(x = cumsum(x = width(x = tr)[seq_len(m1 - 1)]))
        start_codon <- start(x = tr)[m1] + orf_start - m2 - 1
        start_codon <- c(start_codon, start_codon + 2)
        m1 <- min(x = which(x = cumsum(x = width(x = tr)) - orf_end > 0))
        m2 <- max(x = cumsum(x = width(x = tr)[seq_len(m1 - 1)]))
        stop_codon <- start(x = tr)[m1] + orf_end - m2 - 3
        stop_codon <- c(stop_codon, stop_codon + 2)
    } else {
        start_codon <- start(x = tr)[1] + orf_start - 1
        start_codon <- c(start_codon, start_codon + 2)
        if (length(x = tr) == 1) {
            stop_codon <- start(x = tr)[1] + orf_end - 1
            stop_codon <- c(stop_codon - 2, stop_codon)
        } else {
            m1 <- min(x = which(x = cumsum(x = width(x = tr)) - orf_end > 0))
            if (m1 == 1) {
                stop_codon <- start(x = tr)[1] + orf_end - 1
                stop_codon <- c(stop_codon - 2, stop_codon)
            } else {
                m2 <- max(x = cumsum(x = width(x = tr)[seq_len(m1 - 1)]))
                stop_codon <- start(x = tr)[m1] + orf_end - m2 - 3
                stop_codon <- c(stop_codon, stop_codon + 2)
            }
        }
    }
    list(start_codon = start_codon, stop_codon = stop_codon)
}

#' Build start/stop codons for minus-strand transcript
#' @description Internal helper for developAnnotations(). Computes the
#'     genomic coordinates of the start and stop codons for a minus-strand
#'     transcript.
#' @param tr GRanges object with a single transcript.
#' @param orf_start integer value, start of ORF in transcript coordinates.
#' @param orf_end integer value, end of ORF in transcript coordinates.
#' @return list with start_codon and stop_codon integer vectors.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.buildCodonsMinus <- function(tr, orf_start, orf_end) {
    tr <- sort(x = tr, decreasing = TRUE)
    if (orf_start > width(x = tr)[1]) {
        m1 <- min(x = which(x = cumsum(x = width(x = tr)) - orf_start > 0))
        m2 <- max(x = cumsum(x = width(x = tr)[seq_len(m1 - 1)]))
        start_codon <- end(x = tr)[m1] - (orf_start - m2) + 1
        start_codon <- c(start_codon - 2, start_codon)
        m1 <- min(x = which(x = cumsum(x = width(x = tr)) - orf_end > 0))
        m2 <- max(x = cumsum(x = width(x = tr)[seq_len(m1 - 1)]))
        stop_codon <- end(x = tr)[m1] - (orf_end - m2) + 1
        stop_codon <- c(stop_codon, stop_codon + 2)
    } else {
        start_codon <- end(x = tr)[1] - orf_start + 1
        start_codon <- c(start_codon - 2, start_codon)
        if (length(x = tr) == 1) {
            stop_codon <- end(x = tr)[1] - orf_end + 1
            stop_codon <- c(stop_codon, stop_codon + 2)
        } else {
            if (tail(x = cumsum(x = width(x = tr)), n = 1) == orf_end) {
                stop_codon <- start(x = tail(x = tr, n = 1))
                stop_codon <- c(stop_codon, stop_codon + 2)
            } else {
                m1 <- min(x = which(
                    x = cumsum(x = width(x = tr)) - orf_end > 0
                ))
                if (m1 == 1) {
                    stop_codon <- end(x = tr)[1] - orf_end + 1
                    stop_codon <- c(stop_codon, stop_codon + 2)
                } else {
                    m2 <- max(x = cumsum(x = width(x = tr)[seq_len(m1 - 1)]))
                    stop_codon <- end(x = tr)[m1] - (orf_end - m2) + 1
                    stop_codon <- c(stop_codon, stop_codon + 2)
                }
            }
        }
    }
    list(start_codon = start_codon, stop_codon = stop_codon)
}

#' Build CDS exons
#' @description Internal helper for developAnnotations(). Constructs CDS
#'     exon ranges from transcript exons and codon coordinates.
#' @param tr GRanges object with a single transcript.
#' @param start_codon integer vector of length 2 with start codon
#'     coordinates.
#' @param stop_codon integer vector of length 2 with stop codon
#'     coordinates.
#' @return GRanges object with CDS exon ranges.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.buildCdsExons <- function(tr, start_codon, stop_codon) {
    strand_val <- as.character(x = strand(x = tr))[1]
    if (strand_val == "+") {
        cds_start <- min(x = start_codon)
        cds_end <- min(x = stop_codon) - 1
    } else {
        cds_start <- max(x = stop_codon) + 1
        cds_end <- max(x = start_codon)
    }
    CDS <- data.frame(
        seqnames = as.character(x = seqnames(x = tr))[1],
        start = cds_start,
        end = cds_end,
        strand = strand_val
    )
    CDS <- makeGRangesFromDataFrame(df = CDS)
    hits <- findOverlaps(query = CDS, subject = tr, type = "any")
    CDS_exons <- tr[subjectHits(x = hits)]
    CDS_exons$type <- "CDS"
    start(x = CDS_exons)[1] <- start(x = CDS)
    end(x = CDS_exons)[length(x = CDS_exons)] <- end(x = CDS)
    CDS_exons
}

#' Build codon range
#' @description Internal helper for developAnnotations(). Constructs a
#'     GRanges object representing a single codon (start or stop).
#' @param tr GRanges object with a single transcript.
#' @param codon integer vector of length 2 with codon coordinates.
#' @param codon_type character string specifying codon type.
#' @return GRanges object with a codon range.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.buildCodonRange <- function(tr, codon, codon_type) {
    codon_range <- data.frame(
        seqnames = as.character(x = seqnames(x = tr))[1],
        start = codon[1],
        end = codon[2],
        strand = as.character(x = strand(x = tr))[1]
    )
    codon_range <- makeGRangesFromDataFrame(df = codon_range)
    hits <- findOverlaps(
        query = codon_range,
        subject = tr,
        type = "any"
    )
    codon_exons <- tr[subjectHits(x = hits)]
    codon_exons$type <- codon_type
    start(x = codon_exons)[1] <- start(x = codon_range)
    end(x = codon_exons)[1] <- end(x = codon_range)
    codon_exons
}

#' Build UTR ranges
#' @description Internal helper for developAnnotations(). Constructs the
#'     5' and 3' UTR ranges from transcript exons and CDS.
#' @param tr GRanges object with a single transcript.
#' @param CDS_exons GRanges object with CDS exons.
#' @param strand character string, "+" or "-".
#' @return list with five_prime_utr and three_prime_utr GRanges objects.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.buildUtrs <- function(tr, CDS_exons, strand) {
    UTRs <- disjoin(c(tr, CDS_exons))
    hits <- findOverlaps(query = UTRs, subject = CDS_exons, type = "any")
    UTRs <- UTRs[-queryHits(x = hits)]
    if (length(x = UTRs) == 0) {
        return(list(five_prime_utr = UTRs, three_prime_utr = UTRs))
    }
    cds_start <- min(x = start(x = CDS_exons))
    cds_end <- max(x = end(x = CDS_exons))
    if (strand == "+") {
        is_five <- end(x = UTRs) < cds_start
        is_three <- start(x = UTRs) > cds_end
    } else {
        is_five <- start(x = UTRs) > cds_end
        is_three <- end(x = UTRs) < cds_start
    }
    five_prime_utr <- UTRs[is_five]
    three_prime_utr <- UTRs[is_three]
    if (length(x = five_prime_utr) > 0) {
        hits <- findOverlaps(
            query = five_prime_utr, subject = tr, type = "any"
        )
        mcols(x = five_prime_utr) <-
            mcols(x = tr[subjectHits(x = hits)])
        five_prime_utr$type <- "five_prime_utr"
    }
    if (length(x = three_prime_utr) > 0) {
        hits <- findOverlaps(
            query = three_prime_utr, subject = tr, type = "any"
        )
        mcols(x = three_prime_utr) <-
            mcols(x = tr[subjectHits(x = hits)])
        three_prime_utr$type <- "three_prime_utr"
        if (strand == "+") {
            start(x = three_prime_utr) <- start(x = three_prime_utr) + 3
        } else {
            end(x = three_prime_utr) <- end(x = three_prime_utr) - 3
        }
    }
    list(
        five_prime_utr = five_prime_utr,
        three_prime_utr = three_prime_utr
    )
}

#' Check that Perl is available
#' @description Internal helper for scoreSSs(). Stops if Perl is not in
#'     PATH (MaxEntScan requires it).
#' @return invisible(NULL).
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @keywords internal
#' @noRd
.checkPerlAvailable <- function() {
    if (Sys.which("perl") == "") {
        stop(
            "Perl not found in PATH. scoreSSs() requires Perl ",
            "for MaxEntScan. Install Perl (e.g., Strawberry Perl on ",
            "Windows) and ensure it is available in PATH."
        )
    }
    invisible(NULL)
}

#' Assemble the MaxEntScan score table
#' @description Internal helper for scoreSSs(). Builds the final
#'     data.frame with IDs, sequences and scores.
#' @param seqSSs DNAStringSet object with splice site sequences.
#' @param scores character vector with raw MaxEntScan scores.
#' @return data.frame with ss_id, ss_seq and max.ent_score.
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @keywords internal
#' @noRd
.assembleScoreTable <- function(seqSSs, scores) {
    scores <- cbind(
        names(x = seqSSs),
        as.vector(x = as.character(x = seqSSs)),
        scores
    )
    colnames(x = scores) <- c("ss_id", "ss_seq", "max.ent_score")
    scores <- data.frame(scores)
    scores[, 3] <- as.numeric(scores[, 3])
    scores[!duplicated(x = scores), ]
}

#' Build transcript range
#' @description Internal helper for developAnnotations(). Constructs a
#'     single transcript range with type = "transcript".
#' @param tr GRanges object with a single transcript.
#' @return GRanges object with a transcript range.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.buildTranscriptRange <- function(tr) {
    transcript <- tr[1]
    mcols(x = transcript) <-
        mcols(x = transcript)[
            names(x = mcols(x = transcript)) != "exon_number"
        ]
    transcript$type <- "transcript"
    start(x = transcript) <- min(x = start(x = tr))
    end(x = transcript) <- max(x = end(x = tr))
    transcript
}


#' Extract of splice site sequences from a reference genome
#' @description Extraction of splice site sequences from the reference
#'     genome.
#' @param x genomic coordinates of 5' or 3' splice sites. It can be the
#'     name of object of class GRanges (as output of
#'     extractSSsCoordsEEJs() or extractSSsCoordsGTF() function) or
#'     tab-delimited TXT file with the following fields:
#'     i) seqnames - name of chromosome or scaffold with prefix "chr";
#'     ii) start   - start coordinate of splice site;
#'     iii) end    - end coordinate of splice site;
#'     iv) strand  - strand information about splice site;
#'     v-...) ...  - (optionally) any metadata.
#'     If TXT file, it should be in work directory.
#' @param genome character string giving the name of the BSgenome data
#'     package to be used for reference sequence extraction. Default
#'     value is package BSgenome.Hsapiens.UCSC.hg38.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return object of class DNAStringSet containing splice site sequences.
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @examples
#' if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
#'     coords <- GenomicRanges::GRanges("chr1:100000-100008:+")
#'     extractSeqSSs(x = coords)
#' }
#' @export
#' @importFrom BSgenome getSeq
#' @importFrom GenomicRanges makeGRangesFromDataFrame seqnames start
#'     end strand
#' @importFrom methods is
extractSeqSSs <- function(
        x,
        genome = "BSgenome.Hsapiens.UCSC.hg38",
        workDir = NULL) {
    ### Setting the working directory.
    if (is.null(x = workDir)) {
        workDir <- getwd()
    }
    ### Loading of the genomic coordinates of splice sites.
    coordsSSs <- .loadSpliceSiteCoords(x, workDir)
    ### Loading of the reference genome object.
    if (!requireNamespace(genome, quietly = TRUE)) {
        stop(
            "Package ", genome, " is not installed. ",
            "Install it with BiocManager::install('", genome, "')."
        )
    }
    genome_obj <- get(x = genome, envir = asNamespace(genome))
    ### Extraction of splice site sequences from the reference genome.
    seqSSs <- getSeq(x = genome_obj, names = coordsSSs)
    names(x = seqSSs) <- .spliceSiteSeqNames(coordsSSs)
    ### Returning the final object.
    return(seqSSs)
}

#' Load splice site coordinates
#' @description Internal helper for extractSeqSSs(). Returns a GRanges
#'     object either directly (when x is already an object) or by reading
#'     a tab-delimited TXT file.
#' @param x GRanges object or character string with the file name.
#' @param workDir character string with the working directory.
#' @return GRanges object with splice site coordinates.
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @keywords internal
#' @noRd
.loadSpliceSiteCoords <- function(x, workDir) {
    if (is.object(x = x)) {
        return(x)
    }
    frt <- tools::file_ext(x = x)
    if (!frt %in% "txt") {
        stop("Invalid file format")
    }
    coordsSSs <- read.table(
        file = paste(workDir, x, sep = "/"),
        sep = "\t",
        header = TRUE,
        quote = "\"",
        as.is = TRUE
    )
    makeGRangesFromDataFrame(
        df = coordsSSs,
        keep.extra.columns = TRUE
    )
}

#' Build splice site sequence names
#' @description Internal helper for extractSeqSSs(). Builds the sequence
#'     names in the form "chr:start-end_strand".
#' @param coordsSSs GRanges object with splice site coordinates.
#' @return character vector with sequence names.
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @keywords internal
#' @noRd
.spliceSiteSeqNames <- function(coordsSSs) {
    paste(
        paste(seqnames(x = coordsSSs),
            paste(start(x = coordsSSs),
                end(x = coordsSSs),
                sep = "-"
            ),
            sep = ":"
        ),
        strand(x = coordsSSs),
        sep = "_str"
    )
}

#' Extract genomic coordinates of splice sites from exon-exon junctions
#' @description Extraction of the genomic coordinates of splice sites
#'     using experimentally detected exon-exon junctions.
#' @param x character string specifying the name of the tab-delimited TXT
#'     file with experimentally detected exon-exon junctions. This file
#'     must contain the following fields:
#'     i) eej_id     - exon-exon junction ID;
#'     ii) gene_id   - gene ID;
#'     iii) seqnames - name of chromosome or scaffold with prefix "chr";
#'     iv) start     - start coordinate of exon-exon junction (genomic
#'                     position of the last nucleotide in the upstream
#'                     exon);
#'     v) end        - end coordinate of exon-exon junction (genomic
#'                     position of the first nucleotide in the downstream
#'                     exon);
#'     vi) strand    - strand information about exon-exon junction
#'                     location;
#'     vii-...) ...  - (optionally) one or more samples with summarized
#'                     raw RNA-Seq reads spanning exon-exon junction.
#' @param eejDir character string specifying the name of the directory
#'     containing the tab-delimited TXT file(-s) with experimentally
#'     detected exon-exon junctions. NULL by default that mean the
#'     current working directory.
#' @param thr an integer argument. It is a threshold for minimal coverage
#'     (sequencing depth) of exon-exon junction (in raw reads). Default
#'     value is 5 but can be specified by the user.
#' @param workDir character string specifying the path to and name of
#'     working directory. NULL by default that means the current working
#'     directory.
#' @return object of class list containing two GRanges objects with
#'     genomic coordinates of 5' and 3' splice sites.
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @examples
#' ex <- system.file("extdata", package = "huGSVcalleR")
#' extractSSsCoordsEEJs(
#'     x = "test_eejs.txt",
#'     eejDir = NULL,
#'     thr = 5,
#'     workDir = ex
#' )
#' @export
#' @importFrom GenomicRanges makeGRangesFromDataFrame start end strand
#'     start<- end<-
#' @importFrom methods is
extractSSsCoordsEEJs <- function(
        x,
        eejDir = NULL,
        thr = 5,
        workDir = NULL) {
    ### Setting the working directory.
    if (is.null(x = workDir)) {
        workDir <- getwd()
    }
    ### Full path to the file.
    if (is.null(x = eejDir)) {
        path <- workDir
    } else {
        path <- paste(workDir, eejDir, sep = "/")
    }
    ### Retrieving and validation of the file extension.
    frt <- tools::file_ext(x = x)
    if (!frt %in% "txt") {
        stop("Invalid file format")
    }
    ### Loading of exon-exon junctions as an object of class GRanges.
    EEJs <- read.table(
        file = paste(path, x, sep = "/"),
        sep = "\t",
        header = TRUE,
        quote = "\"",
        as.is = TRUE
    )
    if (ncol(x = EEJs) > 7 & !is.null(x = thr)) {
        EEJs <- EEJs[
            rowSums(x = EEJs[, -seq_len(7)] >= thr) >= 1,
        ]
    }
    rownames(x = EEJs) <- NULL
    EEJs <- makeGRangesFromDataFrame(
        df = EEJs,
        keep.extra.columns = TRUE
    )
    ### Converting of genomic coordinates of exon-exon junctions
    #   in genomic coordinates of splice sites.
    .eejToSpliceSites(EEJs)
}

#' Convert exon-exon junction coordinates to splice site coordinates
#' @description Internal helper for extractSSsCoordsEEJs() and
#'     extractSSsCoordsGTF(). Converts the genomic coordinates of
#'     exon-exon junctions into the genomic coordinates of 5' and 3'
#'     splice sites.
#' @param EEJs GRanges object with exon-exon junction coordinates.
#' @return list with fiveSSs and threeSSs GRanges objects.
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @keywords internal
#' @noRd
.eejToSpliceSites <- function(EEJs) {
    fiveSSs <- EEJs
    start(x = fiveSSs[strand(x = fiveSSs) == "+", ]) <-
        start(x = fiveSSs[strand(x = fiveSSs) == "+", ]) - 2
    end(x = fiveSSs[strand(x = fiveSSs) == "+", ]) <-
        start(x = fiveSSs[strand(x = fiveSSs) == "+", ]) + 8
    start(x = fiveSSs[strand(x = fiveSSs) == "-", ]) <-
        end(x = fiveSSs[strand(x = fiveSSs) == "-", ]) - 6
    end(x = fiveSSs[strand(x = fiveSSs) == "-", ]) <-
        start(x = fiveSSs[strand(x = fiveSSs) == "-", ]) + 8
    threeSSs <- EEJs
    start(x = threeSSs[strand(x = threeSSs) == "+", ]) <-
        end(x = threeSSs[strand(x = threeSSs) == "+", ]) - 20
    end(x = threeSSs[strand(x = threeSSs) == "+", ]) <-
        start(x = threeSSs[strand(x = threeSSs) == "+", ]) + 22
    start(x = threeSSs[strand(x = threeSSs) == "-", ]) <-
        start(x = threeSSs[strand(x = threeSSs) == "-", ]) - 2
    end(x = threeSSs[strand(x = threeSSs) == "-", ]) <-
        start(x = threeSSs[strand(x = threeSSs) == "-", ]) + 22
    list(fiveSSs = fiveSSs, threeSSs = threeSSs)
}


#' Extract genomic coordinates of splice sites from annotations
#' @description Extraction of the genomic coordinates of splice sites
#'     using genomic/transcriptomic annotations.
#' @param x character string specifying the name of file with input data.
#'     Accepted formats: GFF2.5/GTF or GFF3 (the gz archive of GFF2.5/GTF
#'     or GFF3 file is also accepted). The file should include models of
#'     genes in general/generic feature format with mandatory field
#'     "transcript_id".
#' @param gtfDir character string specifying the name of the directory
#'     containing the file with input data. NULL by default that mean the
#'     current working directory.
#' @param workDir character string specifying the path to and name of
#'     working directory. NULL by default that means the current working
#'     directory.
#' @return object of class list containing two GRanges objects with
#'     genomic coordinates of 5' and 3' splice sites.
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @examples
#' ex <- system.file("extdata", package = "huGSVcalleR")
#' extractSSsCoordsGTF(
#'     x = "test_annotation.gtf",
#'     gtfDir = NULL,
#'     workDir = ex
#' )
#' @export
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges makeGRangesListFromDataFrame
#'     makeGRangesFromDataFrame
#' @importFrom GenomicRanges start end strand start<- end<-
#'     elementMetadata
#' @importFrom IRanges setdiff
#' @importFrom S4Vectors colnames
#' @importFrom methods as
extractSSsCoordsGTF <- function(x, gtfDir = NULL, workDir = NULL) {
    ### Setting the working directory.
    if (is.null(x = workDir)) {
        workDir <- getwd()
    }
    ### Full path to the file.
    if (is.null(x = gtfDir)) {
        path <- workDir
    } else {
        path <- paste(workDir, gtfDir, sep = "/")
    }
    ### Retrieving and validation of the file extension.
    frt <- tools::file_ext(x = x)
    if (!frt %in% c("gtf", "gff3", "gz")) {
        stop("Invalid file format")
    }
    ### Loading of exon-exon junctions as an object of class GRanges.
    EEJs <- import(con = paste(path, x, sep = "/"))
    EEJs <- .normalizeGtfExons(EEJs)
    ### Converting of genomic coordinates of exon-exon junctions
    #   in genomic coordinates of splice sites.
    SSs <- .eejToSpliceSites(EEJs)
    list(
        fiveSSs = sort(x = SSs$fiveSSs),
        threeSSs = sort(x = SSs$threeSSs)
    )
}

#' Normalize exon annotations from GTF
#' @description Internal helper for extractSSsCoordsGTF(). Extracts exon
#'     records from a GTF/GFF object, converts them to exon-exon junctions,
#'     and normalizes chromosome names to UCSC style.
#' @param EEJs GRanges object imported from GTF/GFF.
#' @return GRanges object with normalized exon-exon junctions.
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @keywords internal
#' @noRd
.normalizeGtfExons <- function(EEJs) {
    EEJs <- EEJs[gsub(
        pattern = ".+exon",
        replacement = "exon",
        x = EEJs$type
    ) == "exon", ]
    if ("transcript_id" %in%
        colnames(x = elementMetadata(x = EEJs)) == "TRUE") {
        splitField <- "transcript_id"
    } else {
        splitField <- "Name"
    }
    EEJs <- makeGRangesListFromDataFrame(
        df = as.data.frame(x = EEJs),
        split.field = splitField,
        keep.extra.columns = FALSE
    )
    EEJs <- unlist(x = setdiff(x = range(x = EEJs), EEJs))
    start(x = EEJs) <- start(x = EEJs) - 1
    end(x = EEJs) <- end(x = EEJs) + 1
    EEJs$transcript_id <- names(x = EEJs)
    names(x = EEJs) <- NULL
    EEJs <- as.data.frame(EEJs, stringsAsFactors = FALSE)
    EEJs$seqnames <- as.character(x = EEJs$seqnames)
    if (length(x = unique(x = nchar(EEJs$seqnames) == 1)) == 2) {
        EEJs$seqnames <- paste("chr", EEJs$seqnames, sep = "")
        if ("chrMT" %in% EEJs$seqnames == "TRUE") {
            EEJs[EEJs$seqnames == "chrMT", ]$seqnames <- "chrM"
        }
        EEJs <- EEJs[
            EEJs$seqnames %in%
                c(paste("chr", c(seq_len(22), "X", "Y", "M"), sep = "")),
        ]
    } else {
        if ("chrMT" %in% EEJs$seqnames == "TRUE") {
            EEJs[EEJs$seqnames == "chrMT", ]$seqnames <- "chrM"
        }
        EEJs <- EEJs[
            EEJs$seqnames %in%
                c(paste("chr", c(seq_len(22), "X", "Y", "M"), sep = "")),
        ]
    }
    EEJs$strand <- as.character(x = EEJs$strand)
    EEJs <- EEJs[!duplicated(x = EEJs), ]
    rownames(x = EEJs) <- NULL
    makeGRangesFromDataFrame(df = EEJs, keep.extra.columns = TRUE)
}


#' Predict amino acid coding changes for single nucleotide variants
#' @description Prediction of amino acid coding changes for single
#'     nucleotide variants.
#' @param vcfDir character string giving the name (or path to and name) of
#'     directory with VCF file(-s). NULL by default, which means the
#'     current working directory will be used instead of the specified
#'     VCF directory.
#' @param vcfFile character string giving the name of input VCF file.
#' @param OrgDb character string giving the name of the OrgDb data package
#'     with genome wide annotations. Default value is package
#'     org.Hs.eg.db.
#' @param TxDb character string giving the name of the TxDb data package
#'     with transcriptome wide annotations. Default value is package
#'     TxDb.Hsapiens.UCSC.hg38.knownGene.
#' @param genome character string giving the name of the BSgenome data
#'     package to be used for reference sequence extraction. Default
#'     value is package BSgenome.Hsapiens.UCSC.hg38.
#' @param proteins character string giving the name of TXT file in
#'     tab-delimited format with protein coding annotations of Ensembl
#'     transcripts. This file must contains the following IDs: Ensembl
#'     transcript IDs, Ensembl gene IDs, Ensembl protein _IDs, PDB IDs,
#'     UniProtKB/Swiss-Prot IDs and UniProtKB/TrEMBL IDs.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return object of class data frame with predicted and annotated amino
#'     acid coding changes for single nucleotide variants.
#' @author Vasily V. Grinev.
#' @examples
#' ready <- all(vapply(
#'     c(
#'         "BSgenome.Hsapiens.UCSC.hg38",
#'         "org.Hs.eg.db",
#'         "TxDb.Hsapiens.UCSC.hg38.knownGene"
#'     ),
#'     requireNamespace,
#'     logical(1),
#'     quietly = TRUE
#' ))
#' if (ready) {
#'     txdb <- getExportedValue(
#'         "TxDb.Hsapiens.UCSC.hg38.knownGene",
#'         "TxDb.Hsapiens.UCSC.hg38.knownGene"
#'     )
#'     genome <- getExportedValue(
#'         "BSgenome.Hsapiens.UCSC.hg38",
#'         "BSgenome.Hsapiens.UCSC.hg38"
#'     )
#'     cds <- GenomicFeatures::cdsBy(txdb, by = "tx")
#'     ok <- vapply(cds, function(gr) {
#'         as.character(BiocGenerics::strand(gr)[1]) == "+" &&
#'             sum(IRanges::width(gr)) >= 3
#'     }, logical(1))
#'     gr <- cds[[which(ok)[1]]]
#'     pos <- BiocGenerics::start(gr)[1]
#'     chr <- as.character(GenomicRanges::seqnames(gr)[1])
#'     ref <- as.character(BSgenome::getSeq(
#'         genome,
#'         GenomicRanges::GRanges(chr, IRanges::IRanges(pos, pos))
#'     ))
#'     alt <- paste(setdiff(c("A", "C", "G", "T"), ref), collapse = ",")
#'     workDir <- tempdir()
#'     writeLines(
#'         c(
#'             "##fileformat=VCFv4.2",
#'             sprintf("##contig=<ID=%s,length=%s>", chr, pos + 10),
#'             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
#'             paste(chr, pos, ".", ref, alt, "10", "PASS", ".", sep = "\t")
#'         ),
#'         file.path(workDir, "coding.vcf")
#'     )
#'     write.table(
#'         data.frame(
#'             transcript_id = "TX",
#'             gene_id = "GENE",
#'             protein_id = "PROT",
#'             pdb_id = NA_character_,
#'             swissprot = NA_character_,
#'             trembl = NA_character_
#'         ),
#'         file = file.path(workDir, "proteins.txt"),
#'         sep = "\t",
#'         quote = FALSE,
#'         row.names = FALSE
#'     )
#'     introduceSNVtoAAs(
#'         vcfDir = NULL,
#'         vcfFile = "coding.vcf",
#'         proteins = "proteins.txt",
#'         workDir = workDir
#'     )
#' }
#' @export
#' @importFrom VariantAnnotation readVcf alt alt<- predictCoding
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom Biostrings DNAStringSetList
#' @importFrom GenomicRanges start seqnames
#' @importFrom BiocGenerics "%in%" match
#' @importFrom AnnotationDbi mapIds
introduceSNVtoAAs <- function(vcfDir = NULL,
    vcfFile,
    OrgDb = "org.Hs.eg.db",
    TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    genome = "BSgenome.Hsapiens.UCSC.hg38",
    proteins,
    workDir = NULL) {
    ### Checking of required packages.
    models <- .loadPreTrainedModels(OrgDb, TxDb, genome)
    ### Setting the working directory and VCF directory.
    if (is.null(x = workDir)) {
        workDir <- getwd()
    }
    if (is.null(x = vcfDir)) {
        vcfDir <- workDir
    } else {
        vcfDir <- paste(workDir, vcfDir, sep = "/")
    }
    ### Loading of the single nucleotide variations, filtered to canonical
    #   chromosomes and with multi-allelic ALT expanded.
    vcf <- .loadCodingVcf(vcfDir, vcfFile)
    ### Loading of the protein annotations as an object of class data frame.
    pr <- .loadProteinAnnotations(workDir, proteins)
    ### Prediction of amino acid coding changes.
    AAs <- predictCoding(
        query = vcf,
        subject = models$txDb,
        seqSource = models$genome
    )
    AAs <- AAs[AAs$CONSEQUENCE == "nonsynonymous", ]
    ### Annotation.
    .buildAaAnnotationFrame(AAs, models$txDb, models$orgDb, pr)
}

#' Load VCF for coding annotation
#' @description Internal helper for introduceSNVtoAAs(). Reads the VCF,
#'     restricts it to canonical chromosomes, and expands multi-allelic
#'     ALT alleles into one record per allele.
#' @param vcfDir character string with the VCF directory.
#' @param vcfFile character string with the VCF file name.
#' @return CollapsedVCF object.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.loadCodingVcf <- function(vcfDir, vcfFile) {
    vcf <- readVcf(file = paste(vcfDir, vcfFile, sep = "/"))
    canonical <- paste0("chr", c(seq_len(22), "X", "Y", "M"))
    vcf <- vcf[seqnames(rowRanges(vcf)) %in% canonical]
    ALT <- as.character(x = unlist(x = rowRanges(x = vcf)$ALT))
    vcf <- vcf[
        rep(
            x = seq_along(vcf),
            times = lengths(x = rowRanges(x = vcf)$ALT)
        ),
    ]
    alt(x = vcf) <- DNAStringSetList(strsplit(x = ALT, split = ""))
    vcf
}

#' Load protein annotations
#' @description Internal helper for introduceSNVtoAAs(). Reads the
#'     tab-delimited protein annotation table and uses its first column
#'     as row names.
#' @param workDir character string with the working directory.
#' @param proteins character string with the protein table file name.
#' @return data.frame with protein annotations.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.loadProteinAnnotations <- function(workDir, proteins) {
    pr <- read.table(
        file = paste(workDir, proteins, sep = "/"),
        sep = "\t",
        header = TRUE,
        quote = "\"",
        as.is = TRUE
    )
    rownames(x = pr) <- pr[, 1]
    pr
}

#' Load pre-trained annotation models
#' @description Internal helper for introduceSNVtoAAs(). Checks that the
#'     required OrgDb, TxDb and BSgenome packages are installed and
#'     returns them as a list.
#' @param OrgDb character string giving the name of the OrgDb package.
#' @param TxDb character string giving the name of the TxDb package.
#' @param genome character string giving the name of the BSgenome package.
#' @return list with txDb, orgDb and genome objects.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.loadPreTrainedModels <- function(OrgDb, TxDb, genome) {
    for (pkg in c(genome, OrgDb, TxDb)) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop(
                "Package '", pkg, "' is required but not installed. ",
                "Install it with BiocManager::install('", pkg, "')."
            )
        }
    }
    list(
        txDb = get(x = TxDb, envir = asNamespace(TxDb)),
        orgDb = get(x = OrgDb, envir = asNamespace(OrgDb)),
        genome = get(x = genome, envir = asNamespace(genome))
    )
}

#' Build amino acid annotation data.frame
#' @description Internal helper for introduceSNVtoAAs(). Adds gene and
#'     protein identifiers, codons, and amino acid substitutions to the
#'     predictCoding() output.
#' @param AAs data.frame from predictCoding().
#' @param txDb TxDb object with transcriptome annotations.
#' @param orgDb OrgDb object with genome-wide annotations.
#' @param pr data.frame with protein annotations.
#' @return data.frame with the final amino acid annotation.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.buildAaAnnotationFrame <- function(AAs, txDb, orgDb, pr) {
    AAs <- .addAaIdentifierColumns(AAs, txDb, orgDb)
    AAs <- .addAaCodonColumns(AAs)
    AAs <- .reorderAaColumns(AAs)
    AAs$seqnames <- as.character(x = AAs$seqnames)
    AAs$strand <- as.character(x = AAs$strand)
    idx <- AAs[, 13] %in% pr[, 1]
    AAs[idx, 17:20] <- pr[AAs[idx, 13], 3:6]
    AAs
}

#' Add identifier columns to the amino acid annotation
#' @description Internal helper for .buildAaAnnotationFrame(). Adds the
#'     SNV identifiers and the Ensembl/Entrez/HUGO gene and transcript
#'     identifiers via mapIds().
#' @param AAs data.frame from predictCoding().
#' @param txDb TxDb object with transcriptome annotations.
#' @param orgDb OrgDb object with genome-wide annotations.
#' @return data.frame with the added identifier columns.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.addAaIdentifierColumns <- function(AAs, txDb, orgDb) {
    AAs$snv_id <- names(x = AAs)
    AAs$snv_ref <- as.character(x = AAs$REF)
    AAs$snv_alt <- as.character(x = unlist(x = AAs$ALT))
    AAs$snv_qual <- AAs$QUAL
    AAs$Ensembl_gene_IDs <- mapIds(
        x = orgDb,
        keys = AAs$GENEID,
        column = "ENSEMBL",
        keytype = "ENTREZID",
        multiVals = "first"
    )
    AAs$Entrez_gene_IDs <- AAs$GENEID
    AAs$HUGO_gene_symbols <- mapIds(
        x = orgDb,
        keys = AAs$GENEID,
        column = "SYMBOL",
        keytype = "ENTREZID",
        multiVals = "first"
    )
    AAs$Ensembl_tx_names <- mapIds(
        x = txDb,
        keys = AAs$TXID,
        column = "TXNAME",
        keytype = "TXID",
        multiVals = "first"
    )
    AAs$Ensembl_tx_names <- substr(
        x = AAs$Ensembl_tx_names,
        start = 1,
        stop = 15
    )
    AAs
}

#' Add codon and amino acid columns
#' @description Internal helper for .buildAaAnnotationFrame(). Adds the
#'     SNV CDS location, reference and alternative codons, the empty
#'     protein/PDB/UniProt placeholders, and the amino acid substitutions.
#' @param AAs data.frame with the identifier columns already added.
#' @return data.frame with the added codon and amino acid columns.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.addAaCodonColumns <- function(AAs) {
    AAs$SNV_CDS_loc <- start(x = AAs$CDSLOC)
    AAs$ref_codon <- as.vector(x = as.character(x = AAs$REFCODON))
    AAs$alt_codon <- as.vector(x = as.character(x = AAs$VARCODON))
    AAs$protein_id <- ""
    AAs$pdb_id <- ""
    AAs$uniprot.swissprot_id <- ""
    AAs$uniprot.trembl_id <- ""
    AAs$SNV_prot_loc <- as.vector(x = unlist(x = AAs$PROTEINLOC))
    AAs$ref_aa <- as.vector(x = as.character(x = AAs$REFAA))
    AAs$alt_aa <- as.vector(x = as.character(x = AAs$VARAA))
    AAs
}

#' Reorder columns in the amino acid annotation
#' @description Internal helper for .buildAaAnnotationFrame(). Reorders the
#'     columns into the final layout expected by the downstream code.
#' @param AAs data.frame with all annotation columns.
#' @return data.frame with the reordered columns.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.reorderAaColumns <- function(AAs) {
    data.frame(AAs)[, c(
        seq(from = 25, to = 28),
        seq_len(5),
        seq(from = 29, to = 42)
    )]
}


#' Introduce single nucleotide variant(-s) to splice site(-s)
#' @description Introduction of single nucleotide variant(-s) into the
#'     splice site sequence(-s).
#' @param x genomic coordinates of 5' or 3' splice sites. It can be the
#'     name of object of class GRanges (as output of
#'     extractSSsCoordsEEJs() or extractSSsCoordsGTF() function) or
#'     tab-delimited TXT file with the following fields:
#'     i) seqnames - name of chromosome or scaffold with prefix "chr";
#'     ii) start   - start coordinate of splice site;
#'     iii) end    - end coordinate of splice site;
#'     iv) strand  - strand information about splice site;
#'     v-...) ...  - (optionally) any metadata.
#'     If TXT file, it should be in work directory.
#' @param vcfDir character string giving the name (or path to and name) of
#'     directory with VCF file(-s). NULL by default, which means the
#'     current working directory will be used instead of the specified
#'     VCF directory.
#' @param vcfFile character string giving the name of input VCF file.
#' @param genome character string giving the name of the BSgenome data
#'     package to be used for reference sequence extraction. Default
#'     value is package BSgenome.Hsapiens.UCSC.hg38.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return object of class list containing extended data frame with
#'     splice sites and DNAStringSet with corrected splice site sequences.
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @examples
#' if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
#'     genome <- getExportedValue(
#'         "BSgenome.Hsapiens.UCSC.hg38",
#'         "BSgenome.Hsapiens.UCSC.hg38"
#'     )
#'     pos <- 100000L
#'     ref <- as.character(BSgenome::getSeq(
#'         genome,
#'         GenomicRanges::GRanges("chr1", IRanges::IRanges(pos, pos))
#'     ))
#'     alt <- setdiff(c("A", "C", "G", "T"), ref)[1]
#'     workDir <- tempdir()
#'     writeLines(
#'         c(
#'             "##fileformat=VCFv4.2",
#'             "##contig=<ID=chr1,length=248956422>",
#'             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
#'             paste("chr1", pos, ".", ref, alt, "10", "PASS", ".", sep = "\t")
#'         ),
#'         file.path(workDir, "splice.vcf")
#'     )
#'     introduceSNVtoSSs(
#'         x = GenomicRanges::GRanges(
#'             "chr1:99995-100010:+",
#'             ss_id = "ss1",
#'             gene_id = "GENE1",
#'             transcript_id = "TX1"
#'         ),
#'         vcfDir = NULL,
#'         vcfFile = "splice.vcf",
#'         workDir = workDir
#'     )
#' }
#' @export
#' @importFrom VariantAnnotation readVcf
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom GenomicRanges start end strand seqnames
#' @importFrom Biostrings DNAStringSet reverseComplement
introduceSNVtoSSs <- function(
        x,
        vcfDir = NULL,
        vcfFile,
        genome = "BSgenome.Hsapiens.UCSC.hg38",
        workDir = NULL) {
    ### Setting the working directory and VCF directory.
    if (is.null(x = workDir)) {
        workDir <- getwd()
    }
    if (is.null(x = vcfDir)) {
        vcfDir <- workDir
    } else {
        vcfDir <- paste(workDir, vcfDir, sep = "/")
    }
    ### Loading of the single nucleotide variations.
    vcf <- .loadSnvAltFromVcf(vcfDir, vcfFile)
    ### Loading of the genomic coordinates of splice sites.
    spl.sites <- .loadSpliceSites(x, workDir)
    ### Find overlaps between splice sites and SNVs.
    hits <- findOverlaps(
        query = spl.sites,
        subject = vcf,
        type = "any",
        ignore.strand = TRUE
    )
    if (length(x = hits) == 0) {
        warning("No matches were found between splice sites and SNVs.")
        return(invisible(NULL))
    }
    spl.sites <- spl.sites[queryHits(x = hits)]
    vcf <- vcf[subjectHits(x = hits)]
    ### Correction of the splice site sequences.
    spl.sites <- .correctSpliceSiteSequences(
        spl.sites, vcf, genome, workDir
    )
    ### Returning the final object.
    .buildIntroducedSpliceSiteResult(spl.sites)
}

#' Load splice site coordinates
#' @description Internal helper for introduceSNVtoSSs(). Loads splice
#'     site coordinates from a GRanges object or a tab-delimited TXT file.
#' @param x GRanges object or character string with the file name.
#' @param workDir character string giving the path to working directory.
#' @return GRanges object with splice site coordinates.
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @keywords internal
#' @noRd
.loadSpliceSites <- function(x, workDir) {
    if (is.object(x = x)) {
        return(x)
    }
    frt <- tools::file_ext(x = x)
    if (!frt %in% "txt") {
        stop("Invalid file format")
    }
    spl.sites <- read.table(
        file = paste(workDir, x, sep = "/"),
        sep = "\t",
        header = TRUE,
        quote = "\"",
        as.is = TRUE
    )
    makeGRangesFromDataFrame(df = spl.sites, keep.extra.columns = TRUE)
}

#' Build the result object for introduceSNVtoSSs()
#' @description Internal helper for introduceSNVtoSSs(). Reorders the
#'     corrected splice site annotation columns and builds the named
#'     DNAStringSet with corrected sequences.
#' @param spl.sites GRanges object with corrected splice site sequences.
#' @return list with annotation and sequences.
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @keywords internal
#' @noRd
.buildIntroducedSpliceSiteResult <- function(spl.sites) {
    anno <- data.frame(spl.sites)
    anno <- anno[, c(
        seq_len(7),
        (ncol(x = anno) - 1):ncol(x = anno),
        seq(from = 8, to = ncol(x = anno) - 2)
    )]
    seqSSs <- DNAStringSet(x = anno$seq_corrected)
    names(x = seqSSs) <- paste(
        paste(seqnames(x = spl.sites),
            paste(start(x = spl.sites),
                end(x = spl.sites),
                sep = "-"
            ),
            sep = ":"
        ),
        strand(x = spl.sites),
        sep = "_str"
    )
    list(annotation = anno, sequences = seqSSs)
}

#' Load SNV ALT alleles from VCF
#' @description Internal helper for introduceSNVtoSSs(). Loads the ALT
#'     alleles from a VCF file and expands multi-allelic records.
#' @param vcfDir character string giving the path to VCF directory.
#' @param vcfFile character string giving the name of input VCF file.
#' @return GRanges object with ALT alleles.
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @keywords internal
#' @noRd
.loadSnvAltFromVcf <- function(vcfDir, vcfFile) {
    vcf <- rowRanges(x = readVcf(
        file = paste(vcfDir, vcfFile, sep = "/")
    ))[, "ALT"]
    ALT <- as.character(x = unlist(x = vcf$ALT))
    vcf <- vcf[rep(x = seq_along(vcf), times = lengths(x = vcf$ALT)), ]
    vcf$ALT <- ALT
    vcf
}

#' Correct splice site sequences
#' @description Internal helper for introduceSNVtoSSs(). Extracts the
#'     reference splice site sequences, substitutes the ALT allele at the
#'     SNV position, and returns the corrected sequences.
#' @param spl.sites GRanges object with splice site coordinates.
#' @param vcf GRanges object with SNV ALT alleles.
#' @param genome character string giving the name of the BSgenome package.
#' @param workDir character string giving the path to working directory.
#' @return GRanges object with corrected splice site sequences.
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @keywords internal
#' @noRd
.correctSpliceSiteSequences <- function(spl.sites, vcf, genome, workDir) {
    spl.sites$seq_original <- extractSeqSSs(
        x = spl.sites,
        genome = genome,
        workDir = workDir
    )
    spl.sites$seq_original <- as.character(x = spl.sites$seq_original)
    spl.sites$seq_corrected <- spl.sites$seq_original
    rev_compl <- as.character(x = strand(x = spl.sites)) == "-"
    spl.sites$seq_corrected[rev_compl] <-
        reverseComplement(
            x = DNAStringSet(x = spl.sites$seq_corrected[rev_compl])
        )
    if (length(x = spl.sites) != length(x = vcf)) {
        stop(
            "The list of splice sites to be corrected and the list ",
            "of relevant SNVs differ in length."
        )
    }
    substr(
        x = spl.sites$seq_corrected,
        start = start(x = vcf) - start(x = spl.sites) + 1,
        stop = start(x = vcf) - start(x = spl.sites) + 1
    ) <- vcf$ALT
    spl.sites$seq_corrected[rev_compl] <-
        reverseComplement(
            x = DNAStringSet(x = spl.sites$seq_corrected[rev_compl])
        )
    spl.sites
}


#' Score of splice sites
#' @description Scoring of the splice sites according to MaxEntScan
#'     models which are available at
#'     http://hollywood.mit.edu/burgelab/maxent or
#'     https://github.com/matthdsm/MaxEntScan.
#' @param x splice site sequences in FASTA format. It can be the name of
#'     object of class DNAStringSet (as output, for example, of
#'     extractSeqSSs() function) or standard FA/FASTA file.
#' @param type integer value 5 or 3 for define of donor or acceptor
#'     splice sites, respectively.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return object of class data frame containing three fields: IDs of
#'     splice sites, splice site sequences and MaxEntScan scores.
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @examples
#' if (nzchar(Sys.which("perl"))) {
#'     sites <- Biostrings::DNAStringSet(c(donor = "CAGGTAAGT"))
#'     local({
#'         old <- setwd(system.file("extdata", package = "huGSVcalleR"))
#'         on.exit(setwd(old))
#'         scoreSSs(x = sites, type = 5)
#'     })
#' }
#' @export
#' @importFrom Biostrings readDNAStringSet width uniqueLetters
#' @importFrom methods is
scoreSSs <- function(x, type, workDir = NULL) {
    if (is.null(x = workDir)) {
        workDir <- getwd()
    }
    seqSSs <- .loadSpliceSiteSequences(x, workDir)
    .checkPerlAvailable()
    .validateSpliceSiteSequences(seqSSs, type)
    scores <- .runMaxEntScan(seqSSs, type, workDir)
    .assembleScoreTable(seqSSs, scores)
}

#' Load splice site sequences
#' @description Internal helper for scoreSSs(). Loads splice site
#'     sequences from a DNAStringSet object or a FASTA file.
#' @param x DNAStringSet object or character string with the file name.
#' @param workDir character string with the working directory.
#' @return DNAStringSet object with splice site sequences.
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @keywords internal
#' @noRd
.loadSpliceSiteSequences <- function(x, workDir) {
    if (is.object(x = x)) {
        return(x)
    }
    frt <- tools::file_ext(x = x)
    if (!frt %in% c("fa", "fasta")) {
        stop("Invalid file format")
    }
    readDNAStringSet(
        filepath = paste(workDir, x, sep = "/"),
        format = frt
    )
}

#' Validate splice site sequences
#' @description Internal helper for scoreSSs(). Checks the splice site type
#'     and the length and content of the input sequences.
#' @param seqSSs DNAStringSet object with splice site sequences.
#' @param type integer value 5 or 3.
#' @return invisible(NULL). Stops on error.
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @keywords internal
#' @noRd
.validateSpliceSiteSequences <- function(seqSSs, type) {
    if (!type %in% c(3, 5)) {
        stop(
            "Invalid type of splice sites. ",
            "The entered type of splice sites must be numeric value 3 ",
            "or 5."
        )
    }
    if (((type == 3) & (unique(x = width(x = seqSSs))) != 23)) {
        stop(
            "Incorrect length of splice sites. ",
            "Acceptor splice site sequences must be 23 in length."
        )
    }
    if (((type == 5) & (unique(x = width(x = seqSSs))) != 9)) {
        stop(
            "Incorrect length of splice sites. ",
            "Donor splice site sequences must be 9 in length."
        )
    }
    if (!all(uniqueLetters(x = seqSSs) %in%
        c("a", "c", "g", "t", "G", "C", "T", "A"))) {
        warning(
            "Some sequences showed characters apart from A, C, G and T. ",
            "Scores can't be calculated for the incorrect records."
        )
    }
    invisible(NULL)
}

#' Run MaxEntScan Perl scripts
#' @description Internal helper for scoreSSs(). Writes the splice site
#'     sequences to a temporary file, runs the appropriate MaxEntScan Perl
#'     script, and parses the output.
#' @param seqSSs DNAStringSet object with splice site sequences.
#' @param type integer value 5 or 3.
#' @param workDir character string giving the path to working directory.
#' @return character vector with raw MaxEntScan scores.
#' @author Vasily V. Grinev, Nadzeya A. Boyeva.
#' @keywords internal
#' @noRd
.runMaxEntScan <- function(seqSSs, type, workDir) {
    script_name <- if (type == 3) "score3.pl" else "score5.pl"
    script_path <- system.file("extdata", script_name, package = "huGSVcalleR")
    if (!nzchar(script_path) || !file.exists(script_path)) {
        stop(
            "MaxEntScan script not found: ", script_name,
            ". Expected in inst/extdata/."
        )
    }
    tmp_file <- tempfile(pattern = "seqSSs_", fileext = ".txt")
    on.exit(unlink(tmp_file), add = TRUE)
    writeLines(as.character(seqSSs), con = tmp_file)
    scores <- system2(
        command = "perl",
        args = c(script_path, tmp_file),
        stdout = TRUE,
        stderr = TRUE
    )
    status <- attr(scores, "status")
    if (!is.null(status) && status != 0) {
        stop(
            "MaxEntScan failed with status ", status, ": ",
            paste(scores, collapse = "\n")
        )
    }
    scores <- scores[!grepl("^\\s*$", scores)]
    if (length(scores) == 0) {
        stop("MaxEntScan returned no scores.")
    }
    sub(pattern = "^[^\t]*\t", replacement = "", x = scores)
}
