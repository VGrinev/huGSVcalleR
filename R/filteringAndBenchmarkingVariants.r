#' Filter of single nucleotide variations
#' @description This function performs multi-stage filtering of single
#'     nucleotide variations given in the standard VCF format. It supports
#'     VCF files produced by all callers of the package: classical tests
#'     (binomial, counting, entropy, fisher, poisson), machine learning
#'     models (ml), Rsubread (fisherRSubread) and GATK. Filters specific
#'     to a particular caller are applied only if the corresponding INFO
#'     or GENO field is present in the VCF.
#' @param vcfDir character string giving the name (or path to and name) of
#'     directory with VCF file(-s). NULL by default, which means the current
#'     working directory will be used instead of the specified VCF directory.
#' @param vcfFile character string giving the name of input VCF file.
#' @param gr character string specifying the name of the tab-delimited TXT file
#'     containing coordinates of genomic location(-s) of interest. The default
#'     value is NULL. If so, this file must contains the following four fields:
#'     i) seqnames - the name of chromosome or scaffold with prefix "chr";
#'     ii) start   - start coordinate of genomic interval of interest;
#'     iii) end    - end coordinate of genomic interval of interest;
#'     iv) strand  - (optionally) strand information about genomic interval
#'                   of interest.
#'     TXT file with genomic location(-s) must be in working directory.
#' @param depth integer value giving the threshold for the minimum sequencing
#'     depth (DP field). Default value is 1.
#' @param score numeric value that specifies the threshold for the QUAL field.
#'     Default value is NULL.
#' @param heterozyg numeric value that specifies the threshold for the total
#'     heterozygosity calculated by the fields MMsum and DP. Default value is 0.
#' @param heterozyg1 numeric value that specifies the threshold for the
#'     heterozygosity of first alternative allele calculated by the fields MM
#'     and DP. Default value is 0.
#' @param vaf numeric value that specifies the threshold for the variant allele
#'     frequency calculated from the AD field (GENO). Default value is NULL.
#' @param qd numeric value that specifies the threshold for the QD field
#'     (Quality by Depth, GATK). Default value is NULL.
#' @param mq numeric value that specifies the threshold for the MQ field
#'     (Mapping Quality, GATK). Default value is NULL.
#' @param sor numeric value that specifies the threshold for the SOR field
#'     (Symmetric Odds Ratio, GATK). Default value is NULL.
#' @param fs numeric value that specifies the threshold for the FS field
#'     (Fisher Strand bias, GATK). Default value is NULL.
#' @param mqrs numeric value that specifies the threshold for the MQRankSum
#'     field (GATK). Default value is NULL.
#' @param rpps numeric value that specifies the threshold for the
#'     ReadPosRankSum field (GATK). Default value is NULL.
#' @param index logical argument, whether to bgzip the output file and generate
#'     a tabix index.
#' @param workDir character string giving the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return filtered VCF file with postfix "filtered".
#' @author Vasily V. Grinev.
#' @examples
#' ex <- system.file("extdata", package = "huGSVcalleR")
#' workDir <- tempdir()
#' file.copy(
#'     file.path(ex, "test_truth.vcf"),
#'     file.path(workDir, "test.vcf"),
#'     overwrite = TRUE
#' )
#' filterVariants(
#'     vcfDir = NULL,
#'     vcfFile = "test.vcf",
#'     depth = 5,
#'     workDir = workDir
#' )
#' @export
#' @importFrom VariantAnnotation readVcf writeVcf
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom S4Vectors queryHits
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom VariantAnnotation geno info fixed "fixed<-"
filterVariants <- function(vcfDir = NULL,
    vcfFile,
    gr = NULL,
    depth = 1,
    vaf = NULL,
    qd = NULL,
    mq = NULL,
    sor = NULL,
    fs = NULL,
    mqrs = NULL,
    rpps = NULL,
    score = NULL,
    heterozyg = 0,
    heterozyg1 = 0,
    index = FALSE,
    workDir = NULL) {
    ### Setting the working directory and VCF directory.
    paths <- .resolveVariantPaths(vcfDir, workDir)
    ### Loading of the single nucleotide variations.
    SNVs <- readVcf(file = paste(paths$vcfDir, vcfFile, sep = "/"))
    ### Detecting available INFO and GENO fields.
    info_fields <- colnames(info(x = SNVs))
    geno_fields <- rownames(geno(x = SNVs))
    ### Multi-stage filtering.
    SNVs <- .applyVariantFilters(
        SNVs = SNVs,
        gr = gr,
        depth = depth,
        vaf = vaf,
        qd = qd,
        mq = mq,
        sor = sor,
        fs = fs,
        mqrs = mqrs,
        rpps = rpps,
        score = score,
        heterozyg = heterozyg,
        heterozyg1 = heterozyg1,
        info_fields = info_fields,
        geno_fields = geno_fields,
        workDir = paths$workDir
    )
    ### Replace "." with "PASS" in the FILTER field.
    SNVs <- .replaceFilterPass(SNVs)
    ### Saving and returning of the final object.
    .writeFilteredVcf(SNVs, paths$vcfDir, vcfFile, index)
    return(SNVs)
}

#' Resolve working directory and VCF directory
#' @description Internal helper for filterVariants(). Returns the absolute
#'     paths to the working directory and the VCF directory.
#' @param vcfDir character string with the VCF directory, or NULL.
#' @param workDir character string with the working directory, or NULL.
#' @return list with workDir and vcfDir.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.resolveVariantPaths <- function(vcfDir, workDir) {
    if (is.null(x = workDir)) {
        workDir <- getwd()
    }
    list(
        workDir = workDir,
        vcfDir = if (is.null(x = vcfDir)) {
            workDir
        } else {
            paste(workDir, vcfDir, sep = "/")
        }
    )
}

#' Apply the requested variant filters
#' @description Internal helper for filterVariants(). Applies each requested
#'     filter in order, provided the required INFO or GENO field is present
#'     in the VCF.
#' @param SNVs CollapsedVCF object.
#' @param gr character string or NULL.
#' @param depth integer threshold.
#' @param vaf numeric threshold or NULL.
#' @param qd numeric threshold or NULL.
#' @param mq numeric threshold or NULL.
#' @param sor numeric threshold or NULL.
#' @param fs numeric threshold or NULL.
#' @param mqrs numeric threshold or NULL.
#' @param rpps numeric threshold or NULL.
#' @param score numeric threshold or NULL.
#' @param heterozyg numeric threshold.
#' @param heterozyg1 numeric threshold.
#' @param info_fields character vector with available INFO fields.
#' @param geno_fields character vector with available GENO fields.
#' @param workDir character string with the working directory.
#' @return filtered CollapsedVCF object.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.applyVariantFilters <- function(
      SNVs, gr, depth, vaf, qd, mq, sor, fs,
      mqrs, rpps, score, heterozyg, heterozyg1,
      info_fields, geno_fields, workDir
) {
    if (!is.null(x = gr)) {
        SNVs <- .filterByGenomicLocation(SNVs, gr, workDir)
    }
    if (!is.null(x = score)) {
        SNVs <- .filterByQual(SNVs, score)
    }
    if (depth > 1 && "DP" %in% info_fields) {
        SNVs <- .filterByDepth(SNVs, depth)
    }
    if (heterozyg > 0 && "MMsum" %in% info_fields) {
        SNVs <- .filterByTotalHeterozygosity(SNVs, heterozyg)
    }
    if (heterozyg1 > 0 && "MM" %in% info_fields) {
        SNVs <- .filterByFirstAlleleHeterozygosity(SNVs, heterozyg1)
    }
    if (!is.null(x = vaf) && "AD" %in% geno_fields) {
        SNVs <- .filterByVaf(SNVs, vaf)
    }
    SNVs <- .applyGatkInfoFilters(
        SNVs = SNVs,
        qd = qd, mq = mq, sor = sor, fs = fs, mqrs = mqrs, rpps = rpps,
        info_fields = info_fields
    )
    SNVs
}

#' Apply GATK-specific INFO filters
#' @description Internal helper for filterVariants(). Applies the
#'     greater-than and less-than INFO filters that are specific to GATK
#'     output, each only when the corresponding field is present.
#' @param SNVs CollapsedVCF object.
#' @param qd numeric threshold or NULL.
#' @param mq numeric threshold or NULL.
#' @param sor numeric threshold or NULL.
#' @param fs numeric threshold or NULL.
#' @param mqrs numeric threshold or NULL.
#' @param rpps numeric threshold or NULL.
#' @param info_fields character vector with available INFO fields.
#' @return filtered CollapsedVCF object.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.applyGatkInfoFilters <- function(
      SNVs, qd, mq, sor, fs, mqrs, rpps,
      info_fields
) {
    if (!is.null(x = qd) && "QD" %in% info_fields) {
        SNVs <- .filterByInfoGreater(SNVs, "QD", qd)
    }
    if (!is.null(x = mq) && "MQ" %in% info_fields) {
        SNVs <- .filterByInfoGreater(SNVs, "MQ", mq)
    }
    if (!is.null(x = sor) && "SOR" %in% info_fields) {
        SNVs <- .filterByInfoLess(SNVs, "SOR", sor)
    }
    if (!is.null(x = fs) && "FS" %in% info_fields) {
        SNVs <- .filterByInfoLess(SNVs, "FS", fs)
    }
    if (!is.null(x = mqrs) && "MQRankSum" %in% info_fields) {
        SNVs <- .filterByInfoGreater(SNVs, "MQRankSum", mqrs)
    }
    if (!is.null(x = rpps) && "ReadPosRankSum" %in% info_fields) {
        SNVs <- .filterByInfoGreater(SNVs, "ReadPosRankSum", rpps)
    }
    SNVs
}

#' Write the filtered VCF
#' @description Internal helper for filterVariants(). Resolves the output
#'     file name and writes the filtered VCF.
#' @param SNVs CollapsedVCF object.
#' @param vcfDir character string with the VCF directory.
#' @param vcfFile character string with the input VCF file name.
#' @param index logical, whether to bgzip and index the output.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.writeFilteredVcf <- function(SNVs, vcfDir, vcfFile, index) {
    writeVcf(
        obj = SNVs,
        filename = paste(
            vcfDir,
            sub(
                pattern = "\\.vcf$",
                replacement = ".filtered.vcf",
                x = vcfFile
            ),
            sep = "/"
        ),
        index = index
    )
    invisible(NULL)
}

#' Filter by genomic location
#' @description Internal helper for filterVariants(). Keeps only the
#'     variants that fall within the specified genomic intervals.
#' @param SNVs CollapsedVCF object.
#' @param gr character string with the intervals file name.
#' @param workDir character string with the working directory.
#' @return filtered CollapsedVCF object.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.filterByGenomicLocation <- function(SNVs, gr, workDir) {
    intervals <- read.table(
        file = paste(workDir, gr, sep = "/"),
        sep = "\t",
        header = TRUE,
        quote = "\"",
        as.is = TRUE
    )
    intervals <- makeGRangesFromDataFrame(df = intervals)
    hits <- findOverlaps(
        query = rowRanges(x = SNVs),
        subject = intervals,
        type = "within",
        ignore.strand = TRUE
    )
    SNVs[queryHits(x = hits)]
}

#' Filter by QUAL
#' @description Internal helper for filterVariants(). Keeps only the variants
#'     with QUAL >= score.
#' @param SNVs CollapsedVCF object.
#' @param score numeric threshold.
#' @return filtered CollapsedVCF object.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.filterByQual <- function(SNVs, score) {
    vals <- rowRanges(x = SNVs)$QUAL
    keep <- !is.na(vals) & vals >= score
    SNVs[keep]
}

#' Filter by sequencing depth
#' @description Internal helper for filterVariants(). Keeps only the variants
#'     with DP >= depth.
#' @param SNVs CollapsedVCF object.
#' @param depth integer threshold.
#' @return filtered CollapsedVCF object.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.filterByDepth <- function(SNVs, depth) {
    vals <- info(x = SNVs)$DP
    keep <- !is.na(vals) & vals >= depth
    SNVs[keep]
}

#' Filter by total heterozygosity
#' @description Internal helper for filterVariants(). Keeps only the variants
#'     with MMsum/DP >= heterozyg.
#' @param SNVs CollapsedVCF object.
#' @param heterozyg numeric threshold.
#' @return filtered CollapsedVCF object.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.filterByTotalHeterozygosity <- function(SNVs, heterozyg) {
    ratio <- info(x = SNVs)$MMsum / info(x = SNVs)$DP
    keep <- !is.na(ratio) & ratio >= heterozyg
    SNVs[keep]
}

#' Filter by first-allele heterozygosity
#' @description Internal helper for filterVariants(). Keeps only the variants
#'     with MM/DP >= heterozyg1.
#' @param SNVs CollapsedVCF object.
#' @param heterozyg1 numeric threshold.
#' @return filtered CollapsedVCF object.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.filterByFirstAlleleHeterozygosity <- function(SNVs, heterozyg1) {
    MM <- vapply(
        strsplit(as.character(info(x = SNVs)$MM), split = ","),
        function(y) max(as.numeric(y), na.rm = TRUE),
        numeric(1)
    )
    ratio <- MM / info(x = SNVs)$DP
    keep <- !is.na(ratio) & ratio >= heterozyg1
    SNVs[keep]
}

#' Filter by variant allele frequency
#' @description Internal helper for filterVariants(). Keeps only the variants
#'     with AD-based allele frequency >= vaf.
#' @param SNVs CollapsedVCF object.
#' @param vaf numeric threshold.
#' @return filtered CollapsedVCF object.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.filterByVaf <- function(SNVs, vaf) {
    ad <- geno(x = SNVs)$AD
    ad_mat <- do.call(rbind, ad)
    if (ncol(ad_mat) < 2) {
        return(SNVs)
    }
    ref <- ad_mat[, 1]
    alt <- ad_mat[, 2]
    af <- alt / (ref + alt)
    af[is.na(x = af)] <- 0
    keep <- !is.na(af) & af >= vaf
    SNVs[keep]
}

#' Filter by INFO field with greater-than comparison
#' @description Internal helper for filterVariants(). Keeps only the variants
#'     with INFO field value >= threshold.
#' @param SNVs CollapsedVCF object.
#' @param field character string with the INFO field name.
#' @param threshold numeric threshold.
#' @return filtered CollapsedVCF object.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.filterByInfoGreater <- function(SNVs, field, threshold) {
    vals <- info(x = SNVs)[[field]]
    keep <- !is.na(vals) & vals >= threshold
    SNVs[keep]
}

#' Filter by INFO field with less-than comparison
#' @description Internal helper for filterVariants(). Keeps only the variants
#'     with INFO field value <= threshold.
#' @param SNVs CollapsedVCF object.
#' @param field character string with the INFO field name.
#' @param threshold numeric threshold.
#' @return filtered CollapsedVCF object.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.filterByInfoLess <- function(SNVs, field, threshold) {
    vals <- info(x = SNVs)[[field]]
    keep <- !is.na(vals) & vals <= threshold
    SNVs[keep]
}

#' Replace "." with "PASS" in FILTER field
#' @description Internal helper for filterVariants(). Replaces the "." value
#'     in the FILTER field with "PASS".
#' @param SNVs CollapsedVCF object.
#' @return CollapsedVCF object with FILTER updated.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.replaceFilterPass <- function(SNVs) {
    fixed_SNVs <- fixed(x = SNVs)
    fixed_SNVs$FILTER[fixed_SNVs$FILTER == "."] <- "PASS"
    fixed(x = SNVs) <- fixed_SNVs
    SNVs
}


#' Benchmark called variants against reference variants
#' @description This function performs benchmarking of called variants against
#'    reference variants based on chromosome, position and allele identity.
#' @param vcfDir character string giving the name (or path to and name) of
#'     directory with VCF file(-s). NULL by default, which means the current
#'     working directory will be used instead of the specified VCF directory.
#' @param ref_vcf character string giving the name of VCF file with reference
#'    variants. Valid data format is *.vcf.
#' @param query_vcf character string giving the name of VCF file with called
#'    variants. Valid data format is *.vcf.
#' @param featureTable logical argument, whether to create the feature table
#'    that links each variant to its predicted classification outcome
#'    (TP - true positive, FP - false positive, FN - false negative). Default
#'    value is FALSE.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return consolidated data frame containing values of benchmarking metrics;
#'     in addition, if featureTable=TRUE, in the working directory is writing
#'     CSV file containing a detailed report of the benchmarking results.
#' @details
#'     The function splits variants into two classes (SNP and INDEL)
#'     based on the length of REF and ALT alleles and reports metrics
#'     separately for each class and for all variants together.
#'     \strong{Precision} is calculated as \code{TP / (TP + FP)}.
#'     \strong{Recall} is calculated as \code{TP / (TP + FN)}.
#'     \strong{F1} is the harmonic mean of precision and recall.
#'     When no variants of a given class were found by the query
#'     (\code{TP + FP = 0}), \strong{Precision} and \strong{F1} are
#'     reported as \code{NA}, because they are undefined in this case.
#'     \strong{Recall} is reported as 0 in this case, because all
#'     truth variants of this class were missed.
#' @author Liudmila S. Varaniuk.
#' @examples
#' ex <- system.file("extdata", package = "huGSVcalleR")
#' workDir <- tempdir()
#' file.copy(
#'     file.path(ex, "test_truth.vcf"),
#'     file.path(workDir, "ref.vcf"),
#'     overwrite = TRUE
#' )
#' file.copy(
#'     file.path(ex, "test_query.vcf"),
#'     file.path(workDir, "query.vcf"),
#'     overwrite = TRUE
#' )
#' benchmarkVariants(
#'     vcfDir = NULL,
#'     ref_vcf = "ref.vcf",
#'     query_vcf = "query.vcf",
#'     workDir = workDir
#' )
#' @export
#' @importFrom data.table fread as.data.table data.table fifelse
benchmarkVariants <- function(
      vcfDir,
      ref_vcf,
      query_vcf,
      featureTable = FALSE,
      workDir = NULL
) {
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
    ### Loading of VCF data.
    truth <- .readVcfAsDataTable(paste(vcfDir, ref_vcf, sep = "/"))
    test <- .readVcfAsDataTable(paste(vcfDir, query_vcf, sep = "/"))
    ### Developing a featureTable.
    if (isTRUE(x = featureTable)) {
        .writeBenchmarkFeatureTable(truth, test, query_vcf, ref_vcf, workDir)
    }
    ### Development and return of the final object.
    rbind(
        SNPs = .getBenchmarkMetrics(
            truth_keys = truth[type == "SNP", key],
            test_keys = test[type == "SNP", key]
        ),
        INDELs = .getBenchmarkMetrics(
            truth_keys = truth[type == "INDEL", key],
            test_keys = test[type == "INDEL", key]
        ),
        ALL = .getBenchmarkMetrics(truth_keys = truth$key, test_keys = test$key)
    )
}

#' Read VCF file into a data.table
#' @description Internal helper for benchmarkVariants(). Reads a VCF file into
#'     a data.table, splits multi-allelic ALT, classifies variant types, and
#'     builds a unique matching key.
#' @param input character string with the VCF file name.
#' @return data.table with the VCF records.
#' @author Liudmila S. Varaniuk.
#' @keywords internal
#' @noRd
.readVcfAsDataTable <- function(input) {
    if (is.character(x = input) && file.exists(x = input)) {
        if (file.size(input) == 0) {
            stop("VCF file is empty: ", input)
        }
        vcf <- fread(
            input,
            skip = "##",
            header = TRUE,
            sep = "\t",
            data.table = TRUE
        )
        if (ncol(vcf) < 8 || nrow(vcf) == 0) {
            stop("VCF file contains no variant records: ", input)
        }
        colnames(x = vcf)[seq_len(8)] <- c(
            "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"
        )
    } else if (is.data.frame(x = input)) {
        vcf <- as.data.table(x = input)
    } else {
        stop("Input data must be a path to a file or data frame.")
    }
    vcf <- vcf[, .(ALT = unlist(x = strsplit(x = ALT, split = ","))),
        by = .(CHROM, POS = as.integer(x = POS), REF)
    ]
    vcf[, type := ifelse(
        nchar(x = REF) == 1 & nchar(x = ALT) == 1,
        "SNP", "INDEL"
    )]
    vcf[, key := paste0(CHROM, ":", POS, "_", REF, ">", ALT)]
    vcf
}

#' Compute benchmark metrics
#' @description Internal helper for benchmarkVariants(). Computes TP, FP, FN,
#'     Precision, Recall and F1 for a given set of truth and test keys.
#' @param truth_keys character vector with truth keys.
#' @param test_keys character vector with test keys.
#' @return data.frame with the metrics.
#' @author Liudmila S. Varaniuk.
#' @keywords internal
#' @noRd
.getBenchmarkMetrics <- function(truth_keys, test_keys) {
    truth_keys <- unique(truth_keys)
    test_keys <- unique(test_keys)
    TP <- sum(test_keys %in% truth_keys)
    FP <- sum(!test_keys %in% truth_keys)
    FN <- sum(!truth_keys %in% test_keys)
    precision <- if ((TP + FP) == 0) NA else TP / (TP + FP)
    recall <- if ((TP + FN) == 0) NA else TP / (TP + FN)
    f1 <- if (is.na(x = precision) || is.na(x = recall) ||
        (precision + recall) == 0) {
        NA
    } else {
        2 * precision * recall / (precision + recall)
    }
    data.frame(
        TP = as.integer(x = TP),
        FP = as.integer(x = FP),
        FN = as.integer(x = FN),
        Precision = round(x = precision, digits = 4),
        Recall = round(x = recall, digits = 4),
        F1 = round(x = f1, digits = 4)
    )
}

#' Write benchmark feature table
#' @description Internal helper for benchmarkVariants(). Builds the feature
#'     table that links each variant to its predicted classification (TP, FP,
#'     or FN) and writes it as a CSV file.
#' @param truth data.table with truth variants.
#' @param test data.table with test variants.
#' @param query_vcf character string with the query VCF name.
#' @param ref_vcf character string with the reference VCF name.
#' @param workDir character string with the working directory.
#' @return invisible(NULL).
#' @author Liudmila S. Varaniuk.
#' @keywords internal
#' @noRd
.writeBenchmarkFeatureTable <- function(
      truth, test, query_vcf, ref_vcf,
      workDir
) {
    merged <- merge(
        truth[, .(key, CHROM, POS, REF, ALT, type)],
        test[, .(key, CHROM, POS, REF, ALT, type)],
        by = "key", all = TRUE, suffixes = c(".truth", ".test")
    )
    merged[, tag := fifelse(
        !is.na(x = CHROM.truth) & !is.na(x = CHROM.test), "TP",
        fifelse(is.na(x = CHROM.truth) & !is.na(x = CHROM.test), "FP", "FN")
    )]
    clean <- data.frame(
        CHROM = fifelse(
            is.na(merged$CHROM.test),
            merged$CHROM.truth,
            merged$CHROM.test
        ),
        POS = fifelse(
            is.na(merged$POS.test),
            merged$POS.truth,
            merged$POS.test
        ),
        REF = merged$REF.test,
        REF.truth = merged$REF.truth,
        ALT = merged$ALT.test,
        type = fifelse(
            is.na(merged$type.test),
            merged$type.truth,
            merged$type.test
        ),
        tag = merged$tag
    )
    outfile <- paste0(
        sub(pattern = "\\.vcf$", replacement = "", x = query_vcf),
        " vs ",
        sub(pattern = "\\.vcf$", replacement = "", x = ref_vcf),
        ", benchmarking results.csv"
    )
    write.csv(
        x = clean,
        file = paste(workDir, outfile, sep = "/"),
        row.names = FALSE
    )
    invisible(NULL)
}
