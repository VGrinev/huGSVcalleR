#' Local alignment of short DNA-Seq reads by seed-and-vote algorithm
#' @description This is a wrapper function that uses the low-level function
#'     align() from package Rsubread to align the short DNA-Seq reads by
#'     seed-and-vote algorithm.
#' @param genome character vector giving the path to and basename of index file.
#' @param fastqDir character vector giving the path to and name of directory
#'     with FASTQ file(-s).
#' @param fastq1 character vector including the name of FASTQ file containing
#'     DNA-Seq reads to be aligned. For paired-end reads, this gives the name
#'     of FASTQ file for first reads in each DNA-Seq library.
#' @param fastq2 character vector including the name of FASTQ file containing
#'     second DNA-Seq reads to be aligned. NULL by default.
#' @param bamDir character vector giving the path to and name of output
#'     directory for BAM file.
#' @param bamFile character vector giving the name of output BAM file.
#' @param orientation character string giving the orientation of the two reads
#'     from the same pair. Default value is "fr" (forward for the first read
#'     and reverse for second one).
#' @param threads numeric value giving the number of threads used for mapping.
#'     Default value is 1.
#' @param SV logical value indicating if structural variants will be detected
#'     during read mapping. Default value is FALSE.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @return BAM and VCF (if SV argument is TRUE) file(-s) are written to the
#'     bamDir directory.
#' @author Vasily V. Grinev
#' @examples
#' genome <- "Reference_Genomes/GRCh38/GRCh38"
#' fastq1 <- c("example_seq.read1.fastq.gz")
#' fastq2 <- c("example_seq.read1.fastq.gz")
#' fastqDir <- "Files_FASTQ"
#' bamDir <- "Files_BAM"
#' bam <- "example_seq.bam"
#' workDir <- "/mnt/data/grinev"
#' al <- alignDNASubread(
#'     genome = genome,
#'     fastqDir = fastqDir, fastq1 = fastq1, fastq2 = fastq2,
#'     bamDir = bamDir,
#'     bamFile = bam,
#'     orientation = "fr",
#'     threads = 28,
#'     SV = FALSE,
#'     workDir = workDir
#' )
#' @export
#' @importFrom Rsubread align
#' @importFrom Rsamtools indexBam
alignDNASubread <- function(
        genome,
        fastqDir, fastq1, fastq2 = NULL,
        bamDir,
        bamFile,
        orientation = "fr",
        threads = 1,
        SV = FALSE,
        workDir = NULL) {
    ### Assignment the basename of index file.
    ref <- paste(workDir, genome, sep = "/")
    ### Assignment the path to and name of FASTQ file(-s) including first reads.
    fq1 <- paste(paste(workDir, fastqDir, sep = "/"), fastq1, sep = "/")
    ### Assignment the path to and name of FASTQ file(-s) including second
    ### reads.
    if (is.null(x = fastq2) == FALSE) {
        fq2 <- paste(paste(workDir, fastqDir, sep = "/"), fastq2, sep = "/")
    } else {
        fq2 <- NULL
    }
    ### Alignment of short DNA-Seq reads against reference genome.
    align <- align(
        index = ref,
        readfile1 = fq1,
        readfile2 = fq2,
        type = "dna",
        input_format = "gzFASTQ",
        output_format = "BAM",
        output_file = paste(
            paste(workDir, bamDir, sep = "/"),
            bamFile,
            sep = "/"
        ),
        PE_orientation = orientation,
        nthreads = threads,
        unique = TRUE,
        detectSV = SV
    )
}


#' Evaluate alignment quality metrics
#' @description Evaluation of the alignment quality metrics for aligned DAN-Seq
#'     and/or RNA-Seq short reads.
#' @param bamDir character string giving the name (or path to and name) of
#'     directory with BAM file(s). NULL by default, which means the current
#'     working directory will be used instead of the specified BAM directory.
#' @param bamFile character vector giving the name of BAM file to be evaluated.
#' @param workDir character string giving the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return object of class list containing all main quality metrics for aligned
#'     short reads.
#' @author Vasily V. Grinev
#' @examples
#' qal <- assessQAlignedReads(
#'     bamDir = "Files_BAM",
#'     bamFile = "SRR21721515.bam",
#'     workDir = "D:/Vasily Grinev"
#' )
#' @export
#' @importFrom Rsamtools scanBam ScanBamParam bamFlagAsBitMatrix
#' @importFrom methods is
assessQAlignedReads <- function(bamDir = NULL,
    bamFile,
    workDir = NULL) {
    ### Setting the working directory.
    if (is.null(x = workDir)) {
        workDir <- getwd()
    }
    ### Full path to the BAM file.
    if (is.null(x = bamDir)) {
        path <- workDir
    } else {
        path <- paste(workDir, bamDir, sep = "/")
    }
    ### Collection of the quality metrics.
    name <- gsub(pattern = ".bam", replacement = "", x = bamFile)
    qa <- scanBam(
        file = paste(path, bamFile, sep = "/"),
        index = paste(paste(path, bamFile, sep = "/"), "bai", sep = "."),
        param = ScanBamParam(what = c("flag", "mapq", "isize"))
    )
    metrics <- .extractAlignmentMetrics(qa)
    rm("qa")
    ### Quality metrics.
    return(list(
        SN = name,
        FN = bamFile,
        TNRs = metrics$TNRs,
        MAPQs = metrics$MAPQs,
        TLENs = metrics$TLENs,
        FLAGs = metrics$FLAGs
    ))
}

#' Extract alignment metrics from scanBam output
#' @description Internal helper for assessQAlignedReads(). Extracts the total
#'     number of reads, mapping quality, template length, and flag bit counts
#'     from the scanBam() output.
#' @param qa list with scanBam() output.
#' @return list with TNRs, MAPQs, TLENs and FLAGs.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.extractAlignmentMetrics <- function(qa) {
    numbers <- length(x = qa[[1]]$flag)
    mapq <- qa[[1]]$mapq
    mapq <- mapq[!is.na(x = mapq)]
    isize <- abs(x = qa[[1]]$isize)
    isize <- isize[!is.na(x = isize)]
    flag <- colSums(x = bamFlagAsBitMatrix(flag = qa[[1]]$flag))
    list(
        TNRs = numbers,
        MAPQs = mapq,
        TLENs = isize,
        FLAGs = flag
    )
}


#' Build a hash table for the reference genome
#' @description This is a wrapper function that uses the low-level function
#'     buildindex() from package Subread to build the hash table for the
#'     reference genome.
#' @param ref_genome character string giving the name of folder for storing
#'     of index for reference genome.
#' @param ref_fasta character string giving the name of folder for storing
#'     of reference sequences.
#' @param index character string giving the basename of created index files.
#' @param fa charater string giving the name of the FASTA file containing all
#'     the reference sequences.
#' @param memory a numeric value specifying the amount of memory (in megabytes)
#'     used for storing the index during read mapping. Default value is 8000 MB.
#' @param workDir character string giving the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return no value is produced but index files are written to the respective
#'     directory.
#' @author Vasily V. Grinev
#' @examples
#' hash <- buildIndexSubread(
#'     ref_genome = "Reference_Genomes",
#'     ref_fasta = "Files_FASTA",
#'     index = "GRCh38",
#'     fa = "GRCh38.fa",
#'     memory = 8000,
#'     workDir = "/mnt/data/grinev"
#' )
#' @export
#' @importFrom Rsubread buildindex
buildIndexSubread <- function(ref_genome,
    ref_fasta,
    index,
    fa,
    memory = 8000,
    workDir = NULL) {
    ### Setting the working directory.
    if (is.null(x = workDir)) {
        workDir <- getwd()
    }
    ### Check for ref_genome directory existing.
    if (!dir.exists(file.path(workDir, ref_genome))) {
        dir.create(file.path(workDir, ref_genome), recursive = TRUE)
    }
    ### Setting the directories.
    refGenome <- paste(workDir, ref_genome, sep = "/")
    refFASTA <- paste(workDir, ref_fasta, sep = "/")
    ### Calculation of hash table.
    idx <- buildindex(
        basename = paste(refGenome, index, sep = "/"),
        reference = paste(refFASTA, fa, sep = "/"),
        gappedIndex = TRUE,
        indexSplit = TRUE,
        memory = memory,
        TH_subread = 24,
        colorspace = FALSE
    )
}


#' Collect and present alignment quality assessment results
#' @description A high-level function for consolidating and presenting the
#'     results of evaluation of the alignment quality metrics for aligned
#'     DAN-Seq and/or RNA-Seq short reads.
#' @param x a character string specifying the name of object developed by the
#'     function assessQAlignedReads().
#' @param output a character string that specifying a postfix for naming all
#'     output report files. The default value is NULL, at which the name of
#'     sample will be used for naming all output files.
#' @param workDir character string specifying the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return a series of text and graphic files.
#' @author Vasily V. Grinev
#' @examples
#' report <- reportQAlResults(
#'     x = qal,
#'     output = NULL,
#'     workDir = "D:/Vasily Grinev"
#' )
#' @export
#' @importFrom grid grid.draw textGrob gpar gTree gList
#' @importFrom gridExtra tableGrob ttheme_minimal
#' @importFrom vioplot vioplot
#' @importFrom grDevices cairo_pdf dev.off rgb
#' @importFrom graphics hist axis title plot rect text segments par mtext
#' @importFrom stats quantile median
reportQAlResults <- function(x,
    output = NULL,
    workDir = NULL) {
    ### Textual report.
    .saveAlTextReport(x, output, workDir)
    ### Graphical report.
    .saveAlGraphicalReport(x, output, workDir)
    ### Returning nothing.
    invisible(NULL)
}

#' Save alignment text report
#' @description Internal helper for reportQAlResults(). Builds the summary
#'     table and writes it as a PDF file.
#' @param x list with alignment metrics.
#' @param output character string with the output postfix, or NULL.
#' @param workDir character string with the working directory.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.saveAlTextReport <- function(x, output, workDir) {
    tab <- .buildAlTable(x)
    file_output <- .resolveAlFileName(
        "Mapping quality metrics of", x$SN, output
    )
    cairo_pdf(
        filename = paste(workDir, file_output, sep = "/"),
        height = nrow(x = tab) / 3 + 1
    )
    gt <- .buildAlTextReportGrob(tab, x$SN)
    grid.draw(gt)
    invisible(dev.off())
    invisible(NULL)
}

#' Build alignment text report grob
#' @description Internal helper for .saveAlTextReport(). Builds the gTree
#'     with the summary table and annotations.
#' @param tab data.frame with alignment metrics.
#' @param sample_name character string with the sample name.
#' @return gTree object.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.buildAlTextReportGrob <- function(tab, sample_name) {
    theme <- ttheme_minimal(
        core = list(
            fg_params = list(hjust = 0, x = 0.03, fontsize = 10),
            bg_params = list(
                fill = rep(c("azure1", "azure3"), each = 1),
                col = "black"
            )
        ),
        colhead = list(
            fg_params = list(col = "white"),
            bg_params = list(fill = "darkcyan", col = "black")
        )
    )
    table <- tableGrob(d = tab, rows = NULL, theme = theme)
    title1 <- textGrob(
        label = "Mapping quality. Summary statistics",
        just = "centre", vjust = -32.4,
        gp = gpar(col = "black", font = 2, fontsize = 12)
    )
    title2 <- textGrob(
        label = sprintf("Sample: %s", sample_name),
        just = "centre", vjust = -33.5,
        gp = gpar(col = "black", fontsize = 11)
    )
    title3 <- textGrob(
        label = sprintf(paste("Date & Time:",
            format(x = Sys.time(), "%Y-%m-%d %X"),
            sep = " "
        )),
        just = "centre", vjust = -34.9,
        gp = gpar(col = "black", fontsize = 10)
    )
    note <- textGrob(
        label = "Generated by assessQAlignedReads() & reportQAlResults()",
        just = "centre", vjust = 47,
        gp = gpar(col = "black", fontsize = 8)
    )
    gTree(children = gList(table, title1, title2, title3, note))
}

#' Build alignment summary table
#' @description Internal helper for reportQAlResults(). Builds the data.frame
#'     with the 25 mapping quality metrics.
#' @param x list with alignment metrics.
#' @return data.frame with two columns: Parameter and Value.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.buildAlTable <- function(x) {
    tab <- data.frame(
        Parameter = c(
            "Sample name", "File name",
            "Number of reads, total", "Number of reads, unique",
            "Number of reads, duplicated", "Number of reads, mapped",
            "Number of reads, unmapped", "Number of reads, paired",
            "Number of reads, properly paired",
            "Number of reads, not properly paired",
            "Number of reads, plus strand", "Number of reads, minus strand",
            "Number of reads, singletons", "Secondary alignments",
            "Supplementary (chimeric) alignments",
            "Mapping quality score, min",
            "Mapping quality score, 1st quartile",
            "Mapping quality score, median",
            "Mapping quality score, 3rd quartile",
            "Mapping quality score, max",
            "Template length, min", "Template length, 1st quartile",
            "Template length, median", "Template length, 3rd quartile",
            "Template length, max"
        ),
        Value = 0
    )
    tab[1, 2] <- x$SN
    tab[2, 2] <- x$FN
    tab[3, 2] <- x$TNRs
    tab[4, 2] <- x$TNRs - x$FLAGs["isDuplicate"][[1]]
    tab[5, 2] <- x$FLAGs["isDuplicate"][[1]]
    tab[6, 2] <- x$TNRs - x$FLAGs["isUnmappedQuery"][[1]]
    tab[7, 2] <- x$FLAGs["isUnmappedQuery"][[1]]
    tab[8, 2] <- x$FLAGs["isPaired"][[1]]
    tab[9, 2] <- x$FLAGs["isProperPair"][[1]]
    tab[10, 2] <- x$FLAGs["isPaired"][[1]] - x$FLAGs["isProperPair"][[1]]
    tab[11, 2] <- x$TNRs - x$FLAGs["isUnmappedQuery"][[1]] -
        x$FLAGs["isMinusStrand"][[1]]
    tab[12, 2] <- x$FLAGs["isMinusStrand"][[1]]
    tab[13, 2] <- x$FLAGs["hasUnmappedMate"][[1]]
    tab[14, 2] <- x$FLAGs["isSecondaryAlignment"][[1]]
    tab[15, 2] <- x$FLAGs["isSupplementaryAlignment"][[1]]
    tab[16:20, 2] <- quantile(
        x = x$MAPQs,
        probs = c(0.1, 0.25, 0.5, 0.75, 0.9)
    )
    tab[21:25, 2] <- quantile(
        x = x$TLENs,
        probs = c(0.1, 0.25, 0.5, 0.75, 0.9)
    )
    tab
}

#' Save graphical alignment report
#' @description Internal helper for reportQAlResults(). Writes all alignment
#'     quality plots to a single PDF file.
#' @param x list with alignment metrics.
#' @param output character string with the output postfix, or NULL.
#' @param workDir character string giving the path to working directory.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.saveAlGraphicalReport <- function(x, output, workDir) {
    file_output <- .resolveAlFileName(
        "Mapping quality plots of", x$SN, output
    )
    cairo_pdf(
        filename = paste(workDir, file_output, sep = "/"),
        onefile = TRUE, width = 8.27, height = 11.69
    )
    par(mfrow = c(3, 2), omi = c(0.79, 0.79, 1, 0.79))
    .plotMappingQualityDistribution(x)
    .plotMappingQualityThreshold(x)
    .plotTemplateLengthDistribution(x)
    .plotTemplateLengthThreshold(x)
    .addAlReportHeader(x)
    invisible(dev.off())
    invisible(NULL)
}

#' Resolve alignment output file name
#' @description Internal helper for reportQAlResults(). Builds the output file
#'     name based on the prefix, sample name or output, and extension.
#' @param prefix character string with the report prefix.
#' @param sample_name character string with the sample name.
#' @param output character string with the output postfix, or NULL.
#' @return character string with the output file name.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.resolveAlFileName <- function(prefix, sample_name, output) {
    if (is.null(x = output)) {
        paste(sprintf("%s %s", prefix, sample_name), "pdf", sep = ".")
    } else {
        paste(sprintf("%s %s", prefix, output), "pdf", sep = ".")
    }
}

#' Plot mapping quality distribution
#' @description Internal helper for reportQAlResults(). Draws the histogram
#'     of mapping quality.
#' @param x list with alignment metrics.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.plotMappingQualityDistribution <- function(x) {
    hist(
        x = x$MAPQs, axes = FALSE, ann = FALSE, border = FALSE,
        xlim = c(0, 40),
        col = rgb(red = 0, green = 158, blue = 115, maxColorValue = 255)
    )
    axis(
        side = 1, at = c(0, 10, 20, 30, 40),
        labels = c(0, 10, 20, 30, 40), las = 1, lwd = 0.8
    )
    axis(
        side = 2,
        at = c(0.0, max(x = table(x$MAPQs)) / 3 * c(1, 2, 3)),
        labels = formatC(
            x = c(0.0, max(x = table(x$MAPQs)) / 3 * c(1, 2, 3)),
            format = "e", digits = 2
        ), lwd = 0.8
    )
    title(main = "Distribution of mapping qualities", cex.main = 1.5)
    title(xlab = "Mapping quality scores", line = 2.5, cex.lab = 1.2)
    title(ylab = "Number of reads", line = 3, cex.lab = 1.2)
    invisible(NULL)
}

#' Plot mapping quality threshold
#' @description Internal helper for reportQAlResults(). Draws the threshold
#'     plot of per-read mapping quality.
#' @param x list with alignment metrics.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.plotMappingQualityThreshold <- function(x) {
    per <- .computeMappingQualityThresholdTable(x$MAPQs)
    .drawMappingQualityThresholdPlot(per)
    .drawMappingQualityThresholdBands()
    .annotateMappingQualityThreshold(x$MAPQs)
    invisible(NULL)
}

#' Compute mapping quality threshold table
#' @description Internal helper for .plotMappingQualityThreshold().
#'     Computes the percentage of reads at or above each MAPQ threshold.
#' @param mapqs numeric vector of per-read mapping qualities.
#' @return data.frame with threshold and percent columns.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.computeMappingQualityThresholdTable <- function(mapqs) {
    per <- data.frame(threshold = 0:40, percent = 0)
    per[1, 2] <- 100
    L <- length(x = mapqs)
    for (i in seq_len(40)) {
        per[i + 1, 2] <- length(x = mapqs[mapqs >= i]) / L * 100
    }
    per
}

#' Draw the mapping quality threshold curve
#' @description Internal helper for .plotMappingQualityThreshold(). Draws
#'     the curve and the axes.
#' @param per data.frame with threshold and percent columns.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.drawMappingQualityThresholdPlot <- function(per) {
    plot(per,
        type = "l", bty = "n",
        xlim = c(0, 40), ylim = c(0, 100), lwd = 1,
        col = rgb(red = 0, green = 114, blue = 178, maxColorValue = 255),
        ann = FALSE, axes = FALSE
    )
    axis(
        side = 1, at = c(0, 10, 20, 30, 40),
        labels = c(0, 10, 20, 30, 40), las = 1, lwd = 0.8
    )
    axis(
        side = 2, at = c(0, 20, 40, 60, 80, 100),
        labels = c(0, 20, 40, 60, 80, 100), las = 1, lwd = 0.8
    )
    title(main = "Thresholding of mapping quality", cex.main = 1.5)
    title(xlab = "Quality threshold", line = 2.5, cex.lab = 1.2)
    title(ylab = "Reads exceeding quality level, %", line = 2.8, cex.lab = 1.2)
    invisible(NULL)
}

#' Draw the mapping quality threshold bands
#' @description Internal helper for .plotMappingQualityThreshold(). Draws
#'     the coloured background bands.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.drawMappingQualityThresholdBands <- function() {
    rect(
        xleft = 0, ybottom = 0, xright = 10, ytop = 100, border = NA,
        col = rgb(240, 228, 66, maxColorValue = 255, alpha = 30)
    )
    rect(
        xleft = 10, ybottom = 0, xright = 20, ytop = 100, border = NA,
        col = rgb(230, 159, 0, maxColorValue = 255, alpha = 40)
    )
    rect(
        xleft = 20, ybottom = 0, xright = 30, ytop = 100, border = NA,
        col = rgb(213, 94, 0, maxColorValue = 255, alpha = 50)
    )
    rect(
        xleft = 30, ybottom = 0, xright = 40, ytop = 100, border = NA,
        col = rgb(213, 94, 0, maxColorValue = 255, alpha = 80)
    )
    invisible(NULL)
}

#' Annotate the mapping quality threshold bands
#' @description Internal helper for .plotMappingQualityThreshold(). Draws
#'     the MAPQ10, MAPQ20 and MAPQ30 labels.
#' @param mapqs numeric vector of per-read mapping qualities.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.annotateMappingQualityThreshold <- function(mapqs) {
    L <- length(x = mapqs)
    for (t in c(10, 20, 30)) {
        text(
            x = t, y = 23, srt = 90, cex = 0.9,
            labels = paste0(
                "MAPQ", t, ", ",
                round(length(mapqs[mapqs >= t]) / L * 100, 1),
                "%"
            )
        )
    }
    invisible(NULL)
}

#' Plot template length distribution
#' @description Internal helper for reportQAlResults(). Draws the violin
#'     plot of template lengths.
#' @param x list with alignment metrics.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.plotTemplateLengthDistribution <- function(x) {
    if (x$TNRs > 1e6) {
        p <- sample(x = log10(x = x$TLENs + 1), size = 1e6)
    } else {
        p <- log10(x = x$TLENs + 1)
    }
    plot(0,
        type = "n", axes = FALSE, ann = FALSE, bty = "n",
        xlim = c(0, 2), ylim = c(0, ceiling(x = max(x = p)))
    )
    vioplot(
        x = p, add = TRUE,
        col = rgb(86, 180, 233, maxColorValue = 255),
        rectCol = rgb(0, 158, 115, maxColorValue = 255),
        pchMed = 0, colMed = NA, lwd = 0.9, wex = 1, frame.plot = FALSE
    )
    segments(
        x0 = 0.9, y0 = median(x = p),
        x1 = 1.1, y1 = median(x = p),
        col = rgb(240, 228, 66, maxColorValue = 255),
        lwd = 3, lend = 2
    )
    axis(side = 1, at = c(0, 1, 2), labels = FALSE, lwd = 0.8)
    axis(
        side = 2,
        at = c(0.0, ceiling(x = max(x = p)) / 4 * seq_len(4)),
        labels = round(
            c(
                0.0,
                ceiling(x = max(x = p)) / 4 * seq_len(4)
            ),
            digits = 2
        ),
        las = 1, lwd = 0.8
    )
    title(main = "Distribution of template lengths", cex.main = 1.5)
    title(xlab = "Standard vioplot", line = 1, cex.lab = 1.2)
    title(
        ylab = expression("Length, log"[10] * " (nucleotides)"),
        line = 2.8, cex.lab = 1.2
    )
    invisible(NULL)
}

#' Plot template length threshold
#' @description Internal helper for reportQAlResults(). Draws the threshold
#'     plot of log10-transformed template lengths.
#' @param x list with alignment metrics.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.plotTemplateLengthThreshold <- function(x) {
    p <- log10(x = x$TLENs + 1)
    per <- .computeTemplateLengthThresholdTable(p)
    .drawTemplateLengthThresholdPlot(per, ceiling(x = max(x = p)))
    .drawTemplateLengthThresholdBands(ceiling(x = max(x = p)))
    invisible(NULL)
}

#' Compute template length threshold table
#' @description Internal helper for .plotTemplateLengthThreshold(). Computes
#'     the percentage of reads at or above each log10 template length
#'     threshold.
#' @param p numeric vector of log10-transformed template lengths.
#' @return data.frame with threshold and percent columns.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.computeTemplateLengthThresholdTable <- function(p) {
    tr <- c(0.0, ceiling(x = max(x = p)) / 100 * seq_len(100))
    per <- data.frame(threshold = tr, percent = 0)
    per[1, 2] <- 100
    L <- length(x = p)
    for (i in seq_len(length(x = tr) - 1)) {
        per[i + 1, 2] <- length(x = p[p >= tr[i]]) / L * 100
    }
    per
}

#' Draw the template length threshold curve
#' @description Internal helper for .plotTemplateLengthThreshold(). Draws
#'     the curve, axes and titles.
#' @param per data.frame with threshold and percent columns.
#' @param max_log10 numeric, ceiling of the maximum log10 template length.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.drawTemplateLengthThresholdPlot <- function(per, max_log10) {
    plot(per,
        type = "l", bty = "n",
        xlim = c(0, max_log10), ylim = c(0, 100), lwd = 1,
        col = rgb(red = 0, green = 114, blue = 178, maxColorValue = 255),
        ann = FALSE, axes = FALSE
    )
    axis(
        side = 1,
        at = c(0.0, max_log10 / 4 * seq_len(4)),
        labels = c(0.0, max_log10 / 4 * seq_len(4)),
        las = 1, lwd = 0.8
    )
    axis(
        side = 2, at = c(0, 20, 40, 60, 80, 100),
        labels = c(0, 20, 40, 60, 80, 100), las = 1, lwd = 0.8
    )
    title(main = "Thresholding of template lengths", cex.main = 1.5)
    title(
        xlab = expression("Length, log"[10] * " (nucleotides)"),
        line = 2.5, cex.lab = 1.2
    )
    title(ylab = "Templates exceeding length, %", line = 2.8, cex.lab = 1.2)
    invisible(NULL)
}

#' Draw the template length threshold bands
#' @description Internal helper for .plotTemplateLengthThreshold(). Draws
#'     the five horizontal coloured background bands.
#' @param max_log10 numeric, ceiling of the maximum log10 template length.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.drawTemplateLengthThresholdBands <- function(max_log10) {
    bounds <- list(
        list(y0 = 0, y1 = 20, alpha = 30),
        list(y0 = 20, y1 = 40, alpha = 40),
        list(y0 = 40, y1 = 60, alpha = 50),
        list(y0 = 60, y1 = 80, alpha = 70),
        list(y0 = 80, y1 = 100, alpha = 90)
    )
    for (b in bounds) {
        rect(
            xleft = 0, ybottom = b$y0,
            xright = max_log10, ytop = b$y1, border = NA,
            col = rgb(213, 94, 0, maxColorValue = 255, alpha = b$alpha)
        )
    }
    invisible(NULL)
}

#' Add alignment report header
#' @description Internal helper for reportQAlResults(). Adds the common
#'     header (title, sample name, date & time, footnote) to the current
#'     graphical device.
#' @param x list with alignment metrics.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.addAlReportHeader <- function(x) {
    mtext(
        text = "Mapping quality. Summary plots",
        side = 3, line = 3.5, cex = 1.1,
        col = "black", font = 2, outer = TRUE
    )
    mtext(
        text = sprintf("Sample: %s", x$SN),
        side = 3, line = 2, cex = 0.95,
        col = "black", outer = TRUE
    )
    mtext(
        text = sprintf(paste("Date & Time:",
            format(x = Sys.time(), "%Y-%m-%d %X"),
            sep = " "
        )),
        side = 3, line = 0.5, cex = 0.8,
        col = "black", outer = TRUE
    )
    note <- "Generated by assessQAlignedReads() & reportQAlResults()"
    mtext(
        text = sprintf(note),
        side = 1, line = 2, cex = 0.7,
        col = "black", outer = TRUE
    )
    invisible(NULL)
}


#' Sort and index BAM file(-s)
#' @description This is a wrapper function that uses the low-level function
#'     sortBam() from package Rsamtools to sort and index BAM file(-s).
#' @param bamDir character vector giving the path to and name of directory
#'     with BAM file(s).
#' @param bamFile character vector giving the name(-s) of BAM file(-s) to be
#'     sorted and indexed.
#' @param byQname logical. It indicats whether the sorted destination file
#'     should be sorted by query-name (TRUE) or by mapping position (FALSE).
#' @param threads numeric value giving the number of threads used for sorting.
#'     Default is one less than the number of logical cores, but at least 1.
#'     If the core count cannot be detected, one thread is used.
#' @param workDir character string giving the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return sorted BAM file(-s). Position-sorted files are also indexed.
#'     Name-sorted files are not indexed, because indexBam() rejects them
#' @author Vasily V. Grinev
#' @examples
#' sortBamFile(
#'     bamDir = "Files_BAM",
#'     bamFile = "SRR21721515.bam",
#'     byQname = TRUE,
#'     threads = 10,
#'     workDir = "D:/Vasily Grinev"
#' )
#' @export
#' @importFrom Rsamtools sortBam indexBam
#' @importFrom parallel detectCores
sortBamFile <- function(bamDir,
    bamFile,
    byQname,
    threads = max(1L, detectCores(logical = TRUE) - 1L, na.rm = TRUE),
    workDir = NULL) {
    ### Setting the working directory.
    if (is.null(x = workDir)) {
        workDir <- getwd()
    }
    ### Full path to the BAM file.
    if (is.null(x = bamDir)) {
        path <- workDir
    } else {
        path <- paste(workDir, bamDir, sep = "/")
    }
    ### Sorting and indexing of generated BAM file(-s).
    for (i in seq_along(bamFile)) {
        sortedBam <- sortBam(
            file = paste(path, bamFile[i], sep = "/"),
            destination = sub(
                "\\.bam$", ".sorted",
                paste(path, bamFile[i], sep = "/")
            ),
            byQname = byQname,
            nThreads = threads
        )
        if (!isTRUE(x = byQname)) {
            invisible(utils::capture.output(
                indexBam(files = sortedBam),
                type = "message"
            ))
        }
    }
}
