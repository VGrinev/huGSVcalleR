#' Extract chromosome and position from a temporary pileup file name
#' @description Internal helper function for parsing temporary pileup file
#'     names of the form "temp_stat_<chr>_<pos>.csv".
#' @param filename character string giving the name of a temporary file.
#' @return a list with two elements: \code{chr} (numeric chromosome code,
#'     where X = 23, Y = 24) and \code{pos} (numeric position).
#' @author Vasily V. Grinev, Dzianis D. Sarnatski, Mikalai M. Yatskou
#' @keywords internal
#' @noRd
.extractChrPos <- function(filename) {
    matches <- str_match(
        string = filename,
        pattern = "temp_stat_chr(\\d+|[XY])_(\\d+)\\.csv"
    )
    chr <- matches[, 2]
    pos <- as.numeric(matches[, 3])
    chr <- ifelse(chr == "X", 23,
        ifelse(chr == "Y", 24, as.numeric(x = chr))
    )
    list(chr = chr, pos = pos)
}


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
#'     Default value is tempdir(), which means the R session temporary
#'     directory will be used.
#' @param postfix character string that will be added to the CSV output file
#'     name as postfix. NULL by default. In this case, the postfix
#'     "PileupResults" will be used.
#' @param threads integer value giving the number of threads (workers) used
#'     for parallel pileup calculation. Default value is 1.
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
#' res <- collectBasesPileup(
#'     bamDir = "Files_BAM",
#'     bamFile = "SRR21721515.bam",
#'     gr = "genomic_intervals.txt",
#'     canChr = TRUE,
#'     depth = 1e6,
#'     baseq = 20,
#'     mapq = 20,
#'     coverage = 10,
#'     genome = "BSgenome.Hsapiens.UCSC.hg38",
#'     tmpDir = tempdir(),
#'     postfix = "PileupResults",
#'     threads = 2,
#'     workDir = "D:/Vasily Grinev"
#' )
#' @export
#' @importFrom data.table data.table fwrite as.data.table
#' @importFrom BiocParallel bplapply SnowParam
#' @importFrom Rsamtools BamFile ScanBamParam seqinfo pileup PileupParam
#' @importFrom GenomicRanges makeGRangesFromDataFrame GRanges
#' @importFrom BSgenome getSeq
#' @importFrom dplyr mutate case_when filter arrange pull select across
#' @importFrom stringr str_match
#' @importFrom purrr map_dfr
#' @importFrom tidyr %>%
#' @importFrom utils read.table
#' @importFrom methods is
#' @importFrom tibble tibble
collectBasesPileup <- function(
        bamDir = NULL,
        bamFile,
        gr = "SeqInfo",
        canChr = TRUE,
        depth = 1e6, baseq = 20, mapq = 20,
        coverage = NULL,
        genome = "BSgenome.Hsapiens.UCSC.hg38",
        tmpDir = tempdir(),
        postfix = NULL,
        threads = 1,
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
    ### Setting the full path to the BAM file and reference to this file.
    bam <- paste(path, bamFile, sep = "/")
    bam_link <- BamFile(file = bam)
    ### Setting the genomic location(-s) of interest.
    GRs <- .buildPileupRegions(gr, canChr, bam_link, workDir)
    ### Loading the reference genome object.
    genome_obj <- .loadReferenceGenome(genome)
    ### Collection of pileup statistics.
    .runPileupInParallel(
        GRs = GRs,
        bam_link = bam_link,
        depth = depth,
        baseq = baseq,
        mapq = mapq,
        coverage = coverage,
        genome_obj = genome_obj,
        tmpDir = tmpDir,
        threads = threads
    )
    ### Consolidation of pileup statistics.
    .consolidatePileupResults(tmpDir, bam, postfix)
}

#' Load a reference genome package
#' @description Internal helper for collectBasesPileup(). Checks that the
#'     BSgenome package is installed and returns its namespace object.
#' @param genome character string with the BSgenome package name.
#' @return the BSgenome object.
#' @author Vasily V. Grinev, Dzianis D. Sarnatski, Mikalai M. Yatskou.
#' @keywords internal
#' @noRd
.loadReferenceGenome <- function(genome) {
    if (!requireNamespace(genome, quietly = TRUE)) {
        stop(
            "Package ", genome, " is not installed. ",
            "Install it with BiocManager::install('", genome, "')."
        )
    }
    get(x = genome, envir = asNamespace(genome))
}

#' Run the pileup computation in parallel
#' @description Internal helper for collectBasesPileup(). Builds the
#'     BiocParallel parameter object and calls .computePileupChunk() for
#'     each genomic chunk.
#' @param GRs GRanges object with genomic regions.
#' @param bam_link BamFile object with the BAM file.
#' @param depth integer, maximum pileup depth.
#' @param baseq integer, minimal base quality.
#' @param mapq integer, minimal mapping quality.
#' @param coverage integer, minimal coverage, or NULL.
#' @param genome_obj BSgenome object with the reference genome.
#' @param tmpDir character string with the temporary directory.
#' @param threads integer, number of workers.
#' @return list returned by bplapply().
#' @author Vasily V. Grinev, Dzianis D. Sarnatski, Mikalai M. Yatskou.
#' @keywords internal
#' @noRd
.runPileupInParallel <- function(GRs, bam_link, depth, baseq, mapq,
    coverage, genome_obj, tmpDir, threads) {
    bpparam <- SnowParam(workers = threads, stop.on.error = TRUE)
    bplapply(
        X = seq_along(GRs),
        FUN = .computePileupChunk,
        GRs = GRs,
        bam_link = bam_link,
        depth = depth,
        baseq = baseq,
        mapq = mapq,
        coverage = coverage,
        genome_obj = genome_obj,
        tmpDir = tmpDir,
        BPPARAM = bpparam
    )
}

#' Build pileup regions
#' @description Internal helper for collectBasesPileup(). Builds the GRanges
#'     object with genomic regions for pileup calculation, either from the
#'     BAM file header (gr = "SeqInfo") or from a user-provided TXT file.
#' @param gr character string with the method or file name.
#' @param canChr logical, remove non-canonical chromosomes?
#' @param bam_link BamFile object with the BAM file.
#' @param workDir character string giving the path to working directory.
#' @return GRanges object with genomic regions.
#' @author Vasily V. Grinev, Dzianis D. Sarnatski, Mikalai M. Yatskou.
#' @keywords internal
#' @noRd
.buildPileupRegions <- function(gr, canChr, bam_link, workDir) {
    if (gr == "SeqInfo") {
        GRs <- as.data.frame(cbind(
            seqnames(seqinfo(x = bam_link)), 1,
            seqlengths(seqinfo(x = bam_link))
        ))
        colnames(x = GRs) <- c("seqnames", "start", "end")
        GRs$start <- as.numeric(x = GRs$start)
        GRs$end <- as.numeric(x = GRs$end)
        if (isTRUE(x = canChr)) {
            GRs <- .filterCanonicalChromosomes(GRs, bam_link)
        }
        rownames(x = GRs) <- NULL
        makeGRangesFromDataFrame(df = GRs)
    } else {
        GRs <- read.table(
            file = paste(workDir, gr, sep = "/"),
            sep = "\t",
            header = TRUE,
            quote = "\"",
            as.is = TRUE
        )
        if (isTRUE(x = canChr)) {
            GRs <- .filterCanonicalChromosomes(GRs, bam_link)
        }
        if (grepl(pattern = "chr", x = seqnames(seqinfo(x = bam_link))[1])) {
            makeGRangesFromDataFrame(df = GRs)
        } else {
            GRs$seqnames <- gsub(
                pattern = "chr",
                replacement = "",
                x = GRs$seqnames
            )
            makeGRangesFromDataFrame(df = GRs)
        }
    }
}

#' Filter canonical chromosomes
#' @description Internal helper for collectBasesPileup(). Keeps only the
#'     records from canonical chromosomes (chr1-chr22, chrX, chrY, chrM).
#' @param GRs data.frame with genomic regions.
#' @param bam_link BamFile object with the BAM file.
#' @return data.frame with canonical chromosomes only.
#' @author Vasily V. Grinev, Dzianis D. Sarnatski, Mikalai M. Yatskou.
#' @keywords internal
#' @noRd
.filterCanonicalChromosomes <- function(GRs, bam_link) {
    if (grepl(pattern = "chr", x = seqnames(seqinfo(x = bam_link))[1])) {
        GRs[nchar(x = GRs$seqnames) <= 5, ]
    } else {
        GRs[nchar(x = GRs$seqnames) == 1, ]
    }
}

#' Compute pileup for a single genomic chunk
#' @description Internal helper for collectBasesPileup(). Computes the
#'     nucleotide-wise coverage for a single genomic chunk and writes the
#'     result to a temporary CSV file.
#' @param i integer index of the genomic chunk.
#' @param GRs GRanges object with genomic regions.
#' @param bam_link BamFile object with the BAM file.
#' @param depth integer, maximum pileup depth.
#' @param baseq integer, minimal base quality.
#' @param mapq integer, minimal mapping quality.
#' @param coverage integer, minimal coverage, or NULL.
#' @param genome_obj BSgenome object with the reference genome.
#' @param tmpDir character string with the temporary directory.
#' @return NULL.
#' @author Vasily V. Grinev, Dzianis D. Sarnatski, Mikalai M. Yatskou.
#' @keywords internal
#' @noRd
.computePileupChunk <- function(i, GRs, bam_link, depth, baseq, mapq,
    coverage, genome_obj, tmpDir) {
    chunk <- GRs[i]
    pileup_data <- pileup(
        file = bam_link,
        scanBamParam = ScanBamParam(which = chunk),
        pileupParam = PileupParam(
            max_depth = depth,
            min_base_quality = baseq,
            min_mapq = mapq,
            distinguish_strands = FALSE,
            include_deletions = FALSE,
            include_insertions = FALSE
        )
    )
    pileup_data$position <- paste(
        pileup_data$seqnames, pileup_data$pos,
        sep = ":"
    )
    temp_stat <- .summarizePileupCounts(pileup_data)
    if (!is.null(x = coverage)) {
        temp_stat <- temp_stat %>%
            filter(rowSums(x = across(c(A, C, G, .data$T))) >= coverage)
    }
    temp_stat <- .assignReferenceNucleotides(temp_stat, genome_obj)
    if (nrow(x = temp_stat) > 0) {
        .writePileupChunk(temp_stat, tmpDir)
    }
    NULL
}

#' Summarize pileup counts by position
#' @description Internal helper for collectBasesPileup(). Aggregates the raw
#'     pileup records into a six-column table per position.
#' @param pileup_data data.frame from pileup().
#' @return data.table with summarized counts.
#' @author Vasily V. Grinev, Dzianis D. Sarnatski, Mikalai M. Yatskou.
#' @keywords internal
#' @noRd
.summarizePileupCounts <- function(pileup_data) {
    temp_stat <- pileup_data %>%
        mutate(reference = 0, A = 0, C = 0, G = 0, T = 0) %>%
        mutate(
            A = case_when(nucleotide == "A" ~ count, TRUE ~ .data$A),
            C = case_when(nucleotide == "C" ~ count, TRUE ~ .data$C),
            G = case_when(nucleotide == "G" ~ count, TRUE ~ .data$G),
            T = case_when(nucleotide == "T" ~ count, TRUE ~ .data$T)
        )
    temp_stat <- as.data.table(temp_stat)
    temp_stat <- temp_stat[,
        .(
            seqnames = seqnames,
            pos = pos,
            reference = sum(x = reference),
            A = sum(x = A),
            C = sum(x = C),
            G = sum(x = G),
            T = sum(x = get("T"))
        ),
        by = position
    ]
    temp_stat[, .(position, seqnames, pos, reference, A, C, G, get("T"))]
}

#' Assign reference nucleotides to pileup positions
#' @description Internal helper for collectBasesPileup(). Retrieves the
#'     reference nucleotide for each position from the BSgenome object.
#' @param temp_stat data.table with pileup statistics.
#' @param genome_obj BSgenome object with the reference genome.
#' @return data.table with reference nucleotides.
#' @author Vasily V. Grinev, Dzianis D. Sarnatski, Mikalai M. Yatskou.
#' @keywords internal
#' @noRd
.assignReferenceNucleotides <- function(temp_stat, genome_obj) {
    ref <- temp_stat %>%
        select(seqnames, pos) %>%
        mutate(start = pos) %>%
        mutate(end = pos) %>%
        select(-pos)
    ref <- makeGRangesFromDataFrame(df = ref)
    temp_stat %>%
        mutate(
            reference = getSeq(
                x = genome_obj,
                names = ref,
                as.character = TRUE
            )
        ) %>%
        select(-pos, -seqnames)
}

#' Write pileup chunk to a temporary file
#' @description Internal helper for collectBasesPileup(). Writes the pileup
#'     statistics for a single chunk to a temporary CSV file.
#' @param temp_stat data.table with pileup statistics.
#' @param tmpDir character string with the temporary directory.
#' @return invisible(NULL).
#' @author Vasily V. Grinev, Dzianis D. Sarnatski, Mikalai M. Yatskou.
#' @keywords internal
#' @noRd
.writePileupChunk <- function(temp_stat, tmpDir) {
    position_label <- temp_stat[1, "position"]
    position_label <- gsub(pattern = ":", replacement = "_", x = position_label)
    tmp_file <- file.path(
        tmpDir,
        paste0("temp_stat_", position_label, ".csv")
    )
    fwrite(x = temp_stat, file = tmp_file)
    invisible(NULL)
}

#' Consolidate pileup results into a single CSV file
#' @description Internal helper for collectBasesPileup(). Merges all
#'     temporary pileup files into a single sorted CSV, removes duplicates,
#'     and cleans up the temporary files.
#' @param tmpDir character string with the temporary directory.
#' @param bam character string with the full path to the BAM file.
#' @param postfix character string with the output postfix, or NULL.
#' @return invisible(NULL).
#' @author Vasily V. Grinev, Dzianis D. Sarnatski, Mikalai M. Yatskou.
#' @keywords internal
#' @noRd
.consolidatePileupResults <- function(tmpDir, bam, postfix) {
    files <- list.files(path = tmpDir, pattern = "^temp_stat.*\\.csv$")
    if (length(x = files) == 0) {
        stop("No temporary pileup files found in: ", tmpDir)
    }
    file_info <- map_dfr(files, function(file) {
        info <- .extractChrPos(file)
        tibble(file = file, chr = info$chr, pos = info$pos)
    })
    files <- file_info %>%
        arrange(chr, pos) %>%
        pull(file)
    if (!is.null(x = postfix)) {
        csv_file <- paste(
            gsub(pattern = ".bam", replacement = "", x = bam),
            postfix, "csv",
            sep = "."
        )
    } else {
        csv_file <- paste(
            gsub(pattern = ".bam", replacement = "", x = bam),
            "PileupResults", "csv",
            sep = "."
        )
    }
    invisible(if (file.exists(csv_file)) {
        file.remove(csv_file)
    })
    for (file in files) {
        full_path <- file.path(tmpDir, file)
        temp_df <- fread(file = full_path)
        if (!file.exists(csv_file)) {
            fwrite(x = temp_df, file = csv_file)
        } else {
            fwrite(x = temp_df, file = csv_file, append = TRUE)
        }
        file.remove(full_path)
    }
    temp_df <- fread(file = csv_file)
    temp_df <- temp_df[!duplicated(x = temp_df), ]
    fwrite(x = temp_df, file = csv_file)
    invisible(NULL)
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
#' @param threads numeric value giving the number of threads used for sorting.
#'     Default value is detectCores(logical = TRUE) - 1.
#' @param overwrite logical value indicating whether an existing output file
#'     with the same name should be overwritten. FALSE by default, in which
#'     case the function stops with an error if the output file already exists.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @return filtered BAM file (optionally, sorted and indexed if index=TRUE)
#'     with postfix "filtered".
#' @author Vasily V. Grinev
#' @examples
#' filterBamFile(
#'     bamDir = "Files_BAM",
#'     bamFile = "SRR21721515.bam",
#'     gr = "genomic_intervals.txt",
#'     flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE),
#'     tlen = 300,
#'     mapq = 20,
#'     index = TRUE,
#'     threads = 10,
#'     overwrite = TRUE,
#'     workDir = "D:/Vasily Grinev"
#' )
#' @export
#' @importFrom Rsamtools scanBamFlag ScanBamParam filterBam
#' @importFrom Rsamtools sortBam indexBam
#' @importFrom GenomicRanges makeGRangesListFromDataFrame
#' @importFrom IRanges IRangesList
#' @importFrom parallel detectCores
#' @importFrom utils read.table
#' @importFrom S4Vectors FilterRules
filterBamFile <- function(bamDir = NULL,
    bamFile,
    gr = NULL,
    flag = scanBamFlag(),
    tlen = NULL,
    mapq = NULL,
    index = FALSE,
    threads = detectCores(logical = TRUE) - 1,
    overwrite = FALSE,
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
    ### Building the filter parameters.
    param <- .buildFilterParam(gr, flag, tlen, mapq, workDir)
    filter <- .buildFilterRules(tlen)
    ### Output file and overwrite check.
    dest <- gsub(
        pattern = "bam",
        replacement = "filtered.bam",
        x = paste(path, bamFile, sep = "/")
    )
    .checkFilteredBamOverwrite(dest, overwrite)
    ### Filtering.
    f <- filterBam(
        file = paste(path, bamFile, sep = "/"),
        index = paste(path, bamFile, sep = "/"),
        filter = filter,
        param = param,
        destination = dest,
        indexDestination = FALSE
    )
    ### Optional sorting and indexing.
    if (isTRUE(x = index)) {
        s <- sortBam(
            file = f,
            destination = gsub(pattern = ".bam", replacement = "", x = f),
            nThreads = threads
        )
        indexBam(files = s)
    }
}

#' Build ScanBamParam for BAM filtering
#' @description Internal helper for filterBamFile(). Builds the ScanBamParam
#'     object with the specified flags, genomic regions, template length and
#'     mapping quality.
#' @param gr character string with the genomic regions file, or NULL.
#' @param flag logical vector with BAM flags.
#' @param tlen integer, maximum template length, or NULL.
#' @param mapq integer, minimum mapping quality, or NULL.
#' @param workDir character string giving the path to working directory.
#' @return ScanBamParam object.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.buildFilterParam <- function(gr, flag, tlen, mapq, workDir) {
    param <- ScanBamParam(flag = flag)
    if (!is.null(x = gr)) {
        GRs <- read.table(
            file = paste(workDir, gr, sep = "/"),
            sep = "\t",
            header = TRUE,
            quote = "\"",
            as.is = TRUE
        )
        GRs$split <- GRs$seqnames
        GRs <- makeGRangesListFromDataFrame(df = GRs, split.field = "split")
        param@which <- IRangesList(GRs)
    }
    if (!is.null(x = tlen)) {
        param@what <- "isize"
    }
    if (!is.null(x = mapq)) {
        param@mapqFilter <- as.integer(x = mapq)
    }
    param
}

#' Build FilterRules for BAM filtering
#' @description Internal helper for filterBamFile(). Builds the FilterRules
#'     object with the template length filter (if requested).
#' @param tlen integer, maximum template length, or NULL.
#' @return FilterRules object.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.buildFilterRules <- function(tlen) {
    if (!is.null(x = tlen)) {
        FilterRules(list(function(x) {
            x$isize <= tlen
        }))
    } else {
        FilterRules()
    }
}

#' Check overwrite flag for filtered BAM file
#' @description Internal helper for filterBamFile(). Removes the output file
#'     if overwrite = TRUE, otherwise stops with an error if the file exists.
#' @param dest character string with the output file name.
#' @param overwrite logical, whether to overwrite.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.checkFilteredBamOverwrite <- function(dest, overwrite) {
    if (file.exists(dest)) {
        if (isTRUE(overwrite)) {
            file.remove(dest)
        } else {
            stop(
                "Output file already exists: ",
                dest,
                ". Set overwrite = TRUE to replace it."
            )
        }
    }
    invisible(NULL)
}


#' Group short reads in BAM file
#' @description This is a wrapper function that uses the low-level function
#'     AddOrReplaceReadGroup() from package GATK to add read group tags to BAM
#'     file of interest. This stage could be unnecessary, if BAM file already
#'     contains read group tags.
#' @param bamDir character string giving the name (or path to and name) of
#'     directory with BAM file(s). NULL by default, which means the current
#'     working directory will be used instead of the specified BAM directory.
#' @param bamFile character string giving the name of input BAM file.
#' @param postfix character string that will be added to the output file name
#'     as postfix. "RG" by default.
#' @param SE if TRUE (by default), standard error determined by function
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
#' @param workDir character string giving the path to and name of working
#'     directory. NULL by default that means the current working directory.
#'     Both bamDir (when provided) and bamFile must be somewhere inside
#'     directory, chosen as work directory.
#' @return a new BAM file with additional information about read groups and
#'     (optionally) auxillary file with standard error log.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @examples
#' groupReadsGATK(
#'     bamDir = "Files_BAM",
#'     bamFile = "SRR21721515.filtered.bam",
#'     postfix = "RG",
#'     SE = TRUE,
#'     RGID = "Bam.Recal.RG",
#'     RGLB = "default.LIBRARY.ID",
#'     RGPL = "default.ILLUMINA",
#'     RGPU = "default.PLATFORM",
#'     RGSM = "default.SAMPLE",
#'     gatk_path = "gatk",
#'     workDir = "D:/Vasily Grinev"
#' )
#' @export
groupReadsGATK <- function(bamDir = NULL,
    bamFile,
    postfix = "RG",
    SE = TRUE,
    RGID = "Bam.Recal.RG",
    RGLB = "default.LIBRARY.ID",
    RGPL = "default.ILLUMINA",
    RGPU = "default.PLATFORM",
    RGSM = "default.SAMPLE",
    gatk_path = "gatk",
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
    ### Naming of output file and settings for saving of standard error.
    output <- .buildReadGroupsOutputName(bamFile, postfix)
    se <- .buildReadGroupsErrorName(bamFile, SE)
    ### Running of command.
    system2(
        command = gatk_path,
        args = c(
            "AddOrReplaceReadGroups",
            paste("-I", paste(path, bamFile, sep = "/"), sep = " "),
            paste("-O", paste(path, output, sep = "/"), sep = " "),
            paste("--RGID", RGID, "--RGLB", RGLB, "--RGPL", RGPL,
                "--RGPU", RGPU, "--RGSM", RGSM,
                sep = " "
            )
        ),
        if (se == "") {
            stderr <- se
        } else {
            stderr <- paste(path, se, sep = "/")
        }
    )
}

#' Build output file name for groupReadsGATK()
#' @description Internal helper for groupReadsGATK(). Builds the output BAM
#'     file name from the input BAM file name and postfix.
#' @param bamFile character string giving the name of input BAM file.
#' @param postfix character string that will be added as postfix.
#' @return character string with the output file name.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.buildReadGroupsOutputName <- function(bamFile, postfix) {
    paste(
        gsub(pattern = ".bam", replacement = "", x = bamFile),
        postfix, "bam",
        sep = "."
    )
}

#' Build standard error file name for groupReadsGATK()
#' @description Internal helper for groupReadsGATK(). Builds the standard
#'     error log file name, or returns "" when SE = FALSE.
#' @param bamFile character string giving the name of input BAM file.
#' @param SE logical, save standard error?
#' @return character string with the standard error file name, or "".
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.buildReadGroupsErrorName <- function(bamFile, SE) {
    if (SE == TRUE) {
        paste(
            gsub(pattern = ".bam", replacement = "", x = bamFile),
            "STDERR.RG.TXT",
            sep = "."
        )
    } else {
        ""
    }
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
#'     Both fastaDir (when provided) and faFile must be somewhere inside
#'     directory, chosen as working directory.
#' @return index files in working directory.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @examples
#' indexGenomeGATK(
#'     fastaDir = "Files_FASTA",
#'     faFile = "Homo_sapiens.GRCh38.dna.primary_assembly.fa",
#'     fai = TRUE,
#'     dict = TRUE,
#'     gatk_path = "gatk",
#'     workDir = "D:/Vasily Grinev"
#' )
#' @export
#' @importFrom Rsamtools indexFa
indexGenomeGATK <- function(fastaDir = NULL,
    faFile,
    fai = TRUE,
    dict = TRUE,
    gatk_path = "gatk",
    workDir = NULL) {
    ### Setting the working directory.
    if (is.null(x = workDir)) {
        workDir <- getwd()
    }
    ### Full path to the FASTA/FA file.
    if (is.null(x = fastaDir)) {
        path <- workDir
    } else {
        path <- paste(workDir, fastaDir, sep = "/")
    }
    ### Building the .fai index.
    if (fai == TRUE) {
        indexFa(file = paste(path, faFile, sep = "/"))
    }
    ### Building the .dict index.
    if (dict == TRUE) {
        system2(
            command = gatk_path,
            args = c(
                "CreateSequenceDictionary",
                "-R",
                paste(path, faFile, sep = "/")
            )
        )
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
#'     default, which means that the name of input VCF file will be used.
#'     Important note: to be used correctly by the GATK, the VCF file and its
#'     index must have the same name.
#' @param gatk_path character string giving the path to internal or external
#'     (system) GATK executive script. This path can be obtained by running
#'     Find.Gatk() or settled manually. Default value is "gatk", i.e. function
#'     will uses GATK systemly from $PATH.
#' @param workDir character string giving the path to and name of working
#'     directory. NULL by default that means the current working directory.
#'     Both vcfDir (when provided) and vcfFile must be somewhere inside
#'     directory, chosen as working directory.
#' @return .tbi file in working directory.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @examples
#' indexVcfGATK(
#'     vcfDir = "Files_VCF",
#'     vcfFile = "00-All.vcf.gz",
#'     output_name = NULL,
#'     gatk_path = "gatk",
#'     workDir = "D:/Vasily Grinev"
#' )
#' @export
indexVcfGATK <- function(vcfDir = NULL,
    vcfFile = NULL,
    output_name = NULL,
    gatk_path = "gatk",
    workDir = NULL) {
    ### Setting the working directory.
    if (is.null(x = workDir)) {
        workDir <- getwd()
    }
    ### Full path to the VCF file.
    if (is.null(x = vcfDir)) {
        path <- workDir
    } else {
        path <- paste(workDir, vcfDir, sep = "/")
    }
    ### Running of command.
    if (is.null(x = output_name)) {
        system2(
            command = gatk_path,
            args = c(
                "IndexFeatureFile -I",
                paste(path, vcfFile, sep = "/")
            )
        )
    } else {
        system2(
            command = gatk_path,
            args = c(
                "IndexFeatureFile -I",
                paste(path, vcfFile, sep = "/"),
                "-O", paste(path, output_name, sep = "/")
            )
        )
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
#' @param SE if TRUE (by default), standard error determined by function
#'     MarkDuplicates() will be writen in file STDERR.MARKDUP.TXT.
#' @param gatk_path character string giving the path to internal or external
#'     (system) GATK executive script. This path can be obtained by running
#'     Find.Gatk() or settled manually. Default value is "gatk", i.e. function
#'     will uses GATK systemly from $PATH.
#' @param workDir character string giving the path to and name of working
#'     directory. NULL by default that means the current working directory.
#'     Both bamDir (when provided) and bamFile must be somewhere inside
#'     directory, chosen as working directory.
#' @return a new BAM file with processed duplicated reads, a file with metrics
#'     of input BAM file and (optionally) a file with standard error log.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @examples
#' markDuplicatesGATK(
#'     bamDir = "Files_BAM",
#'     bamFile = "SRR21721515.RG.bam",
#'     postfix = "markedDuplicates",
#'     removeDupls = "no",
#'     SE = TRUE,
#'     gatk_path = "gatk",
#'     workDir = "D:/Vasily Grinev"
#' )
#' @export
markDuplicatesGATK <- function(bamDir = NULL,
    bamFile,
    postfix = "markedDuplicates",
    removeDupls = "no",
    SE = TRUE,
    gatk_path = "gatk",
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
    ### Settings of run.
    args_run <- .buildMarkDuplicatesArgs(path, bamFile, postfix, removeDupls)
    ### Settings for saving of standard error.
    se <- ""
    if (SE == TRUE) {
        se <- paste(
            gsub(pattern = ".bam", replacement = "", x = bamFile),
            "STDERR.MARKDUP.TXT",
            sep = "."
        )
    }
    ### Running of command.
    system2(
        command = gatk_path,
        args = c("MarkDuplicates", args_run),
        if (se == "") {
            stderr <- se
        } else {
            stderr <- paste(path, se, sep = "/")
        }
    )
}

#' Build MarkDuplicates arguments
#' @description Internal helper for markDuplicatesGATK(). Builds the GATK
#'     MarkDuplicates arguments: input, output, metrics and the
#'     duplicate-removal flags.
#' @param path character string with the BAM directory.
#' @param bamFile character string with the BAM file name.
#' @param postfix character string with the output postfix.
#' @param removeDupls character string ("no", "yes" or "seq").
#' @return character vector with GATK MarkDuplicates arguments.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.buildMarkDuplicatesArgs <- function(path, bamFile, postfix, removeDupls) {
    input <- paste("-I", paste(path, bamFile, sep = "/"), sep = " ")
    output <- paste("-O",
        paste(path,
            paste(
                gsub(pattern = ".bam", replacement = "", x = bamFile),
                postfix, "bam",
                sep = "."
            ),
            sep = "/"
        ),
        sep = " "
    )
    metrics <- paste("-M",
        paste(path,
            paste(
                gsub(
                    pattern = ".bam",
                    replacement = "", x = bamFile
                ),
                postfix, "metrics.txt",
                sep = "."
            ),
            sep = "/"
        ),
        sep = " "
    )
    remove.d <- .resolveRemoveDuplicates(removeDupls)
    c(input, output, metrics, remove.d)
}

#' Resolve removeDuplicates flag
#' @description Internal helper for markDuplicatesGATK(). Converts the
#'     user-provided removeDupls value to the corresponding GATK arguments.
#' @param removeDupls character string ("no", "yes" or "seq").
#' @return character string with GATK arguments.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.resolveRemoveDuplicates <- function(removeDupls) {
    if (removeDupls == "no") {
        paste("--REMOVE_DUPLICATES", "FALSE")
    } else if (removeDupls == "yes") {
        paste("--REMOVE_DUPLICATES", "TRUE")
    } else if (removeDupls == "seq") {
        paste(
            "--REMOVE_DUPLICATES", "FALSE",
            "--REMOVE_SEQUENCING_DUPLICATES", "TRUE"
        )
    } else {
        paste("--REMOVE_DUPLICATES", "FALSE")
    }
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
#'     genome in format FASTA or FA.
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
#'     All vcfDir, fastaDir, bamDir, as vcfFile, bamFile, fastaFile must be
#'     somewhere inside directory, chosen as work directory.
#' @return if createBQSR=TRUE, a new BQSR table in working directory;
#'     if recalibrate=TRUE, a new BAM file with recalibrated sequencing quality
#'     scores of bases in working directory; if SE=TRUE, a file with standart
#'     error log in working directory.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @examples
#' recalibrateBamGATK(
#'     bamDir = "Files_BAM",
#'     bamFile = "SRR21721515.markedDuplicates.bam",
#'     vcfDir = "Files_VCF",
#'     vcfFile = "00-All.vcf.gz",
#'     fastaDir = "Files_FASTA",
#'     faFile = "Homo_sapiens.GRCh38.dna.primary_assembly.fa",
#'     createBQSR = TRUE,
#'     recalibrate = TRUE,
#'     nameBQSR = "SRR21721515.bqsr.table",
#'     postfix = "Recalibrated",
#'     SE = TRUE,
#'     gatk_path = "gatk",
#'     workDir = "D:/Vasily Grinev"
#' )
#' @export
recalibrateBamGATK <- function(
        bamDir = NULL,
        bamFile,
        vcfDir = NULL,
        vcfFile,
        fastaDir = NULL,
        faFile,
        createBQSR = TRUE,
        recalibrate = TRUE,
        nameBQSR,
        postfix = "Recalibrated",
        SE = TRUE,
        gatk_path = "gatk",
        workDir = NULL) {
    ### Setting the working directory.
    if (is.null(x = workDir)) {
        workDir <- getwd()
    }
    ### Full paths to the input files.
    path <- if (is.null(x = bamDir)) {
        workDir
    } else {
        paste(workDir, bamDir, sep = "/")
    }
    path2 <- if (is.null(x = vcfDir)) {
        workDir
    } else {
        paste(workDir, vcfDir, sep = "/")
    }
    path3 <- if (is.null(x = fastaDir)) {
        workDir
    } else {
        paste(workDir, fastaDir, sep = "/")
    }
    ### Generation of BQSR table.
    if (isTRUE(x = createBQSR)) {
        .runBaseRecalibrator(
            path, bamFile, path2, vcfFile, path3, faFile,
            gatk_path
        )
    }
    ### Recalibration process.
    if (isTRUE(x = recalibrate)) {
        bqsr_name <- .resolveBqsrName(bamFile, createBQSR, nameBQSR)
        se <- if (SE) {
            paste(
                gsub(pattern = ".bam", replacement = "", x = bamFile),
                "STDERR.RECALIBRATION.TXT",
                sep = "."
            )
        } else {
            ""
        }
        .runApplyBQSR(path, bamFile, bqsr_name, postfix, se, gatk_path)
    }
}

#' Run BaseRecalibrator
#' @description Internal helper for recalibrateBamGATK(). Runs the GATK
#'     BaseRecalibrator tool to generate the BQSR table.
#' @param path character string with the BAM directory.
#' @param bamFile character string with the BAM file name.
#' @param path2 character string with the VCF directory.
#' @param vcfFile character string with the VCF file name.
#' @param path3 character string with the FASTA directory.
#' @param faFile character string with the FASTA file name.
#' @param gatk_path character string with the GATK executable.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.runBaseRecalibrator <- function(path, bamFile, path2, vcfFile, path3,
    faFile, gatk_path) {
    bqsr_name <- paste(
        gsub(pattern = ".bam", replacement = "", x = bamFile),
        "bqsr.table",
        sep = "."
    )
    system2(
        command = gatk_path,
        args = c(
            "BaseRecalibrator",
            paste("-I", paste(path, bamFile, sep = "/"), sep = " "),
            paste("-O", paste(path, bqsr_name, sep = "/"), sep = " "),
            paste("--known-sites", paste(path2, vcfFile, sep = "/"),
                sep = " "
            ),
            paste("-R", paste(path3, faFile, sep = "/"), sep = " ")
        )
    )
    invisible(NULL)
}

#' Resolve BQSR table name
#' @description Internal helper for recalibrateBamGATK(). Returns the name of
#'     the BQSR table to use for recalibration.
#' @param bamFile character string with the BAM file name.
#' @param createBQSR logical, whether the BQSR table was created by the
#'     function itself.
#' @param nameBQSR character string with the name of a pre-computed BQSR
#'     table, or NULL.
#' @return character string with the BQSR table name.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.resolveBqsrName <- function(bamFile, createBQSR, nameBQSR) {
    if (isTRUE(x = createBQSR)) {
        paste(
            gsub(pattern = ".bam", replacement = "", x = bamFile),
            "bqsr.table",
            sep = "."
        )
    } else {
        nameBQSR
    }
}

#' Run ApplyBQSR
#' @description Internal helper for recalibrateBamGATK(). Runs the GATK
#'     ApplyBQSR tool to produce the recalibrated BAM file.
#' @param path character string with the BAM directory.
#' @param bamFile character string with the BAM file name.
#' @param bqsr_name character string with the BQSR table name.
#' @param postfix character string with the output postfix.
#' @param se character string with the standard error file, or "".
#' @param gatk_path character string with the GATK executable.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.runApplyBQSR <- function(path, bamFile, bqsr_name, postfix, se, gatk_path) {
    output <- paste(path,
        paste(
            gsub(pattern = ".bam", replacement = "", x = bamFile),
            postfix, "bam",
            sep = "."
        ),
        sep = "/"
    )
    system2(
        command = gatk_path,
        args = c(
            "ApplyBQSR",
            paste("-bqsr", paste(path, bqsr_name, sep = "/"), sep = " "),
            paste("-I", paste(path, bamFile, sep = "/"), sep = " "),
            paste("-O", output, sep = " ")
        ),
        if (se == "") {
            stderr <- se
        } else {
            stderr <- paste(path, se, sep = "/")
        }
    )
    invisible(NULL)
}
