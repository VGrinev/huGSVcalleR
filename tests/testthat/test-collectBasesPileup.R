test_that("collectBasesPileup writes A C G T counts", {
    skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")
    bam <- system.file("extdata", "test.bam", package = "huGSVcalleR")
    bai <- system.file("extdata", "test.bam.bai", package = "huGSVcalleR")
    skip_if(bam == "" || bai == "", "test BAM not available")
    td <- tempfile("pileup-")
    dir.create(td)
    file.copy(bam, file.path(td, "test.bam"), overwrite = TRUE)
    file.copy(bai, file.path(td, "test.bam.bai"), overwrite = TRUE)
    collectBasesPileup(
        bamDir = NULL,
        bamFile = "test.bam",
        gr = "SeqInfo",
        baseq = 0,
        mapq = 0,
        coverage = NULL,
        threads = 1,
        postfix = "example",
        workDir = td
    )
    csv <- file.path(td, "test.example.csv")
    expect_true(file.exists(csv))
    tab <- data.table::fread(csv)
    expect_true(all(c("position", "reference", "A", "C", "G", "T") %in%
        colnames(tab)))
    expect_false("T_count" %in% colnames(tab))
    expect_gt(nrow(tab), 0)
    expect_true(all(nchar(as.character(tab$reference)) == 1))
    expect_true(all(tab$A >= 0 & tab$C >= 0 & tab$G >= 0 & tab$T >= 0))
    expect_gt(sum(tab$A + tab$C + tab$G + tab$T), 0)
})

test_that("parallel pileup workers use the master summary", {
    pile <- data.frame(
        seqnames = "chr1",
        pos = 10L,
        position = "chr1:10",
        nucleotide = c("A", "C", "G", "T"),
        count = c(3L, 1L, 2L, 8L)
    )
    bp <- BiocParallel::SnowParam(
        workers = 2L,
        exportglobals = FALSE,
        stop.on.error = TRUE
    )
    res <- BiocParallel::bplapply(1:2, function(i, summarize, pile) {
        tab <- summarize(pile)
        base_cols <- c("A", "C", "G", "T")
        depth_ok <- rowSums(tab[, base_cols, with = FALSE]) >= 10
        filtered <- tab[which(depth_ok)]
        c(
            paste(colnames(filtered), collapse = ","),
            as.character(nrow(filtered))
        )
    }, summarize = huGSVcalleR:::.summarizePileupCounts, pile = pile,
        BPPARAM = bp)
    expect_true(all(vapply(res, function(x) {
        grepl("T", x[1], fixed = TRUE) && !grepl("T_count", x[1], fixed = TRUE)
    }, logical(1))))
    expect_true(all(vapply(res, function(x) as.integer(x[2]) > 0L, logical(1))))
})
