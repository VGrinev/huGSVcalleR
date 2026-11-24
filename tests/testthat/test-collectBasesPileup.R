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
