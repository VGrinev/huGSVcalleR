test_that("fisher p-values are deterministic and match fisher.test", {
    counts <- data.frame(
        reference = c("A", "A", "A", "A"),
        A = c(10L, 8L, 5L, 0L),
        C = c(0L, 2L, 5L, 10L),
        G = c(0L, 0L, 0L, 0L),
        T = c(0L, 0L, 0L, 0L)
    )
    first <- huGSVcalleR:::snp_fisher_cpp_parallel(counts)
    second <- huGSVcalleR:::snp_fisher_cpp_parallel(counts)
    expect_equal(first, second)
    expect_equal(first[1], 1)
    expected <- vapply(2:4, function(i) {
        ref <- counts$A[i]
        alt <- counts$C[i]
        total <- ref + alt
        tab <- matrix(c(alt, 0, ref, total), nrow = 2)
        fisher.test(tab)$p.value
    }, numeric(1))
    expect_equal(as.numeric(first)[2:4], expected, tolerance = 1e-8)
    expect_true(all(first > 0))
})

test_that("entropy p-values use the natural log and 3 degrees of freedom", {
    counts <- data.frame(
        reference = c("A", "A"),
        A = c(10L, 5L),
        C = c(0L, 5L),
        G = c(0L, 0L),
        T = c(0L, 0L)
    )
    got <- huGSVcalleR:::snp_entropy_cpp_parallel(counts, p_value = TRUE)
    entropy <- -2 * 0.5 * log(0.5)
    expected <- c(1, pchisq(2 * 10 * entropy, df = 3, lower.tail = FALSE))
    expect_equal(as.numeric(got), expected, tolerance = 1e-8)
})
