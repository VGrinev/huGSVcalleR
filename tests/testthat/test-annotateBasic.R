test_that("annotateBasic writes basic annotation TXT", {
    vcf <- system.file("extdata", "test_annotation.vcf", package = "huGSVcalleR")
    genes <- system.file("extdata", "ref_genes.txt", package = "huGSVcalleR")
    skip_if(vcf == "" || genes == "", "test data not available")
    td <- tempdir()
    file.copy(vcf, file.path(td, "test.vcf"), overwrite = TRUE)
    file.copy(genes, file.path(td, "ref_genes.txt"), overwrite = TRUE)
    annotateBasic(
        vcfDir = NULL, vcfFile = "test.vcf",
        ref_genes = "ref_genes.txt",
        output = "test", workDir = td
    )
    expect_true(file.exists(
        file.path(td, "test.basic annotation.txt")
    ))
    anno <- read.delim(file.path(td, "test.basic annotation.txt"))
    expect_true(all(c(
        "snv_id", "seqnames", "start", "end", "strand",
        "ref", "alt", "pos_depth", "score", "gene_name"
    ) %in%
        colnames(anno)))
    expect_equal(anno$gene_name, c("GENE1", "GENE2"))
    expect_equal(anno$start, c(100, 300))
})
