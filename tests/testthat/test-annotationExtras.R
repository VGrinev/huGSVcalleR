test_that("scoreSSs scores a canonical donor site", {
    skip_if(Sys.which("perl") == "", "perl is not available")
    sites <- Biostrings::DNAStringSet(c(donor = "CAGGTAAGT"))
    scores <- local({
        old <- setwd(system.file("extdata", package = "huGSVcalleR"))
        on.exit(setwd(old))
        scoreSSs(x = sites, type = 5)
    })
    expect_equal(scores$ss_seq, "CAGGTAAGT")
    expect_equal(scores$max.ent_score, 10.86)
})

test_that("extractSeqSSs returns the hg38 bases", {
    skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")
    seqs <- extractSeqSSs(
        x = GenomicRanges::GRanges("chr1:100000-100008:+")
    )
    expect_s4_class(seqs, "DNAStringSet")
    expect_equal(unname(as.character(seqs)), "CACTAAGCA")
})

test_that("introduceSNVtoSSs substitutes the alternate base", {
    skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")
    genome <- getExportedValue(
        "BSgenome.Hsapiens.UCSC.hg38",
        "BSgenome.Hsapiens.UCSC.hg38"
    )
    pos <- 100000L
    ref <- as.character(BSgenome::getSeq(
        genome,
        GenomicRanges::GRanges("chr1", IRanges::IRanges(pos, pos))
    ))
    alt <- setdiff(c("A", "C", "G", "T"), ref)[1]
    td <- tempfile("splice-")
    dir.create(td)
    writeLines(
        c(
            "##fileformat=VCFv4.2",
            "##contig=<ID=chr1,length=248956422>",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            paste("chr1", pos, ".", ref, alt, "10", "PASS", ".", sep = "\t")
        ),
        file.path(td, "splice.vcf")
    )
    res <- introduceSNVtoSSs(
        x = GenomicRanges::GRanges(
            "chr1:99995-100010:+",
            ss_id = "ss1",
            gene_id = "GENE1",
            transcript_id = "TX1"
        ),
        vcfDir = NULL,
        vcfFile = "splice.vcf",
        workDir = td
    )
    expect_equal(substr(res$annotation$seq_original, 6, 6), ref)
    expect_equal(substr(res$annotation$seq_corrected, 6, 6), alt)
    expect_equal(
        unname(as.character(res$sequences)),
        res$annotation$seq_corrected
    )
})

test_that("annotateDbSNP copies the matching dbSNP allele", {
    ex <- system.file("extdata", package = "huGSVcalleR")
    td <- tempfile("dbsnp-")
    dir.create(td)
    file.copy(
        file.path(ex, "test_annotation.vcf"),
        file.path(td, "test.vcf"),
        overwrite = TRUE
    )
    file.copy(file.path(ex, "ref_genes.txt"), td, overwrite = TRUE)
    annotateBasic(
        vcfDir = NULL,
        vcfFile = "test.vcf",
        ref_genes = "ref_genes.txt",
        output = "example",
        workDir = td
    )
    anno <- read.delim(
        file.path(td, "example.basic annotation.txt"),
        sep = "\t"
    )
    expect_equal(anno$gene_name, c("GENE1", "GENE2"))
    snv <- anno[1, ]
    fields <- c(
        RS = "Integer", dbSNPBuildID = "Integer", SAO = "Integer",
        PM = "Integer", FREQ = "String", COMMON = "Integer",
        CLNHGVS = "String", CLNVI = "String", CLNORIGIN = "String",
        CLNSIG = "String", CLNDISDB = "String", CLNDN = "String",
        CLNREVSTAT = "String", CLNACC = "String"
    )
    header <- c(
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1,length=1000000>",
        paste0(
            "##INFO=<ID=", names(fields),
            ",Number=1,Type=", fields,
            ",Description=\"x\">"
        ),
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    )
    info <- paste(
        "RS=7", "dbSNPBuildID=150", "SAO=0", "PM=0", "FREQ=0.01",
        "COMMON=1", "CLNHGVS=.", "CLNVI=.", "CLNORIGIN=.", "CLNSIG=.",
        "CLNDISDB=.", "CLNDN=.", "CLNREVSTAT=.", "CLNACC=.",
        sep = ";"
    )
    record <- paste(
        snv$seqnames, snv$start, ".", snv$ref, snv$alt,
        ".", "PASS", info,
        sep = "\t"
    )
    writeLines(c(header, record), file.path(td, "dbsnp.vcf"))
    gz <- Rsamtools::bgzip(
        file.path(td, "dbsnp.vcf"),
        dest = file.path(td, "dbsnp.vcf.gz"),
        overwrite = TRUE
    )
    Rsamtools::indexTabix(gz, format = "vcf")
    annotateDbSNP(
        annoDir = NULL,
        annoFile = "example.basic annotation.txt",
        output = "example",
        dbDir = NULL,
        dbFile = "dbsnp.vcf.gz",
        chrToNC = FALSE,
        chrID = "unused.txt",
        workDir = td
    )
    out <- read.delim(file.path(td, "example.dbSNP annotation.txt"))
    expect_equal(out$RS[1], 7)
    expect_true(is.na(out$RS[2]))
})

test_that("developAnnotations writes a GTF when GenomeInfoDbData is present", {
    skip_if_not_installed("GenomeInfoDbData")
    ex <- system.file("extdata", "test_annotation.gtf", package = "huGSVcalleR")
    td <- tempfile("anno-")
    dir.create(td)
    file.copy(ex, file.path(td, "test_annotation.gtf"), overwrite = TRUE)
    write.table(
        data.frame(transcript_id = "TX1", orf_start = 1L, orf_end = 30L),
        file = file.path(td, "orfs.txt"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    developAnnotations(
        gtfDir = NULL,
        gtfFile = "test_annotation.gtf",
        orf = "orfs.txt",
        gtf = "example.anno.gtf",
        src = "example",
        organism = "Homo sapiens",
        sqlite = "example.anno.sqlite",
        workDir = td
    )
    gtf <- rtracklayer::import(file.path(td, "example.anno.gtf"))
    expect_true("gene" %in% gtf$type)
    expect_true(file.exists(file.path(td, "example.anno.sqlite")))
})

test_that("introduceSNVtoAAs stops when the OrgDb package is missing", {
    skip_if(requireNamespace("org.Hs.eg.db", quietly = TRUE))
    td <- tempfile("aa-")
    dir.create(td)
    writeLines(
        c(
            "##fileformat=VCFv4.2",
            "##contig=<ID=chr1,length=1000>",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            paste("chr1", 100, ".", "A", "G", "10", "PASS", ".", sep = "\t")
        ),
        file.path(td, "coding.vcf")
    )
    write.table(
        data.frame(
            transcript_id = "TX",
            gene_id = "GENE",
            protein_id = "PROT",
            pdb_id = NA_character_,
            swissprot = NA_character_,
            trembl = NA_character_
        ),
        file = file.path(td, "proteins.txt"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    expect_error(
        introduceSNVtoAAs(
            vcfDir = NULL,
            vcfFile = "coding.vcf",
            proteins = "proteins.txt",
            workDir = td
        ),
        "org.Hs.eg.db"
    )
})
