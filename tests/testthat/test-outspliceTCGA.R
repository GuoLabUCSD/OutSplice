test_that("outspliceTCGA function works", {
    junction <- system.file("extdata", "TCGA_HNSC_junctions.txt.gz", package = "OutSplice")
    gene_expr <- system.file("extdata", "TCGA_HNSC_genes_normalized.txt.gz", package = "OutSplice")
    rawcounts <- system.file("extdata", "Total_Rawcounts.txt", package = "OutSplice")
    dir <- paste0(tempdir(), "/")
    results_matrixTCGA <- outspliceTCGA(junction, gene_expr, rawcounts, "outspliceTCGA_unit_test", dir, filterSex = TRUE, annotation = "org.Hs.eg.db", TxDb = "TxDb.Hsapiens.UCSC.hg19.knownGene", offsets_value = 0.00001, correction_setting = "fdr", p_value = 0.05)
    expect_no_error(results_matrixTCGA)
    expect_equal(unname(unlist(results_matrixTCGA["chr1:150483674-150483933", "Num_OE_Outliers"])), 28)
})
