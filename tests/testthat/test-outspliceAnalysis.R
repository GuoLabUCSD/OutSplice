test_that("outsplice function runs", {
  dir <- paste0(tempdir(), '/')
  results_matrix <- outspliceAnalysis('../../inst/extdata/HNSC_junctions.txt.gz', '../../inst/extdata/HNSC_genes_normalized.txt.gz', '../../inst/extdata/Total_Rawcounts.txt', '../../inst/extdata/HNSC_pheno_table.txt', 'outsplice_unit_test', dir, filterSex=TRUE, annotation='org.Hs.eg.db', TxDb='TxDb.Hsapiens.UCSC.hg19.knownGene', offsets_value=0.00001, correction_setting='fdr', p_value=0.05)
  expect_no_error(results_matrix)
})
