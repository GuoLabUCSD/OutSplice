test_that("outsplice function runs", {
  dir <- paste0(tempdir(), '/')
  results_matrix <- outsplice('../../inst/extdata/HNSC_junctions.txt', '../../inst/extdata/HNSC_genes_normalized.txt', '../../inst/extdata/Total_Rawcounts.txt', '../../inst/extdata/HNSC_pheno_table.txt', 'outsplice_unit_test', dir, filterSex=TRUE, annotation='org.Hs.eg.db', TxDb='TxDb.Hsapiens.UCSC.hg19.knownGene', offsets_value=0.00001, correction_setting='fdr', p_value=0.05)
  expect_no_error(results_matrix)
})
