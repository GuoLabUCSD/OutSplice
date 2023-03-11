test_that("outspliceTCGA function works", {
  dir <- paste0(tempdir(), '/')
  results_matrixTCGA <- outspliceTCGA('../../inst/extdata/TCGA_HNSC_junctions.txt', '../../inst/extdata/TCGA_HNSC_genes_normalized.txt', '../../inst/extdata/Total_Rawcounts.txt', 'outspliceTCGA_unit_test', dir, filterSex=TRUE, annotation='org.Hs.eg.db', TxDb='TxDb.Hsapiens.UCSC.hg19.knownGene', offsets_value=0.00001, correction_setting='fdr', p_value=0.05)
  expect_no_error(results_matrixTCGA)
})
