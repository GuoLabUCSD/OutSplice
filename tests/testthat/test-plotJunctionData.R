test_that("Junction Plotting Works", {
  pdf_output <- paste0(tempdir(), '/', 'ecm1_expression.pdf')
  junc_plots <- plotJunctionData('../../inst/extdata/OutSplice_Example_2023-01-06.RDa', NUMBER=1, junctions="chr1:150483674-150483933", tail=NULL, p_value = 0.05, GENE=FALSE, SYMBOL=NULL, makepdf=TRUE, pdffile = pdf_output, tumcol='red', normcol='blue')
  expect_no_error(junc_plots)
})
