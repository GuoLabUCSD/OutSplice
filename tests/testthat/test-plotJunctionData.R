test_that("Junction Plotting Works for User Specified Junctions", {
    data_file <- system.file("extdata", "OutSplice_Example_2023-01-06.RDa", package = "OutSplice")
    pdf_output <- paste0(tempdir(), "/", "ecm1_expression.pdf")
    junc_plots <- plotJunctionData(data_file, NUMBER = 1, junctions = "chr1:150483674-150483933", tail = NULL, p_value = 0.05, GENE = FALSE, SYMBOL = NULL, makepdf = TRUE, pdffile = pdf_output, tumcol = "red", normcol = "blue")
    expect_no_error(junc_plots)
})

test_that("Junction Plotting Works for Over-Expressed Tail Values", {
    data_file <- system.file("extdata", "OutSplice_Example_2023-01-06.RDa", package = "OutSplice")
    pdf_output <- paste0(tempdir(), "/", "top2_OE.pdf")
    junc_plots <- plotJunctionData(data_file, NUMBER = 2, tail = "RIGHT", p_value = 0.05, GENE = FALSE, SYMBOL = NULL, makepdf = TRUE, pdffile = pdf_output, tumcol = "red", normcol = "blue")
    expect_no_error(junc_plots)
})

test_that("Junction Plotting Works for Under-Expressed Tail Values", {
    data_file <- system.file("extdata", "OutSplice_Example_2023-01-06.RDa", package = "OutSplice")
    pdf_output <- paste0(tempdir(), "/", "top2_UE.pdf")
    junc_plots <- plotJunctionData(data_file, NUMBER = 2, tail = "LEFT", p_value = 0.05, GENE = FALSE, SYMBOL = NULL, makepdf = TRUE, pdffile = pdf_output, tumcol = "red", normcol = "blue")
    expect_no_error(junc_plots)
})

test_that("Junction Plotting Works for Under-Expressed Tail Values", {
    data_file <- system.file("extdata", "OutSplice_Example_2023-01-06.RDa", package = "OutSplice")
    pdf_output <- paste0(tempdir(), "/", "ecm1_gene_expression.pdf")
    junc_plots <- plotJunctionData(data_file, p_value = 0.05, GENE = TRUE, SYMBOL = "ECM1", makepdf = TRUE, pdffile = pdf_output, tumcol = "red", normcol = "blue")
    expect_no_error(junc_plots)
})
