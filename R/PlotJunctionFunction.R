# Plot Junction Expression Data

#' Create Bar and Waterfall plots of raw junction expression, overall gene expression, and junction expression normalized by gene expression for splicing events found by OutSplice.
#'
#' @title Create Plots of Junction Expression
#' @param data_file An R data file containing OutSplice Output.
#' @param NUMBER An integer indicating the number of junctions to plot. This can be top number of junctions (over or under expressed), or can be specific junctions in a list. Default is 1.
#' @param junctions A vector of user specified junctions that should be plotted. Default is NULL.
#' @param tail A character string indicating either "RIGHT" to plot the top over expressed junctions, or "LEFT" to plot the top under expressed junctions. Default is NULL.
#' @param p_value Set the alpha value for the significance threshold.
#' @param GENE A boolean indicating whether to plot junctions by a specific gene. TRUE means you will pick all the junctions mapping to a certain gene. FALSE means you do not pick based on the gene. Default is NULL.
#' @param SYMBOL The HGNC gene symbol of the gene to be graphed. Default is NULL.
#' @param makepdf A boolean specifying whether or not to save plots to a PDF. Default is FALSE.
#' @param pdffile A character string giving the file path to the desired pdf. Default is NULL.
#' @param tumcol A character string defining the color of the tumor samples in the plots. Default is red.
#' @param normcol A character string defining the color of the normal samples in the plots. Default is blue.
#' @return NULL. Displays or saves a pdf containing waterfall plots of junction expression.
#' @examples
#' data_file <- system.file("extdata", "OutSplice_Example_2023-01-06.RDa", package = "OutSplice")
#' ecm1_junc <- "chr1:150483674-150483933"
#' pdf <- "ecm1_expression.pdf"
#' pdf_output <- paste0(tempdir(), "/", pdf)
#' plotJunctionData(data_file, NUMBER = 1, junctions = ecm1_junc, tail = NULL, p_value = 0.05, GENE = FALSE, SYMBOL = NULL, makepdf = TRUE, pdffile = pdf_output, tumcol = "red", normcol = "blue")
#' message("Output is located at: ", pdf_output)
#' @references
#' Cancer Genome Atlas Network. Comprehensive genomic characterization of head and neck squamous cell carcinomas. Nature. 2015 Jan 29;517(7536):576-82. doi: 10.1038/nature14129. PMID: 25631445; PMCID: PMC4311405.
#'
#' Guo T, Sakai A, Afsari B, Considine M, Danilova L, Favorov AV, Yegnasubramanian S, Kelley DZ, Flam E, Ha PK, Khan Z, Wheelan SJ, Gutkind JS, Fertig EJ, Gaykalova DA, Califano J. A Novel Functional Splice Variant of AKT3 Defined by Analysis of Alternative Splice Expression in HPV-Positive Oropharyngeal Cancers. Cancer Res. 2017 Oct 1;77(19):5248-5258. doi: 10.1158/0008-5472.CAN-16-3106. Epub 2017 Jul 21. PMID: 28733453; PMCID: PMC6042297.
#'
#' Liu C, Guo T, Sakai A, Ren S, Fukusumi T, Ando M, Sadat S, Saito Y, Califano JA. A novel splice variant of LOXL2 promotes progression of human papillomavirus-negative head and neck squamous cell carcinoma. Cancer. 2020 Feb 15;126(4):737-748. doi: 10.1002/cncr.32610. Epub 2019 Nov 13. PMID: 31721164.
#'
#' Liu C, Guo T, Xu G, Sakai A, Ren S, Fukusumi T, Ando M, Sadat S, Saito Y, Khan Z, Fisch KM, Califano J. Characterization of Alternative Splicing Events in HPV-Negative Head and Neck Squamous Cell Carcinoma Identifies an Oncogenic DOCK5 Variant. Clin Cancer Res. 2018 Oct 15;24(20):5123-5132. doi: 10.1158/1078-0432.CCR-18-0752. Epub 2018 Jun 26. PMID: 29945995; PMCID: PMC6440699.
#'
#' M. F. Ochs, J. E. Farrar, M. Considine, Y. Wei, S. Meshinchi, and R. J. Arceci. Outlier analysis and top scoring pair for integrated data analysis and biomarker discovery. IEEE/ACM Trans Comput Biol Bioinform, 11: 520-32, 2014. PMCID: PMC4156935
#' @export
plotJunctionData <- function(data_file, NUMBER = 1, junctions = NULL, tail = NULL, p_value = 0.05, GENE = FALSE, SYMBOL = NULL, makepdf = FALSE, pdffile = NULL, tumcol = "red", normcol = "blue") {
    stopifnot("R Data File does not exist. Check path to file." = file.exists(data_file))
    if (!NUMBER > 0) {
        stop("NUMBER must be an integer greater than 0.")
    }
    if (!is.null(SYMBOL)) {
        if (GENE == FALSE) stop("Need to set the GENE argument to TRUE.")
    }
    suppressPackageStartupMessages(load(data_file))
    ## if you want to make a pdf, this will be specified.  Stop/error if not specified
    ## the number of junctions to be ploted should be equal or greater than the number of significant junctions that exist
    if (makepdf == TRUE) {
        if (is.null(pdffile)) {
            stop("Need to specify PDF file path name. pdffile=...")
        }
        if (isTRUE(tail) && tail == "RIGHT") {
            if (sum(FisherAnalyses[, "FisherP2"] < p_value) < NUMBER) {
                stop("Number of significant junctions is less than the number of junctions specified.")
            }
        } else if (isTRUE(tail) && tail == "LEFT") {
            if (sum(FisherAnalyses[, "FisherP1"] < p_value) < NUMBER) {
                stop("Number of significant junctions is less than the number of junctions specified.")
            }
        }
        pdf(pdffile)
    }
    ## Define which junctions are going to be graphed
    ## if you do not define which specific junctions you want to graph
    if (is.null(junctions)) {
        if (GENE == FALSE) {
            if (is.null(tail)) {
                junctions <- row.names(junc.RPM)[seq_len(NUMBER)]
            } else if (tail == "RIGHT") {
                toplist90 <- subset_overEx(FisherAnalyses, p_value)
                junctions <- row.names(toplist90)[seq_len(NUMBER)]
                message("Top ", NUMBER, " overexpressed junctions graphed.")
            } else if (tail == "LEFT") {
                toplist10 <- subset_underEx(FisherAnalyses, p_value)
                junctions <- row.names(toplist10)[seq_len(NUMBER)]
                message("Top ", NUMBER, " underexpressed junctions graphed.")
            }
        } else if (GENE == TRUE) {
            if (is.null(SYMBOL)) stop("Need to specify gene of interest by gene symbol. SYMBOL='...'")
            junctions <- names(geneAnnotations)[grep(SYMBOL, geneAnnotations$SYMBOL)]
            message(length(junctions), " junction(s) for gene ", SYMBOL, " found.")
        }
    }

    if (length(junctions) == 0) stop("No junctions found.")
    count <- 0
    par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
    for (j in junctions)
    {
        ## first barplot is for raw junction data
        message(j)
        phenotypes <- as.data.frame(pheno)
        samples <- order(phenotypes$pheno, junc.RPM.norm[j, ], junc.RPM.norm[j, ])
        barplot(log2(junc.RPM[j, samples] + 1), cex.names = 0.3, las = 2, col = ifelse(phenotypes[samples, "pheno"] == "Tumor", tumcol, normcol))
        title(sprintf(
            "Junction Expression (log)\nskip? %s del? %s? ins? %s",
            geneAnnotations[j, ]$skipping, geneAnnotations[j, ]$deletions,
            geneAnnotations[j, ]$insertions
        ), cex.main = 0.7)

        ## this plot is for normalized junction expression - with gene_expr.
        barplot(log2(junc.RPM.norm[j, samples] + 1) - log2(NORM.gene_expr.norm[j] + 1),
            las = 2, cex.names = 0.3,
            col = ifelse(phenotypes[samples, "pheno"] == "Tumor", tumcol, normcol)
        )
        title("Junction expression \nNormalized by gene_expr (log)", cex.main = 0.7)

        ## this barplot is for gene_expr of gene expression
        barplot(log2(gene_expr[j, samples] + 1),
            las = 2, cex.names = 0.3,
            col = ifelse(phenotypes[samples, "pheno"] == "Tumor", tumcol, normcol)
        )
        title(sprintf("%s \ngene_expr gene expression", geneAnnotations[j, ]$SYMBOL), cex.main = 0.7)

        plot.new()

        count <- count + 1
    }



    if (makepdf == TRUE) dev.off()
}
