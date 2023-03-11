#Function to Calculate Splicing Burden Based on the Number of Fisher Significant Events in Each Sample

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
#' pdf <- 'ecm1_expression.pdf'
#' pdf_output <- paste0(tempdir(), '/', pdf)
#' plotJunctionData(data_file, NUMBER=1, junctions=ecm1_junc, tail=NULL, p_value = 0.05, GENE=FALSE, SYMBOL=NULL, makepdf=TRUE, pdffile = pdf_output, tumcol='red', normcol='blue')
#' print(paste0("Output is located at: ", pdf_output))
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
calcBurden <- function(junc.Outliers, FisherAnalyses, p_value) {

  #Calculate Splicing Burden for Significant Events Over-Expressed in Tumors
  Over_Expressed_Junctions_Outlier_Calls <- as.data.frame(junc.Outliers$TumorOverExpression)
  row.names(Over_Expressed_Junctions_Outlier_Calls) <- sub('\\.', ':', rownames(Over_Expressed_Junctions_Outlier_Calls))
  row.names(Over_Expressed_Junctions_Outlier_Calls) <- sub('\\.', '-', rownames(Over_Expressed_Junctions_Outlier_Calls))

  OE_pvalue <- subset(FisherAnalyses, select = FisherP2)
  Over_Expressed_Junctions_Outlier_Calls <- merge(OE_pvalue, Over_Expressed_Junctions_Outlier_Calls, by='row.names', all.x=TRUE)
  rownames(Over_Expressed_Junctions_Outlier_Calls) <- Over_Expressed_Junctions_Outlier_Calls[, 1]
  Over_Expressed_Junctions_Outlier_Calls <- subset(Over_Expressed_Junctions_Outlier_Calls, select = -Row.names)
  Over_Expressed_Junctions_Outlier_Calls <- na.omit(Over_Expressed_Junctions_Outlier_Calls)
  Over_Expressed_Junctions_Outlier_Calls <- subset(Over_Expressed_Junctions_Outlier_Calls, Over_Expressed_Junctions_Outlier_Calls$FisherP2 < p_value)
  Over_Expressed_Junctions_Outlier_Calls <- subset(Over_Expressed_Junctions_Outlier_Calls, select = -FisherP2)

  OutlierNumberOver = colSums(Over_Expressed_Junctions_Outlier_Calls)
  all_oe_results_df = data.frame(OutlierNumberOver)

  #Calculate Splicing Burden for Significant Events Under-Expressed in Tumors
  Under_Expressed_Junctions_Outlier_Calls <- as.data.frame(junc.Outliers$TumorUnderExpression)
  row.names(Under_Expressed_Junctions_Outlier_Calls) <- sub('\\.', ':', rownames(Under_Expressed_Junctions_Outlier_Calls))
  row.names(Under_Expressed_Junctions_Outlier_Calls) <- sub('\\.', '-', rownames(Under_Expressed_Junctions_Outlier_Calls))

  UE_pvalue <- subset(FisherAnalyses, select = FisherP1)
  Under_Expressed_Junctions_Outlier_Calls <- merge(UE_pvalue, Under_Expressed_Junctions_Outlier_Calls, by='row.names', all.x=TRUE)
  rownames(Under_Expressed_Junctions_Outlier_Calls) <- Under_Expressed_Junctions_Outlier_Calls[, 1]
  Under_Expressed_Junctions_Outlier_Calls <- subset(Under_Expressed_Junctions_Outlier_Calls, select = -Row.names)
  Under_Expressed_Junctions_Outlier_Calls <- na.omit(Under_Expressed_Junctions_Outlier_Calls)
  Under_Expressed_Junctions_Outlier_Calls <- subset(Under_Expressed_Junctions_Outlier_Calls, Under_Expressed_Junctions_Outlier_Calls$FisherP1 < p_value)
  Under_Expressed_Junctions_Outlier_Calls <- subset(Under_Expressed_Junctions_Outlier_Calls, select = -FisherP1)

  OutlierNumberUnder = colSums(Under_Expressed_Junctions_Outlier_Calls)
  all_ue_results_df = data.frame(OutlierNumberUnder)

  #Total Splicing Burden as the sum of Over + Under Expressed Events
  total_results = transform(merge(all_oe_results_df, all_ue_results_df, by = 0, all = TRUE), TotalOutliers = OutlierNumberOver + OutlierNumberUnder)
  rownames(total_results) <- total_results[, 1]
  convert_total <- subset(total_results, select = c(-OutlierNumberOver, -OutlierNumberUnder, -Row.names))
  total_burden <- as.data.frame(convert_total)

  return(splice_burden <- cbind(all_oe_results_df, all_ue_results_df, total_burden))
}
