% File src/library/base/man/load.Rd
\name{PlotJunctionData}
\alias{PlotJunctionData}
\title{Create Plots of Junction Expression}
\description{
  Create Bar and Waterfall plots of raw junction expression, overall gene expression, and junction expression normalized by gene expression.
}
\usage{
PlotJunctionData(data_file, NUMBER=1, junctions=NULL, tail=NULL, p_value = 0.05, GENE=FALSE, SYMBOL=NULL, makepdf=FALSE, pdffile = NULL, tumcol='red', normcol='blue')
}
\arguments{
  \item{data_file}{An R data file containing OutSplice Output.}
  \item{NUMBER}{An integer indicating the number of junctions to plot. This can
    be top number of junctions (over or under expressed), or can be specific junctions in a list. Default is 1.}
  \item{junctions}{A vector of user specified junctions that should be plotted. Default is NULL.}
  \item{tail}{A character string indicating either "RIGHT" to plot the top over
    expressed junctions, or "LEFT" to plot the top under expressed junctions. Default is NULL.}
  \item{p_value}{Set the alpha value for the significance threshold.}
  \item{GENE}{A boolean indicating whether to plot junctions by a specific gene.
    TRUE means you will pick all the junctions mapping to a certain gene. FALSE
    means you do not pick based on the gene. Default is NULL.}
  \item{SYMBOL}{The HGNC gene symbol of the gene to be graphed. Default is NULL.}
  \item{makepdf}{A boolean specifying whether or not to save plots to a PDF. Default is FALSE.}
  \item{pdffile}{A character string giving the file path to the desired pdf. Default is NULL.}
  \item{tumcol}{A character string defining the color of the tumor samples in the plots. Default is red.}
  \item{normcol}{A character string defining the color of the normal samples in the plots. Default is blue.}
}
\examples{
data_file <- system.file("extdata", "OutSplice_Example_2023-01-06.RDa", package = "OutSplice")
ecm1_junc <- "chr1:150483674-150483933"
pdf <- 'ecm1_expression.pdf'
pdf_output <- paste0(tempdir(), '/', pdf)

PlotJunctionData(data_file, NUMBER=1, junctions=ecm1_junc, tail=NULL, p_value = 0.05, GENE=FALSE, SYMBOL=NULL, makepdf=TRUE, pdffile = pdf_output, tumcol='red', normcol='blue')

print(paste0("Output is located at: ", pdf_output))

}
\value{
NULL. Displays or saves a pdf containing waterfall plots of junction expression.
}
\references{

Cancer Genome Atlas Network. Comprehensive genomic characterization of head and neck squamous cell carcinomas. Nature. 2015 Jan 29;517(7536):576-82. doi: 10.1038/nature14129. PMID: 25631445; PMCID: PMC4311405.

Guo T, Sakai A, Afsari B, Considine M, Danilova L, Favorov AV, Yegnasubramanian S, Kelley DZ, Flam E, Ha PK, Khan Z, Wheelan SJ, Gutkind JS, Fertig EJ, Gaykalova DA, Califano J. A Novel Functional Splice Variant of AKT3 Defined by Analysis of Alternative Splice Expression in HPV-Positive Oropharyngeal Cancers. Cancer Res. 2017 Oct 1;77(19):5248-5258. doi: 10.1158/0008-5472.CAN-16-3106. Epub 2017 Jul 21. PMID: 28733453; PMCID: PMC6042297.

Liu C, Guo T, Sakai A, Ren S, Fukusumi T, Ando M, Sadat S, Saito Y, Califano JA. A novel splice variant of LOXL2 promotes progression of human papillomavirus-negative head and neck squamous cell carcinoma. Cancer. 2020 Feb 15;126(4):737-748. doi: 10.1002/cncr.32610. Epub 2019 Nov 13. PMID: 31721164.

Liu C, Guo T, Xu G, Sakai A, Ren S, Fukusumi T, Ando M, Sadat S, Saito Y, Khan Z, Fisch KM, Califano J. Characterization of Alternative Splicing Events in HPV-Negative Head and Neck Squamous Cell Carcinoma Identifies an Oncogenic DOCK5 Variant. Clin Cancer Res. 2018 Oct 15;24(20):5123-5132. doi: 10.1158/1078-0432.CCR-18-0752. Epub 2018 Jun 26. PMID: 29945995; PMCID: PMC6440699.

M. F. Ochs, J. E. Farrar, M. Considine, Y. Wei, S. Meshinchi, and R. J.
Arceci. Outlier analysis and top scoring pair for integrated data analysis and
biomarker discovery. IEEE/ACM Trans Comput Biol Bioinform, 11: 520-32, 2014. PMCID: PMC4156935

}
\keyword{file}
