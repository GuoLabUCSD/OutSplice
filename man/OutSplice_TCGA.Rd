% File src/library/base/man/load.Rd
\name{OutSplice_TCGA}
\alias{OutSplice_TCGA}
\title{Alternative Splicing Analysis for TCGA Data}
\description{
  Analyze differential splicing events between tumor and normal samples for TCGA formatted datasets.
  Examples of TCGA file formats can be viewed here (https://gdac.broadinstitute.org/)
}
\usage{
OutSplice_TCGA(junction, gene_expr, rawcounts, output_file_prefix, dir, filterSex=TRUE, genome = 'Homo.sapiens', annotation = 'org.Hs.eg.db', TxDb = 'TxDb.Hsapiens.UCSC.hg19.knownGene', offsets_value = 0.00001, correction_setting='fdr', p_value=0.05)
}
\arguments{
  \item{junction}{A character string giving the path to a tab separated text file with raw juntion counts.}
  \item{gene_expr}{A character string giving the path to a tab separated file with normalized gene
    expression data.}
  \item{rawcounts}{A character string giving the path to a tab separated text file with the reads
    per million counts for each sample.}
  \item{output_file_prefix}{A character string giving the name of the prefix the
    user would like to use for the output data file.}
  \item{dir}{A character string giving the path to the directory the user would
    like to save output to.}
  \item{filterSex}{A boolean representing whether or not to include junctions
    found on the sex chromosomes.}
  \item{genome}{A character string giving the name of the Bioconductor library
    the user would like to use containing the genome object.}
  \item{annotation}{A connection or a character string giving the name of the
    Bioconductor library the user would like to use containing the genome wide annotation.}
  \item{TxDb}{A character string giving the name of the Bioconductor library the
    user would like to use that will expose the annotation database as a TxDb object.}
  \item{offsets_value}{The minimum expression value needed to call an event an outlier after
    normalizing event expression with gene expression.}
  \item{correction_setting}{Option to designate how to correct significance.}
  \item{p_value}{Set the alpha value for the significance threshold.}
}
\examples{
## peform differential splicing analysis

#Get Path to Example Files and establish variables

junction <- system.file("extdata", "TCGA_HNSC_junctions.txt", package="OutSplice")
gene_expr <- system.file("extdata", "TCGA_HNSC_genes_normalized.txt", package="OutSplice")
rawcounts <- system.file("extdata", "Total_Rawcounts.txt", package="OutSplice")
output_file_prefix = "TCGA_OutSplice_Example"
dir = paste0(tempdir(), '/')

#Run the OutSplice FUnction for TCGA data

OutSplice_TCGA(junction, gene_expr, rawcounts, output_file_prefix, dir, filterSex=TRUE, genome = 'Homo.sapiens', annotation = 'org.Hs.eg.db', TxDb = 'TxDb.Hsapiens.UCSC.hg19.knownGene', offsets_value=0.00001, correction_setting='fdr', p_value=0.05)

print(paste0("Output is located at: ", dir))

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
