## analyze junctions from RNA based sequencing data

#' Analyzes differential splicing events between tumor and normal samples.
#'
#' @title Alternative Splicing Analyses
#' @param junction A character string giving the path to a tab separated text file with raw junction counts. The first column includes all junctions to be looked at by OutSplice (Ex: chr1: 1-100). Each proceeding column is a sample with the raw count information for each corresponding junction. The header row contains the name of the junction column, and the names of the samples.
#' @param gene_expr A character string giving the path to a tab separated file with normalized gene expression data. The first column are the entrez ids for each gene, and each proceeding column should be a sample with the normalized expression values for each gene. The file header row contains the name of the entrez id column, and the names of the samples.
#' @param rawcounts A character string giving the path to a tab separated text file with the reads per million counts for each sample. This file can either include a row with the total counts per sample, or multiple rows with raw counts per gene per sample that will be summed automatically by OutSplice. The first column includes user defined row names and the subsequent columns are the sample's rawcount information. The header row contains the name of the row names column and the names of the samples.
#' @param sample_labels A character string giving the path to a tab separated text file with a matrix of tumor and normal labels (T/F) for each sample. The first column should include the names of the samples, and the second column should include "T" for tumors and "F" for normals. The header row contains user defined column names.
#' @param output_file_prefix A character string giving the name of the prefix the user would like to use for the output data file.
#' @param dir A character string giving the path to the directory the user would like to save output to.
#' @param filterSex A boolean representing whether or not to include junctions found on the sex chromosomes.
#' @param annotation A connection or a character string giving the name of the Bioconductor library the user would like to use containing the genome wide annotation.
#' @param TxDb A character string giving the name of the Bioconductor library the user would like to use that will expose the annotation database as a TxDb object.
#' @param offsets_value The minimum expression value needed to call an event an outlier after normalizing event expression with gene expression.
#' @param correction_setting Option to designate how to correct significance. The available options are: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", and "none".
#' @param p_value Set the alpha value for the significance threshold.
#' @return Matrix of junction events containing the number of under/over-expressed outliers in the tumor group (outRank1/outRank2), the Fisher p-value for under/over-expressed events (FisherP1/FisherP2), and vectors that can be used to order and rank the under/over expressed events (var1/var2): Ex) FisherAnalyses[FisherAnalyses[,'var2'],]
#' @examples
#' junction <- system.file("extdata", "HNSC_junctions.txt.gz", package = "OutSplice")
#' gene_expr <- system.file("extdata", "HNSC_genes_normalized.txt.gz", package = "OutSplice")
#' rawcounts <- system.file("extdata", "Total_Rawcounts.txt", package = "OutSplice")
#' sample_labels <- system.file("extdata", "HNSC_pheno_table.txt", package = "OutSplice")
#' output_file_prefix <- "OutSplice_Example"
#' TxDb_hg19 <- "TxDb.Hsapiens.UCSC.hg19.knownGene"
#' dir <- paste0(tempdir(), "/")
#' outspliceAnalysis(junction, gene_expr, rawcounts, sample_labels, output_file_prefix, dir, filterSex = TRUE, annotation = "org.Hs.eg.db", TxDb = TxDb_hg19, offsets_value = 0.00001, correction_setting = "fdr", p_value = 0.05)
#' print(paste0("Output is located at: ", dir))
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
outspliceAnalysis <- function(junction, gene_expr, rawcounts, sample_labels, output_file_prefix, dir, filterSex = TRUE, annotation = "org.Hs.eg.db", TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene", offsets_value = 0.00001, correction_setting = "fdr", p_value = 0.05) {
  checkOSArgs(junction, gene_expr, rawcounts, output_file_prefix, dir)
  stopifnot("Sample Labels File does not exist. Check path to file." = file.exists(sample_labels))

  date <- Sys.Date()
  ## this is non-log transformed data and includes pheno (RAW, in RPM)
  print("Loading data")
  all.junc <- read.table(file = junction, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  all.gene_expr <- read.table(file = gene_expr, header = TRUE, row.names = 1, sep = "\t")
  rawcounts <- read.table(file = rawcounts, sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  samps.labels <- read.table(file = sample_labels, sep = "\t", row.names = 1, header = TRUE, stringsAsFactors = FALSE)

  colnames(rawcounts) <- gsub("\\.", "-", colnames(rawcounts))

  all.samples <- colnames(all.junc)

  all.samples <- gsub("\\.", "-", all.samples)
  junction.names <- all.junc[, 1]
  # remove duplicates
  all.junc <- all.junc[!duplicated(junction.names), ]
  junction.names <- all.junc[, 1]

  rownames(all.junc) <- junction.names
  colnames(all.junc) <- all.samples

  # remove first column
  all.junc <- all.junc[, -1]


  expression.samples <- gsub("\\.", "-", colnames(all.gene_expr))
  colnames(all.gene_expr) <- expression.samples

  # Order Columns by Name
  samps.labels_df <- as.data.frame((t(samps.labels)))
  samps.labels_df <- samps.labels_df[, order(colnames(samps.labels_df))]

  all.junc <- all.junc[, order(colnames(all.junc))]

  all.gene_expr <- all.gene_expr[, order(colnames(all.gene_expr))]

  rawcounts <- rawcounts[, order(colnames(rawcounts))]

  # Get Vector of labels
  samps.labels <- (t(samps.labels))
  samps.labels <- samps.labels[, order(colnames(samps.labels))]

  # infer phenotype from sample names
  all.samples <- intersect(colnames(all.junc), colnames(all.gene_expr))
  all.samples <- intersect(all.samples, colnames(rawcounts))
  all.samples <- intersect(all.samples, colnames(samps.labels_df))

  pheno <- samps.labels
  names(pheno) <- all.samples
  # include only tumor or normal
  pheno <- pheno[pheno == TRUE | pheno == FALSE]
  pheno[pheno == TRUE] <- "Tumor"
  pheno[pheno == FALSE] <- "Normal"
  print(table(pheno))
  #
  if (sum(pheno == "Normal") < 10) {
    stop("Too few normal samples")
  }
  # subset only primary tumor 01 or normal 11
  all.samples <- names(pheno)
  all.junc <- all.junc[, all.samples]
  all.gene_expr <- all.gene_expr[, all.samples]
  rawcounts <- rawcounts[, all.samples]
  pheno <- pheno[all.samples]

  ## change from char to numeric
  n <- vapply(all.gene_expr, as.numeric, numeric(nrow(all.gene_expr)))
  rownames(n) <- rownames(all.gene_expr)
  all.gene_expr <- n
  n <- vapply(all.junc, as.numeric, numeric(nrow(all.junc)))
  rownames(n) <- rownames(all.junc)
  all.junc <- n
  remove(n)

  ## normalize Junctions to get RPM
  junc.RPM <- getRPM(rawcounts, all.junc, all.samples)

  ###### NEW FILTER: See if tumors have any outliers ###########################
  ## PHENO should have 'Normal' or 'Tumor' calls where Tumor ==1, Normal ==0, and names of each sample associated
  PHENO <- pheno == "Tumor"
  PHENO <- as.numeric(PHENO)
  names(PHENO) <- names(pheno)

  print("run the ogsa function for pre filtering")
  ## get function

  test2 <- dotheogsa(Sample.data = junc.RPM, PHENO = PHENO, offsets = 0.1, dir = dir)
  has.outliers <- test2[, "outRankTumor1"] > 1 | test2[, "outRankTumor2"] > 1
  junc.RPM <- junc.RPM[has.outliers, ]
  #
  #  ### Filter all genes on the X and Y chromosomes
  #
  if (filterSex) {
    junc.RPM <- junc.RPM[grep("chr[XY]", row.names(junc.RPM), value = TRUE, invert = TRUE), ]
  }

  # ##############################################################
  print("get the genomic information for all the junctions")

  geneAnnot <- getGenomicInfo(junc.RPM, annotation, TxDb)

  print("deletions")
  print("insertions")
  print("skipping")

  geneAnnot <- geneAnnot[apply(cbind(
    geneAnnot$deletions, geneAnnot$insertions,
    geneAnnot$skipping
  ), 1, any), ]
  junc.RPM <- junc.RPM[names(geneAnnot), ]

  ##################################################################
  print("remove all that map to 'NA' no gene name, and assign gene expression from gene_expr")
  junc.RPM <- junc.RPM[!is.na(geneAnnot$SYMBOL), ]
  geneAnnot <- geneAnnot[row.names(junc.RPM), ]

  ## remove row names for unknown genes containing "?" unknown genes

  all.gene_expr <- all.gene_expr[!grepl("\\?", row.names(all.gene_expr)), ]

  ## get gene IDs
  gene_exprEntrezID <- row.names(all.gene_expr)

  print("align with gene_expr data")

  junctionGenegene_expr <- getExpressions(geneAnnot, all.gene_expr, gene_exprEntrezID)

  # remove junctions without gene expression
  junc.RPM <- junc.RPM[row.names(junctionGenegene_expr), ]


  print("subset removing any genes without normalization")
  junc.RPM2 <- junc.RPM[!(apply(junctionGenegene_expr, 1, sum)) == 0, ]
  dim(junc.RPM2)
  # ## filters junctions
  junc.RPM.original <- junc.RPM
  junc.RPM <- junc.RPM2

  print("Perform normalization using gene_expr values")

  junc.RPM.norm <- normalizeJunctions(junc.RPM2, junctionGenegene_expr)

  ### perform outlier analysis with OGSA
  ############### Use OGSA for outlier ranking #################

  print("run the ogsa function")

  FisherAnalyses <- dotheogsa(Sample.data = junc.RPM.norm, PHENO = PHENO, offsets = offsets_value, Fisher = TRUE, correction = correction_setting, dir = dir)
  ## use default offset=0.001 for normalized data
  # create lists of top genes
  topgenelist10.ogsa <- FisherAnalyses[, "var1"]

  topgenelist90.ogsa <- FisherAnalyses[, "var2"]


  ## Use Fisher to subset candidates that have Fisher test p value <0.05, or cutoff can be adjusted
  ## Underexpression in tumors
  toplist10 <- subset_underEx(FisherAnalyses, topgenelist10.ogsa, p_value)

  ## Overexpression in tumors
  toplist90 <- subset_overEx(FisherAnalyses, topgenelist90.ogsa, p_value)

  junctions <- c(names(toplist90), names(toplist10))
  pvalues <- c(toplist90, toplist10)

  ### get outier calls##
  junc.Outliers <- dotheogsa(Sample.data = junc.RPM.norm, PHENO = PHENO, offsets = offsets_value, Fisher = TRUE, correction = correction_setting, outliers = TRUE, dir = dir)

  ## median of normal expression within the normalized data
  NORM.gene_expr.norm <- apply(junc.RPM.norm[, pheno == "Normal"], 1, median)

  ## aggregate the data
  gene_expr <- junctionGenegene_expr
  junc.RPM.norm <- junc.RPM.norm
  geneAnnotations <- geneAnnot
  geneAnnot <- geneAnnot[junctions]
  ASE.type <- cbind(geneAnnot$skipping, geneAnnot$insertions, geneAnnot$deletions)
  colnames(ASE.type) <- c("skipping", "insertions", "deletions")
  row.names(ASE.type) <- junctions

  ## Calculate Splicing Burden
  splice_burden <- calcBurden(junc.Outliers, FisherAnalyses, p_value)

  ## save output file
  save(junc.RPM, gene_expr, junc.RPM.norm, pvalues, pheno, FisherAnalyses, geneAnnotations, ASE.type, NORM.gene_expr.norm, junc.Outliers, splice_burden, file = paste0(dir, output_file_prefix, "_", date, ".RDa"))

  # Write Files
  write.table(ASE.type, file = paste0(dir, output_file_prefix, "_", date, "_", "event_types.txt"), sep = "\t", quote = FALSE, col.names = NA)
  write.table(FisherAnalyses, file = paste0(dir, output_file_prefix, "_", date, "_", "FisherAnalyses.txt"), sep = "\t", quote = FALSE, col.names = NA)
  write.table(as.data.frame(junc.Outliers$TumorOverExpression), file = paste0(dir, output_file_prefix, "_", date, "_", "TumorOverExpression.txt"), sep = "\t", quote = FALSE, col.names = NA)
  write.table(as.data.frame(junc.Outliers$TumorUnderExpression), file = paste0(dir, output_file_prefix, "_", date, "_", "TumorUnderExpression.txt"), sep = "\t", quote = FALSE, col.names = NA)
  write.table(splice_burden, file = paste0(dir, output_file_prefix, "_", date, "_", "splice_burden.txt"), sep = "\t", quote = FALSE, col.names = NA)

  annotations_df <- annoGR2DF(geneAnnotations)
  write.table(annotations_df, file = paste0(dir, output_file_prefix, "_", date, "_", "gene_annotations.txt"), sep = "\t", quote = FALSE, col.names = NA)

  return(FisherAnalyses)
}
