## analyze junctions from RNA based sequencing data

#' Analyzes differential splicing events between tumor and normal samples.
#'
#' @title Alternative Splicing Analyses
#' @param junction A character string giving the path to a tab separated text file with raw junction counts. One column should include all of the junctions to be looked at by OutSplice (Ex: chr1: 1-100). Each proceeding column is a sample with the raw count information for each corresponding junction. The header row contains the name of the junction column, and the names of the samples.
#' @param gene_expr A character string giving the path to a tab separated file with normalized gene expression data. One column should include all of the entrez ids for each gene, and each proceeding column should be a sample with the normalized expression values for each gene. The file header row contains the name of the entrez id column, and the names of the samples.
#' @param rawcounts A character string giving the path to a tab separated text file with the reads per million counts for each sample. This file can either include a row with the total counts per sample, or multiple rows with raw counts per gene per sample that will be summed automatically by OutSplice. One of the columns should include the user defined row names and the subsequent columns are the sample's rawcount information. The header row contains the name of the row names column and the names of the samples.
#' @param sample_labels A character string giving the path to a tab separated text file with a matrix of tumor and normal labels (T/F) for each sample. One of the columns should include the names of the samples, and the other column should include "T" for tumors and "F" for normals. The header row contains user defined column names.
#' @param saveOutput A boolean representing whether or not to save the results to an R data file and tab separated files. Default is FALSE. Optional.
#' @param output_file_prefix A character string giving the name of the prefix the user would like to use for the output data file. Default is NULL. Optional.
#' @param dir A character string giving the path to the directory the user would like to save output to. Default is NULL. Optional.
#' @param filterSex A boolean representing whether or not to include junctions found on the sex chromosomes. Default is TRUE. Optional.
#' @param annotation A connection or a character string giving the name of the Bioconductor library the user would like to use containing the genome wide annotation. Default is "org.Hs.eg.db". Optional.
#' @param TxDb A character string giving the name of the Bioconductor library the user would like to use that will expose the annotation database as a TxDb object. Default is "TxDb.Hsapiens.UCSC.hg38.knownGene". Optional.
#' @param offsets_value The minimum expression value needed to call an event an outlier after normalizing event expression with gene expression. Default is 0.00001. Optional.
#' @param correction_setting Option to designate how to correct significance. The available options are: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", and "none". Default is "fdr". Optional.
#' @param p_value Set the alpha value for the significance threshold. Default is 0.05. Optional.
#' @param use_junc_col An integer indicating which column in the junction matrix contains the junction regions in your matrices. Default is 1. Optional
#' @param use_gene_col An integer indicating which column in the gene_expr matrix contains the entrez ids of your genes. Default is 1. Optional
#' @param use_rc_col An integer indicating which column in the rawcounts matrix contains the row names. Default is 1. Optional
#' @param use_labels_col An integer indicating which column in the sample_labels matrix contains the sample names. Default is 1. Optional
#' @return A list containing the below data.
#' \itemize{
#'     \item FisherAnalyses: Data Frame of junction events containing the number of under/over-expressed outliers in the tumor group (Num_UE_Outliers/Num_OE_Outliers), the Fisher p-value for under/over-expressed events (FisherP1/FisherP2), and a ranking of the under/over expressed events (UE_Rank/OE_Rank)
#'     \item ASE.type: significant junction events labeled by type (skipping, insertion, or deletion)
#'     \item geneAnnotations: object containing gene names corresponding to each junction region
#'     \item junc.Outliers: list containing the logical matrices TumorOverExpression and TumorUnderExpression. "True" indicates an over-expressed event in TumorOverExpression, or an under-expressed event in TumorUnderExpression.
#'     \item junc.RPM: junction counts in reads per million following a division of the junction counts input by the total rawcounts for each sample
#'     \item junc.RPM.norm: junction counts normalized by each event's total gene expression value
#'     \item gene_expr: gene expression values for each junction event
#'     \item splice_burden: matrix containing the number of Fisher-P significant over-expressed, under-expressed, and total number of outliers per sample
#'     \item NORM.gene_expr.norm: Median of junction data normalized by gene expression for normal samples only (Used for Junction Plotting Only)
#'     \item pheno: Phenotypes of Samples (Tumor or Normal)
#'     \item pvalues: Junction Fisher P-values
#' }
#' @examples
#' junction <- system.file("extdata", "HNSC_junctions.txt.gz", package = "OutSplice")
#' gene_expr <- system.file("extdata", "HNSC_genes_normalized.txt.gz", package = "OutSplice")
#' rawcounts <- system.file("extdata", "Total_Rawcounts.txt", package = "OutSplice")
#' sample_labels <- system.file("extdata", "HNSC_pheno_table.txt", package = "OutSplice")
#' output_file_prefix <- "OutSplice_Example"
#' TxDb_hg19 <- "TxDb.Hsapiens.UCSC.hg19.knownGene"
#' dir <- paste0(tempdir(), "/")
#' results <- outspliceAnalysis(junction, gene_expr, rawcounts, sample_labels, saveOutput = TRUE, output_file_prefix, dir, filterSex = TRUE, annotation = "org.Hs.eg.db", TxDb = TxDb_hg19, offsets_value = 0.00001, correction_setting = "fdr", p_value = 0.05)
#' message("Output is located at: ", dir)
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
outspliceAnalysis <- function(junction, gene_expr, rawcounts, sample_labels,
                              saveOutput = FALSE, output_file_prefix = NULL,
                              dir = NULL, filterSex = TRUE, annotation = "org.Hs.eg.db",
                              TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
                              offsets_value = 0.00001, correction_setting = "fdr", p_value = 0.05,
                              use_junc_col = 1, use_gene_col = 1, use_rc_col = 1, use_labels_col = 1) {
    checkOSArgs(junction, gene_expr, rawcounts)
    stopifnot("Sample Labels File does not exist. Check path to file." = file.exists(sample_labels))

    if (saveOutput) {
        checkDirArgs(output_file_prefix, dir)
    }

    date <- Sys.Date()

    ## this is non-log transformed data and includes pheno (RAW, in RPM)
    message("Loading data")
    all.junc <- read.table(file = junction, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    all.gene_expr <- read.table(file = gene_expr, header = TRUE, row.names = use_gene_col, sep = "\t")
    rawcounts <- read.table(
        file = rawcounts, sep = "\t", header = TRUE,
        row.names = use_rc_col, stringsAsFactors = FALSE
    )
    samps.labels <- read.table(
        file = sample_labels, sep = "\t", row.names = use_labels_col,
        header = TRUE, stringsAsFactors = FALSE
    )

    colnames(rawcounts) <- gsub("\\.", "-", colnames(rawcounts))

    all.samples <- colnames(all.junc)

    all.samples <- gsub("\\.", "-", all.samples)
    junction.names <- all.junc[, use_junc_col]
    # remove duplicates
    all.junc <- all.junc[!duplicated(junction.names), ]
    junction.names <- all.junc[, use_junc_col]

    rownames(all.junc) <- junction.names
    colnames(all.junc) <- all.samples

    # remove first column
    all.junc <- all.junc[, -use_junc_col]

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
    message(paste0(capture.output(table(pheno)), collapse = "\n"))
    #
    if (sum(pheno == "Normal") < 10) {
        stop("Too few normal samples")
    }
    ## Run Data Processing Function
    data <- processMatrices1(
        pheno, all.junc, all.samples, all.gene_expr, rawcounts,
        filterSex, annotation, TxDb
    )

    ## get gene IDs
    gene_exprEntrezID <- row.names(data$all.gene_expr)

    message("align with gene_expr data")

    junctionGenegene_expr <- getExpressions(data$geneAnnot, data$all.gene_expr, gene_exprEntrezID)

    ## Run Second Data Processing Function

    final_results <- processMatrices2(
        data$junc.RPM, junctionGenegene_expr, data$PHENO,
        offsets_value, correction_setting, p_value,
        data$pheno, data$geneAnnot, saveOutput, output_file_prefix, dir, date
    )

    return(final_results)
}
