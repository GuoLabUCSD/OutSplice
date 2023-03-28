## analyze TCGA junctions from TCGA based sequencing data
## TCGA Firehose pipeline

#' Analyze differential splicing events between tumor and normal samples for TCGA formatted datasets. Examples of TCGA file formats can be viewed here (https://gdac.broadinstitute.org/)
#'
#' @title Alternative Splicing Analysis for TCGA Data
#' @param junction A character string giving the path to a tab separated text file with raw juntion counts.
#' @param gene_expr A character string giving the path to a tab separated file with normalized gene expression data.
#' @param rawcounts A character string giving the path to a tab separated text file with the reads per million counts for each sample.
#' @param output_file_prefix A character string giving the name of the prefix the user would like to use for the output data file.
#' @param dir A character string giving the path to the directory the user would like to save output to.
#' @param filterSex A boolean representing whether or not to include junctions found on the sex chromosomes.
#' @param annotation A connection or a character string giving the name of the Bioconductor library the user would like to use containing the genome wide annotation.
#' @param TxDb A character string giving the name of the Bioconductor library the user would like to use that will expose the annotation database as a TxDb object.
#' @param offsets_value The minimum expression value needed to call an event an outlier after normalizing event expression with gene expression.
#' @param correction_setting Option to designate how to correct significance. The available options are: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", and "none".
#' @param p_value Set the alpha value for the significance threshold.
#' @return Data Frame of junction events containing the number of under/over-expressed outliers in the tumor group (Num_UE_Outliers/Num_OE_Outliers), the Fisher p-value for under/over-expressed events (FisherP1/FisherP2), and a ranking of the under/over expressed events (UE_Rank/OE_Rank). This function will also output tab sepaerated text files and an R data file with the following:
#' \itemize{
#'     \item FisherAnalyses: Data Frame of the above return value
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
#' junction <- system.file("extdata", "TCGA_HNSC_junctions.txt.gz", package = "OutSplice")
#' gene_expr <- system.file("extdata", "TCGA_HNSC_genes_normalized.txt.gz", package = "OutSplice")
#' rawcounts <- system.file("extdata", "Total_Rawcounts.txt", package = "OutSplice")
#' output_file_prefix <- "TCGA_OutSplice_Example"
#' dir <- paste0(tempdir(), "/")
#' outspliceTCGA(junction, gene_expr, rawcounts, output_file_prefix, dir, filterSex = TRUE, annotation = "org.Hs.eg.db", TxDb = "TxDb.Hsapiens.UCSC.hg19.knownGene", offsets_value = 0.00001, correction_setting = "fdr", p_value = 0.05)
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
outspliceTCGA <- function(junction, gene_expr, rawcounts, output_file_prefix, dir, filterSex = TRUE, annotation = "org.Hs.eg.db", TxDb = "TxDb.Hsapiens.UCSC.hg19.knownGene", offsets_value = 0.00001, correction_setting = "fdr", p_value = 0.05) {
    checkOSArgs(junction, gene_expr, rawcounts, output_file_prefix, dir)

    date <- Sys.Date()

    ## this is non-log transformed data and includes pheno (RAW, in RPM)
    message("Loading data")
    all.junc <- read.table(file = junction, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    all.gene_expr <- read.table(file = gene_expr, header = TRUE, row.names = 1, sep = "\t")
    rawcounts <- read.table(file = rawcounts, sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

    colnames(rawcounts) <- substr(colnames(rawcounts), 1, 15)
    colnames(rawcounts) <- gsub("\\.", "-", colnames(rawcounts))

    # remove header; rename row names and junctions
    all.junc <- all.junc[-1, ] # remove header
    all.gene_expr <- all.gene_expr[-1, ] # remove header

    all.samples <- colnames(all.junc)
    all.samples <- substr(all.samples, 1, 15)
    all.samples <- gsub("\\.", "-", all.samples)
    junction.names <- all.junc[, 1]
    # remove duplicates
    all.junc <- all.junc[!duplicated(junction.names), ]
    junction.names <- all.junc[, 1]
    j <- paste0(vapply(strsplit(junction.names, split = ":"), function(x) {
        x[[1]]
    }, character(1)), ":", vapply(strsplit(junction.names, split = ":"), function(x) {
        x[[2]]
    }, character(1)), "-", vapply(strsplit(junction.names, split = ":"), function(x) {
        x[[4]]
    }, character(1)))

    rownames(all.junc) <- j
    colnames(all.junc) <- all.samples
    # remove first column
    all.junc <- all.junc[, -1]

    expression.samples <- substr(colnames(all.gene_expr), 1, 15)
    expression.samples <- gsub("\\.", "-", expression.samples)
    colnames(all.gene_expr) <- expression.samples

    # infer phenotype from sample names
    all.samples <- intersect(colnames(all.junc), colnames(all.gene_expr))
    all.samples <- intersect(all.samples, colnames(rawcounts))
    pheno <- substr(all.samples, 14, 15)
    names(pheno) <- all.samples
    # include only primary tumor 01 or normal 11
    pheno <- pheno[pheno == "01" | pheno == "11"]
    pheno[pheno == "01"] <- "Tumor"
    pheno[pheno == "11"] <- "Normal"
    message(paste0(capture.output(table(pheno)), collapse = "\n"))

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

    message("run the ogsa function for pre filtering")
    ## get function
    test2 <- dotheogsa(Sample.data = junc.RPM, PHENO = PHENO, offsets = 0.1, dir = dir)
    has.outliers <- test2[, "Num_UE_Outliers"] > 1 | test2[, "Num_OE_Outliers"] > 1
    junc.RPM <- junc.RPM[has.outliers, ]

    ### Filter all genes on the X and Y chromosomes

    if (filterSex) {
        junc.RPM <- junc.RPM[grep("chr[XY]", row.names(junc.RPM), value = TRUE, invert = TRUE), ]
    }
    ##############################################################
    message("get the genomic information for all the junctions")

    geneAnnot <- getGenomicInfo(junc.RPM, annotation, TxDb)

    message("deletions")
    message("insertions")
    message("skipping")

    geneAnnot <- geneAnnot[apply(cbind(
        geneAnnot$deletions, geneAnnot$insertions,
        geneAnnot$skipping
    ), 1, any), ]
    junc.RPM <- junc.RPM[names(geneAnnot), ]

    ##################################################################
    message("remove all that map to 'NA' no gene name, and assign gene expression from gene_expr")
    junc.RPM <- junc.RPM[!is.na(geneAnnot$SYMBOL), ]
    geneAnnot <- geneAnnot[row.names(junc.RPM), ]

    ## remove row names for unknown genes containing "?" unknown genes

    all.gene_expr <- all.gene_expr[!grepl("\\?", row.names(all.gene_expr)), ]

    ## get gene names
    gene_exprgenenames <- vapply(strsplit(row.names(all.gene_expr), split = "\\|"), function(x) {
        x[1]
    }, character(1))
    gene_exprEntrezID <- vapply(strsplit(row.names(all.gene_expr), split = "\\|"), function(x) {
        x[2]
    }, character(1))

    message("align with gene_expr data")

    junctionGenegene_expr <- getExpressions(geneAnnot, all.gene_expr, gene_exprEntrezID)

    ## remove junctions without gene expression
    junc.RPM <- junc.RPM[row.names(junctionGenegene_expr), ]


    message("subset removing any genes without normalization")
    junc.RPM2 <- junc.RPM[!(apply(junctionGenegene_expr, 1, sum)) == 0, ]
    dim(junc.RPM2)
    # ## filters junctions
    junc.RPM.original <- junc.RPM
    junc.RPM <- junc.RPM2

    message("Perform normalization using gene_expr values")

    junc.RPM.norm <- normalizeJunctions(junc.RPM2, junctionGenegene_expr)

    ### perform outlier analysis with OGSA
    ############### Use OGSA for outlier ranking #################

    message("run the ogsa function")

    FisherAnalyses <- dotheogsa(Sample.data = junc.RPM.norm, PHENO = PHENO, offsets = offsets_value, Fisher = TRUE, correction = "fdr", dir = dir)
    ## use default offset=0.001 for normalized data
    # create lists of top genes

    ## Use Fisher to subset candidates that have Fisher test p value <0.05, or cutoff can be adjusted
    ## Underexpression in tumors
    toplist10 <- subset_underEx(FisherAnalyses, p_value)

    ## Overexpression in tumors
    toplist90 <- subset_overEx(FisherAnalyses, p_value)

    junctions <- c(rownames(toplist90), rownames(toplist10))
    pvalues <- c(toplist90$FisherP2, toplist10$FisherP1)

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
