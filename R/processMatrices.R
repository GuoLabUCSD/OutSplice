# Shared Sub FUnction for outSpliceAnalysis and outSpliceTCGA for overall data processing

processMatrices1 <- function(pheno, all.junc, all.samples, all.gene_expr, rawcounts, filterSex, annotation, TxDb) {
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
    test2 <- dotheogsa(Sample.data = junc.RPM, PHENO = PHENO, offsets = 0.1)
    #View(test2)
    has.outliers <- test2[, "Num_UE_Outliers"] > 1 | test2[, "Num_OE_Outliers"] > 1
    junc.RPM <- junc.RPM[has.outliers, ]

    ### Filter all genes on the X and Y chromosomes

    if (filterSex) {
        junc.RPM <- junc.RPM[grep("chr[XY]", row.names(junc.RPM), value = TRUE, invert = TRUE), ]
    }
    ##############################################################
    message("get the genomic information for all the junctions")

    geneAnnot <- getGenomicInfo(junc.RPM, annotation, TxDb)
    #View(as.data.frame(geneAnnot))

    #Filter events without event type
    geneAnnot <- geneAnnot[apply(cbind(
        geneAnnot$deletions, geneAnnot$insertions,
        geneAnnot$skipping
    ), 1, any), ]

    #test_this <- as.data.frame(geneAnnot)
    #View(test_this)
    junc.RPM <- junc.RPM[names(geneAnnot), ]

    ##################################################################
    message("remove all that map to 'NA' no gene name, and assign gene expression from gene_expr")
    junc.RPM <- junc.RPM[!is.na(geneAnnot$SYMBOL), ]

    geneAnnot <- geneAnnot[row.names(junc.RPM), ]

    ## remove row names for unknown genes containing "?" unknown genes

    all.gene_expr <- all.gene_expr[!grepl("\\?", row.names(all.gene_expr)), ]

    os_data <- list(
        "junc.RPM" = junc.RPM, "geneAnnot" = geneAnnot,
        "all.gene_expr" = all.gene_expr, "PHENO" = PHENO, "pheno" = pheno
    )

    return(os_data)
}

# Second Shared Sub Function for outSpliceAnalysis and outSpliceTCGA for overall data processing

processMatrices2 <- function(junc.RPM, junctionGenegene_expr, PHENO, offsets_value,
                             correction_setting, p_value, pheno, geneAnnot, saveOutput, output_file_prefix, dir, date) {
    # remove junctions without gene expression
    junc.RPM <- junc.RPM[row.names(junctionGenegene_expr), ]

    message("subset removing any genes without normalization")
    junc.RPM2 <- junc.RPM[!(apply(junctionGenegene_expr, 1, sum)) == 0, ]
    dim(junc.RPM2)
    # ## filters junctions
    junc.RPM.original <- junc.RPM
    junc.RPM <- junc.RPM2
    #View(junc.RPM)

    message("Perform normalization using gene_expr values")

    junc.RPM.norm <- normalizeJunctions(junc.RPM2, junctionGenegene_expr)

    ### perform outlier analysis with OGSA
    ############### Use OGSA for outlier ranking #################

    message("run the ogsa function")

    FisherAnalyses <- dotheogsa(
        Sample.data = junc.RPM.norm, PHENO = PHENO,
        offsets = offsets_value, Fisher = TRUE, correction = correction_setting
    )
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
    junc.Outliers <- dotheogsa(
        Sample.data = junc.RPM.norm, PHENO = PHENO,
        offsets = offsets_value, Fisher = TRUE, correction = correction_setting, outliers = TRUE
    )

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

    ## calculate splicing burden
    splice_burden <- calcBurden(junc.Outliers, FisherAnalyses, p_value)

    ## optionally save output if specified
    if (saveOutput) {
        saveResults(
            junc.RPM, gene_expr, junc.RPM.norm, pvalues, pheno, FisherAnalyses,
            geneAnnotations, ASE.type, NORM.gene_expr.norm, junc.Outliers,
            splice_burden, output_file_prefix, date, dir
        )
    }

    ## create return value
    final_results <- list(
        "junc.RPM" = junc.RPM, "gene_expr" = gene_expr, "junc.RPM.norm" = junc.RPM.norm,
        "pvalues" = pvalues, "pheno" = pheno, "FisherAnalyses" = FisherAnalyses,
        "geneAnnotations" = geneAnnotations, "ASE.type" = ASE.type,
        "NORM.gene_expr.norm" = NORM.gene_expr.norm, "junc.Outliers" = junc.Outliers,
        "splice_burden" = splice_burden
    )
    return(final_results)
}
