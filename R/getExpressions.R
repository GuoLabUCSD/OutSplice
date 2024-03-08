# Collect Gene Expression values from normalized expression matrix

getExpressions <- function(geneAnnot, all.gene_expr, gene_exprEntrezID) {
    ## collect gene_expr values
    # initialize matrix of values for each junction
    junctionGenegene_expr <- matrix(0,
        nrow = length(geneAnnot), ncol = ncol(all.gene_expr),
        dimnames = list(names(geneAnnot), colnames(all.gene_expr))
    )

    geneAnnot$ENTREZID -> genes2Junc_ENTREZ
    names(genes2Junc_ENTREZ) <- names(geneAnnot)

    ## shows how many junctions aligned to a single gene"
    length(grep(genes2Junc_ENTREZ, pattern = ";", invert = TRUE, value = TRUE))

    ### Make sure the genes are ordered by lower to higher ENTREZID
    genes2Junc_ENTREZ <- vapply(genes2Junc_ENTREZ, function(x)
        paste(sort(as.numeric(strsplit(x, ";")[[1]])), collapse = ";"), character(1))

    ## this selects just the first gene that each junction aligns to
    genes2Junc_ENTREZ <- vapply(strsplit(genes2Junc_ENTREZ, ";"), function(x) {
        x[1]
    }, character(1))

    ## fill in the matrix
    no.gene_expr <- vector(length = length(genes2Junc_ENTREZ))
    for (g in seq_along(genes2Junc_ENTREZ)) {
        if (!isTRUE(intersect(gene_exprEntrezID, genes2Junc_ENTREZ[g]) > 0)) {
            ## for gene where there is no gene_expr data, skip it
            no.gene_expr[g] <- TRUE
            next
        }
        junctionGenegene_expr[names(genes2Junc_ENTREZ)[g], ] <-
            all.gene_expr[which(gene_exprEntrezID == genes2Junc_ENTREZ[g]), ]
    }

    ################################################################
    ## Remove genes without gene_expr data, not in true genes
    #print(no.gene_expr)
    junctionGenegene_expr[!no.gene_expr, ] -> junctionGenegene_expr

    return(junctionGenegene_expr)
}
