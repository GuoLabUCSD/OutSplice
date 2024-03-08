# Optional Funciton to save output to tab separated text files

saveResults <- function(junc.RPM, gene_expr, junc.RPM.norm, pvalues, pheno, FisherAnalyses,
                        geneAnnotations, ASE.type, NORM.gene_expr.norm, junc.Outliers,
                        splice_burden, output_file_prefix, date, dir) {
    ## save output file
    save(junc.RPM, gene_expr, junc.RPM.norm, pvalues, pheno, FisherAnalyses, geneAnnotations,
        ASE.type, NORM.gene_expr.norm, junc.Outliers, splice_burden,
        file = paste0(dir, output_file_prefix, "_", date, ".RDa")
    )

    # Write Files
    write.table(ASE.type, file = paste0(
        dir, output_file_prefix, "_",
        date, "_", "event_types.txt"
    ), sep = "\t", quote = FALSE, col.names = NA)
    write.table(FisherAnalyses, file = paste0(
        dir, output_file_prefix, "_",
        date, "_", "FisherAnalyses.txt"
    ), sep = "\t", quote = FALSE, col.names = NA)
    write.table(as.data.frame(junc.Outliers$TumorOverExpression), file = paste0(
        dir, output_file_prefix, "_",
        date, "_", "TumorOverExpression.txt"
    ), sep = "\t", quote = FALSE, col.names = NA)
    write.table(as.data.frame(junc.Outliers$TumorUnderExpression), file = paste0(
        dir, output_file_prefix, "_",
        date, "_", "TumorUnderExpression.txt"
    ), sep = "\t", quote = FALSE, col.names = NA)
    write.table(splice_burden, file = paste0(
        dir, output_file_prefix, "_",
        date, "_", "splice_burden.txt"
    ), sep = "\t", quote = FALSE, col.names = NA)

    write.table(x = data.frame(geneAnnotations), file = paste0(
        dir, output_file_prefix, "_",
        date, "_", 'gene_annotations.txt'
    ), sep = '\t', col.names = TRUE, row.names=FALSE, quote = FALSE)
}
