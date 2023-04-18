# Get Genomic Information for the Junctions and Determine the Event Type

getGenomicInfo <- function(junc.RPM, annotation, TxDb) {
    # create GenomicRanges object for junctions
    chr <- vapply(strsplit(row.names(junc.RPM), split = ":"), function(x) {
        x[[1]]
    }, character(1))
    start <- as.numeric(vapply(
        strsplit(row.names(junc.RPM), split = "[:-]"),
        function(x) {
            x[[2]]
        }, character(1)
    ))
    end <- as.numeric(vapply(
        strsplit(row.names(junc.RPM), split = "[:-]"),
        function(x) {
            x[[3]]
        }, character(1)
    ))


    geneAnnot <- GRanges(
        seqnames = Rle(chr),
        IRanges(start = start, end = end)
    )
    names(geneAnnot) <- row.names(junc.RPM)

    geneAnnotAll <- GRanges(
        seqnames = Rle(chr),
        IRanges(start = start, end = end)
    )
    names(geneAnnotAll) <- row.names(junc.RPM)

    # get GenomicRanges object with genes for whole genome
    gn <- genes(get(TxDb))

    gSymbol <- select(get(annotation),
        keys = as.character(gn$gene_id),
        columns = c("SYMBOL"), keytype = "ENTREZID"
    )
    gn$SYMBOL <- gSymbol$SYMBOL
    gn$ENTREZID <- gSymbol$ENTREZID
    rm(gSymbol, chr, end, start)
    gc()

    # find symbols and ENTREZID for our junctions
    overlap <- findOverlaps(geneAnnot, gn)

    geneSYMBOLS <- tapply(gn$SYMBOL[subjectHits(overlap)],
        queryHits(overlap), paste,
        collapse = ";"
    )

    # add to genome ranges object
    geneAnnot$SYMBOL <- NA
    geneAnnot$SYMBOL[as.numeric(names(geneSYMBOLS))] <- geneSYMBOLS

    geneENTREZID <- tapply(gn$ENTREZID[subjectHits(overlap)],
        queryHits(overlap), paste,
        collapse = ";"
    )

    geneAnnot$ENTREZID <- NA
    geneAnnot$ENTREZID[as.numeric(names(geneENTREZID))] <- geneENTREZID

    ### Make sure the genes are ordered by lower to higher ENTREZID
    geneAnnot$ENTREZID <- vapply(geneAnnot$ENTREZID, function(x)
        paste(sort(as.numeric(strsplit(x, ";")[[1]])), collapse = ";"), character(1))

    # get known exons
    en <- exons(get(TxDb))

    message("filtering for putative junctions, skipping, insertion and deletion events based on known exons.")
    # annotate skipping events as those which overlap with more than two exons
    exo <- findOverlaps(
        GRanges(seqnames(geneAnnot),
            ranges = IRanges(
                start = start(geneAnnot) + 1,
                end = end(geneAnnot) - 1
            ),
            strand = strand(geneAnnot)
        ), en,
        type = "within"
    )
    exoCount <- table(factor(queryHits(exo))) > 0
    geneAnnot$skipping <- FALSE
    geneAnnot$skipping[as.numeric(names(exoCount))] <- exoCount
    message("skipping")

    # insertion events start or end outside of known exons
    geneAnnot.start <- GRanges(seqnames(geneAnnot),
        ranges = IRanges(
            start = start(geneAnnot),
            end = start(geneAnnot)
        ),
        strand = strand(geneAnnot)
    )
    geneAnnot.end <- GRanges(seqnames(geneAnnot),
        ranges = IRanges(
            start = end(geneAnnot),
            end = end(geneAnnot)
        ),
        strand = strand(geneAnnot)
    )

    geneAnnot$insertions <- !overlapsAny(geneAnnot.start, en) |
        !overlapsAny(geneAnnot.end, en)

    message("insertions")

    # deletion events occur w/in an exon but not at its start or end
    en.start <- GRanges(seqnames(en),
        ranges = IRanges(
            start = start(en),
            end = start(en)
        ),
        strand = strand(en)
    )
    en.end <- GRanges(seqnames(en),
        ranges = IRanges(
            start = end(en),
            end = end(en)
        ),
        strand = strand(en)
    )
    geneAnnot$deletions <-
        (!overlapsAny(geneAnnot.end, en.start) & overlapsAny(geneAnnot.end, en)) |
            (!overlapsAny(geneAnnot.start, en.end) & overlapsAny(geneAnnot.start, en))

    message("deletions")

    return(geneAnnot)
}
