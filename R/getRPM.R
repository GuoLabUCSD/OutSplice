# Function to normalize by Total Read Counts and remove junctions with no difference between tumor and normal

getRPM <- function(rawcounts, all.junc, all.samples) {
    ## get data in RPM
    message("convert to RPM")
    totalrawcount <- colSums(rawcounts)
    names(totalrawcount) <- colnames(rawcounts)
    junc.RPM <- apply(all.junc[, all.samples], 1, function(x) {
        x / totalrawcount[all.samples] * 1000000
    })
    junc.RPM <- t(junc.RPM)

    message("filter the putative junctions")
    ### filter genes by cut off of
    fcCutoff <- 10
    Cutoff.ratio <- 1 - (1 / fcCutoff)
    # enforce overall Fold change
    logFC <- apply(junc.RPM, 1, function(x) {
        (max(x) - min(x)) / max(x)
    }) > Cutoff.ratio
    junc.RPM <- junc.RPM[which(logFC), ]

    rm(logFC)
    gc()

    return(junc.RPM)
}
