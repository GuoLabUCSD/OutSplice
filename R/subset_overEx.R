# Function to subset significant over-expressed junction events

subset_overEx <- function(FisherAnalyses, p_value) {
    toplist90 <- FisherAnalyses[order(FisherAnalyses$OE_Rank), ]
    if (sum(toplist90[, "FisherP2"] < p_value) > 0) {
        toplist90 <- toplist90[toplist90[, "FisherP2"] < p_value, ]
        toplist90 <- subset(toplist90, select = FisherP2)
    } else {
        toplist90 <- NULL
    }
    return(toplist90)
}
