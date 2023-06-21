# Function to subset significant under-expressed junction events

subset_underEx <- function(FisherAnalyses, p_value) {
    toplist10 <- FisherAnalyses[order(FisherAnalyses$UE_Rank), ]
    if (sum(toplist10[, "FisherP1"] < p_value) > 0) {
        toplist10 <- toplist10[toplist10[, "FisherP1"] < p_value, ]
        toplist10 <- subset(toplist10, select = FisherP1)
    } else {
        toplist10 <- NULL
    }
    return(toplist10)
}
