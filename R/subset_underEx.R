# Function to subset significant under-expressed junction events

subset_underEx <- function(FisherAnalyses, topgenelist10.ogsa, p_value) {
  toplist10 <- FisherAnalyses[topgenelist10.ogsa, ]
  if (sum(toplist10[, "FisherP1"] < p_value) > 0) {
    toplist10 <- toplist10[toplist10[, "FisherP1"] < p_value, ]
    toplist10 <- toplist10[, "FisherP1"]
  } else {
    toplist10 <- vector()
    print("no under expression outliers")
  }
  return(toplist10)
}
