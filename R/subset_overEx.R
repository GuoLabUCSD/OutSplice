# Function to subset significant over-expressed junction events

subset_overEx <- function(FisherAnalyses, topgenelist90.ogsa, p_value) {
  toplist90 <- FisherAnalyses[topgenelist90.ogsa, ]
  if (sum(toplist90[, "FisherP2"] < p_value) > 0) {
    toplist90 <- toplist90[toplist90[, "FisherP2"] < p_value, ]
    toplist90 <- toplist90[, "FisherP2"]
  } else {
    toplist90 <- vector()
    print("no over expression outliers")
  }
  return(toplist90)
}
