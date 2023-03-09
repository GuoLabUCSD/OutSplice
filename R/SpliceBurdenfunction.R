#Function to Calculate Splicing Burden Based on the Number of Fisher Significant Events in Each Sample

calcBurden <- function(junc.Outliers, FisherAnalyses, p_value) {

  #Calculate Splicing Burden for Significant Events Over-Expressed in Tumors
  Over_Expressed_Junctions_Outlier_Calls <- as.data.frame(junc.Outliers$TumorOverExpression)
  row.names(Over_Expressed_Junctions_Outlier_Calls) <- sub('\\.', ':', rownames(Over_Expressed_Junctions_Outlier_Calls))
  row.names(Over_Expressed_Junctions_Outlier_Calls) <- sub('\\.', '-', rownames(Over_Expressed_Junctions_Outlier_Calls))

  OE_pvalue <- subset(FisherAnalyses, select = FisherP2)
  Over_Expressed_Junctions_Outlier_Calls <- merge(OE_pvalue, Over_Expressed_Junctions_Outlier_Calls, by='row.names', all.x=TRUE)
  rownames(Over_Expressed_Junctions_Outlier_Calls) <- Over_Expressed_Junctions_Outlier_Calls[, 1]
  Over_Expressed_Junctions_Outlier_Calls <- subset(Over_Expressed_Junctions_Outlier_Calls, select = -Row.names)
  Over_Expressed_Junctions_Outlier_Calls <- na.omit(Over_Expressed_Junctions_Outlier_Calls)
  Over_Expressed_Junctions_Outlier_Calls <- subset(Over_Expressed_Junctions_Outlier_Calls, Over_Expressed_Junctions_Outlier_Calls$FisherP2 < p_value)
  Over_Expressed_Junctions_Outlier_Calls <- subset(Over_Expressed_Junctions_Outlier_Calls, select = -FisherP2)

  OutlierNumberOver = colSums(Over_Expressed_Junctions_Outlier_Calls)
  all_oe_results_df = data.frame(OutlierNumberOver)

  #Calculate Splicing Burden for Significant Events Under-Expressed in Tumors
  Under_Expressed_Junctions_Outlier_Calls <- as.data.frame(junc.Outliers$TumorUnderExpression)
  row.names(Under_Expressed_Junctions_Outlier_Calls) <- sub('\\.', ':', rownames(Under_Expressed_Junctions_Outlier_Calls))
  row.names(Under_Expressed_Junctions_Outlier_Calls) <- sub('\\.', '-', rownames(Under_Expressed_Junctions_Outlier_Calls))

  UE_pvalue <- subset(FisherAnalyses, select = FisherP1)
  Under_Expressed_Junctions_Outlier_Calls <- merge(UE_pvalue, Under_Expressed_Junctions_Outlier_Calls, by='row.names', all.x=TRUE)
  rownames(Under_Expressed_Junctions_Outlier_Calls) <- Under_Expressed_Junctions_Outlier_Calls[, 1]
  Under_Expressed_Junctions_Outlier_Calls <- subset(Under_Expressed_Junctions_Outlier_Calls, select = -Row.names)
  Under_Expressed_Junctions_Outlier_Calls <- na.omit(Under_Expressed_Junctions_Outlier_Calls)
  Under_Expressed_Junctions_Outlier_Calls <- subset(Under_Expressed_Junctions_Outlier_Calls, Under_Expressed_Junctions_Outlier_Calls$FisherP1 < p_value)
  Under_Expressed_Junctions_Outlier_Calls <- subset(Under_Expressed_Junctions_Outlier_Calls, select = -FisherP1)

  OutlierNumberUnder = colSums(Under_Expressed_Junctions_Outlier_Calls)
  all_ue_results_df = data.frame(OutlierNumberUnder)

  #Total Splicing Burden as the sum of Over + Under Expressed Events
  total_results = transform(merge(all_oe_results_df, all_ue_results_df, by = 0, all = TRUE), TotalOutliers = OutlierNumberOver + OutlierNumberUnder)
  rownames(total_results) <- total_results[, 1]
  convert_total <- subset(total_results, select = c(-OutlierNumberOver, -OutlierNumberUnder, -Row.names))
  total_burden <- as.data.frame(convert_total)

  return(splice_burden <- cbind(all_oe_results_df, all_ue_results_df, total_burden))
}
