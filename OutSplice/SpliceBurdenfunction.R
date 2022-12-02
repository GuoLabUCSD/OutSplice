#Function to Calculate Splicing Burden Based on the Number of Fisher Significant Events in Each Sample

CalcBurden <- function(junc.Outliers, ASE.type, FisherAnalyses, p_value) {

  library(dplyr)

  #Calculate Splicing Burden for Significant Events Over-Expressed in Tumors
  Over_Expressed_Junctions_Outlier_Calls <- as.data.frame(junc.Outliers$TumorOverExpression)

  ASE.typedf <- as.data.frame(ASE.type)
  row.names(ASE.typedf) <- sub('\\.', ':', rownames(ASE.typedf))
  row.names(ASE.typedf) <- sub('\\.', '-', rownames(ASE.typedf))

  Over_Expressed_Calls_with_Type <- merge(ASE.typedf, Over_Expressed_Junctions_Outlier_Calls, by='row.names', all.x=TRUE)
  rownames(Over_Expressed_Calls_with_Type) <- Over_Expressed_Calls_with_Type[, 1]
  Over_Expressed_Calls_with_Type <- subset(Over_Expressed_Calls_with_Type, select = -Row.names)
  Over_Expressed_Calls_with_Type <- na.omit(Over_Expressed_Calls_with_Type)
  Over_Expressed_Calls_with_Type <- subset(Over_Expressed_Calls_with_Type, select = c(-skipping, -insertions, -deletions))

  OE_pvalue <- subset(FisherAnalyses, select = FisherP2)
  Over_Expressed_Calls_with_Type <- merge(OE_pvalue, Over_Expressed_Calls_with_Type, by='row.names', all.x=TRUE)
  rownames(Over_Expressed_Calls_with_Type) <- Over_Expressed_Calls_with_Type[, 1]
  Over_Expressed_Calls_with_Type <- subset(Over_Expressed_Calls_with_Type, select = -Row.names)
  Over_Expressed_Calls_with_Type <- na.omit(Over_Expressed_Calls_with_Type)
  Over_Expressed_Calls_with_Type <- subset(Over_Expressed_Calls_with_Type, Over_Expressed_Calls_with_Type$FisherP2 < p_value)
  Over_Expressed_Calls_with_Type <- subset(Over_Expressed_Calls_with_Type, select = -FisherP2)
  
  OutlierNumberOver = colSums(Over_Expressed_Calls_with_Type)
  all_oe_results_df = data.frame(OutlierNumberOver)

  #Calculate Splicing Burden for Significant Events Under-Expressed in Tumors
  Under_Expressed_Junctions_Outlier_Calls <- as.data.frame(junc.Outliers$TumorUnderExpression)

  ASE.typedf <- as.data.frame(ASE.type)
  row.names(ASE.typedf) <- sub('\\.', ':', rownames(ASE.typedf))
  row.names(ASE.typedf) <- sub('\\.', '-', rownames(ASE.typedf))

  Under_Expressed_Calls_with_Type <- merge(ASE.typedf, Under_Expressed_Junctions_Outlier_Calls, by='row.names', all.x=TRUE)
  rownames(Under_Expressed_Calls_with_Type) <- Under_Expressed_Calls_with_Type[, 1]
  Under_Expressed_Calls_with_Type <- subset(Under_Expressed_Calls_with_Type, select = -Row.names)
  Under_Expressed_Calls_with_Type <- na.omit(Under_Expressed_Calls_with_Type)
  Under_Expressed_Calls_with_Type <- subset(Under_Expressed_Calls_with_Type, select = c(-skipping, -insertions, -deletions))

  UE_pvalue <- subset(FisherAnalyses, select = FisherP1)
  Under_Expressed_Calls_with_Type <- merge(UE_pvalue, Under_Expressed_Calls_with_Type, by='row.names', all.x=TRUE)
  rownames(Under_Expressed_Calls_with_Type) <- Under_Expressed_Calls_with_Type[, 1]
  Under_Expressed_Calls_with_Type <- subset(Under_Expressed_Calls_with_Type, select = -Row.names)
  Under_Expressed_Calls_with_Type <- na.omit(Under_Expressed_Calls_with_Type)
  Under_Expressed_Calls_with_Type <- subset(Under_Expressed_Calls_with_Type, Under_Expressed_Calls_with_Type$FisherP1 < p_value)
  Under_Expressed_Calls_with_Type <- subset(Under_Expressed_Calls_with_Type, select = -FisherP1)
  
  OutlierNumberUnder = colSums(Under_Expressed_Calls_with_Type)
  all_ue_results_df = data.frame(OutlierNumberUnder)
  
  #Total Splicing Burden as the sum of Over + Under Expressed Events
  total_results = transform(merge(all_oe_results_df, all_ue_results_df, by = 0, all = TRUE), TotalOutliers = OutlierNumberOver + OutlierNumberUnder)
  rownames(total_results) <- total_results[, 1]
  convert_total <- subset(total_results, select = c(-OutlierNumberOver, -OutlierNumberUnder, -Row.names))
  total_burden <- as.data.frame(convert_total)

  detach("package:dplyr")

  return(splice_burden <- cbind(all_oe_results_df, all_ue_results_df, total_burden))
}
