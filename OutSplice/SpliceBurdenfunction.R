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

  Number_overexpressed_per_tumor <- sapply(Over_Expressed_Calls_with_Type, table)
  Number_overexpressed_per_tumor <- as.data.frame(Number_overexpressed_per_tumor)

  transposed_Over <- t(Number_overexpressed_per_tumor)
  colnames(transposed_Over) <- c('NoOutlierNumber', 'OutlierNumberOver')
  transposed_Over <- subset(transposed_Over, select = -NoOutlierNumber)
  Over_Expressed_Splicing_Burden <- as.data.frame(transposed_Over)

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

  Number_underexpressed_per_tumor <- sapply(Under_Expressed_Calls_with_Type, table)
  Number_underexpressed_per_tumor <- as.data.frame(Number_underexpressed_per_tumor)

  transposed_Under <- t(Number_underexpressed_per_tumor)
  colnames(transposed_Under) <- c('NoOutlierNumber', 'OutlierNumberUnder')
  transposed_Under <- subset(transposed_Under, select = -NoOutlierNumber)
  Under_Expressed_Splicing_Burden <- as.data.frame(transposed_Under)

  #Total Splicing Burden as the sum of Over + Under Expressed Events
  Total_Outliers <- transform(merge(transposed_Over, transposed_Under, by = 0, all = TRUE), TotalOutliers = OutlierNumberOver + OutlierNumberUnder)
  rownames(Total_Outliers) <- Total_Outliers[, 1]
  transposed_Total <- subset(Total_Outliers, select = c(-OutlierNumberOver, -OutlierNumberUnder, -Row.names))
  Total_Splicing_Burden <- as.data.frame(transposed_Total)

  detach("package:dplyr")

  return(splice_burden <- cbind(Over_Expressed_Splicing_Burden, Under_Expressed_Splicing_Burden, Total_Splicing_Burden))
}
