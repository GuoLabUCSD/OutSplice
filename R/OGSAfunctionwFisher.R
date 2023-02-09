# -- Function: OGSA for ranking fusions or junctions within the data###
## First version Nov 14,2014
## This function ranks samples based on the number of outliers present in tumor samples, and then by median of expression in tumors
## Add Fisher test to compare outlier calls in tumors and normals. December 17, 2014
## Add FDR correction as a parameter, April 17, 2014.  default, no correction.
# Otherwise corection==name of correction method, i.e. "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
# if getting outlier calls for output, outliers=T

dotheogsa<-function(Sample.data, PHENO, offsets=0.001, Fisher=F, correction=NA, outliers=F, dir=NA)
  ## default value for offsets is 0.001
#   PHENO<-pheno$classes=='Tumor'
# PHENO<-as.numeric(PHENO)
# names(PHENO)<-row.names(pheno)
  ##PHENO should have 'Normal' or 'Tumor' calls where Tumor ==1, Normal ==0, and names of each sample associated
{
# libraries on which OGSA depends
#library('limma')

# data required by the package
# source code
#source(paste0(dir, 'OGSA/copaStat.R'))
#source(paste0(dir, 'OGSA/outCount.R'))
#source(paste0(dir, 'OGSA/outRank.R'))
#source(paste0(dir, 'OGSA/outCallRank.R'))
#source(paste0(dir, 'OGSA/outCallTib.R'))
#source(paste0(dir, 'OGSA/copaInt.R'))
#print('yes!')

## Define tumors and normals
## take fus.RPM (final filtered list) and pheno

### run outlier analysis for true distribution of tumors and normals

if (Fisher==T){
## Initialize a matrix so program will call outliers in normals
Sample.data.normals<-Sample.data[,PHENO==0]
no.of.normals<-length(PHENO[PHENO==0])
no.of.tumors<-length(PHENO[PHENO==1])
# repeat the normal samples so they will be tested as if they are a set of tumors
Sample.data.normals<-cbind(Sample.data.normals, Sample.data.normals)
PHENO.N<-PHENO[PHENO==0]
PHENO.N[PHENO.N==0]<-1
PHENO2<-c(PHENO[PHENO==0], PHENO.N)
dataExprs2<- list(list(Sample.data.normals, PHENO2))
}

  PHENO<-t(data.frame(PHENO))

dataExprs <- list(list(Sample.data,PHENO))


## get expression of junction in tumors
Sample.data.tumor<-Sample.data[,PHENO==1] ## subset to expression in tumors only
SampleTotals<-apply(Sample.data.tumor, 1, median) ## median of tumor expression

## Rank first by number of outliers identified.  Then if there are samples with the same outlier, order by the median of tumor expression

## smallest to largest, copa10 (underexpression in tumors)
## now do the outlier ranking with tail='left' for if tumor is less than normal
outRankExprs1 <- outCallRank(dataExprs,names=c('Expr'),tail='left',corr=T,offsets=offsets)$Expr
### get number of outliers in the tumor group
outRankTumor1 <- apply(outRankExprs1,1,sum)
SampleTotals1<-order(SampleTotals, decreasing=T)
## ranks first by # of outliers, then by smaller level of tumor expression
var1<-order(outRankTumor1, SampleTotals1, decreasing=T)

if (Fisher==T){
## Calculate outliers in normal tissue
outRankExprs1.normals<-outCallRank(dataExprs2,names=c('Expr'),tail='left',corr=T,offsets=offsets)$Expr
outRankNormal1 <- apply(outRankExprs1.normals,1,sum)

## Perform Fisher's exact test
FisherP1<-matrix(NA, nrow=nrow(Sample.data), ncol=1, dimnames=list(c(row.names(Sample.data)), 'FisherP1'))

for (i in 1:nrow(Sample.data)){
  ## Make contingency tables
  test<-matrix(c(outRankNormal1[i],no.of.normals-outRankNormal1[i], outRankTumor1[i],no.of.tumors-outRankTumor1[i]), nrow=2, ncol=2, dimnames=list(c('outliers','rest'), c('Normal', 'Tumor')))
  ## Get p-value from fisher test.  Will do two sided test
  FisherP1[i]<-fisher.test(test, alternative='two.sided')$p.value
}

}

###### largest to smallest, copa90 (overexpression in tumors)#####################
outRankExprs2 <- outCallRank(dataExprs,names=c('Expr'),tail='right',corr=T,offsets=offsets)$Expr
outRankTumor2 <- apply(outRankExprs2,1,sum)
var2<-order(outRankTumor2, SampleTotals, decreasing = T)


if (Fisher==T) {
## Calculate outliers in normal tissue
outRankExprs2.normals<-outCallRank(dataExprs2,names=c('Expr'),tail='right',corr=T,offsets=offsets)$Expr
outRankNormal2 <- apply(outRankExprs2.normals,1,sum)

## Perform Fisher's exact test
FisherP2<-matrix(NA, nrow=nrow(Sample.data), ncol=1, dimnames=list(c(row.names(Sample.data)), 'FisherP2'))

for (i in 1:nrow(Sample.data)){
  ## Make contingency tables
  test<-matrix(c(outRankNormal2[i],no.of.normals-outRankNormal2[i], outRankTumor2[i],no.of.tumors-outRankTumor2[i]), nrow=2, ncol=2, dimnames=list(c('outliers','rest'), c('Normal', 'Tumor')))
  ## Get p-value from fisher test.  Will do two sided test
  FisherP2[i]<-fisher.test(test, alternative='two.sided')$p.value
}
}

##output data for the function!
if (outliers==F){

if (Fisher==T){
  ## add FDR correction
  if (!is.na(correction)){
  FisherP1<-p.adjust(FisherP1, method = correction)
  FisherP2<-p.adjust(FisherP2, method = correction)
  }
return(RankingData<-cbind(outRankTumor1, var1, FisherP1, outRankTumor2, var2, FisherP2))
} else   {return(RankingData<-cbind(outRankTumor1, var1, outRankTumor2, var2))
}

}else {return(OutlierCalls<-list("TumorUnderExpression" = outRankExprs1, "TumorOverExpression" = outRankExprs2))}



}
