#  Outlier Gene Set Analysis
# 
#  cite with PMID: 24277953
#  Ochs et al IEEE/ACM Trans Comput Biol Bioinform 2013
#
#  copaInt - counts outliers by Tibshirani-Hastie method
#      by calling outCount after setting up list
#      or by rank outlier method by calling outRank
# 
#  Inputs: dataList - set of matrices of molecular data
#          phenotype - vector of 1 for case, 0 for control
#          tails - vector equal to number of matrices with
#                  values 'left' or 'right' for where to find
#                  outliers
#          thres = alpha value
#          method = 'Tibshirani', 'Rank'
#          corr = whether to correct for normal outliers
#          offsets = vector equal to number of matrices which
#                    sets minimum value relative to normal to
#                    call outlier (corrected rank only)
# 
#  Output: a vector with outlier counts by gene
# 
#  v 1.0 - Michael Ochs, TCNJ, 15 Jan 2014

# source("outRank.R")
# source("outCount.R")

copaInt <- function(dataList, phenotype, tails, thres = 0.05, method='Tibshirani', corr=FALSE, offsets=NULL) {

    nType <- length(dataList)
    nGene <- dim(dataList[[1]])[1]
    
    # reorder the columns in all data matrices
    # to match
    for (i in 1:nType) {
        tempMat <- dataList[[i]]
        temp <- match(names(phenotype),colnames(tempMat))
        tempMat <- tempMat[,temp]
        dataList[[i]] <- tempMat
    }
    
    # reorder the rows in all data matrices
    # to match
    geneNames <- rownames(dataList[[1]])
    missingData <- rep(FALSE,length(geneNames))
    for (i in 1:nType) {
        tempMat <- dataList[[i]]
        temp <- match(geneNames,rownames(tempMat))
        tempMat <- tempMat[temp,]
        dataList[[i]] <- tempMat
        missingData[is.na(tempMat[,1])] <- TRUE
    }
    
    
    # remove rows that are not in all matrices
    # for counting, all genes must be in each matrix
    for (i in 1:nType) {
        tempMat <- dataList[[i]]
        tempMat <- tempMat[!missingData,]
        dataList[[i]] <- tempMat
    }
    nGene <- dim(tempMat)[1]
    outCts <- rep(0,nGene)
    names(outCts) <- rownames(tempMat)
    
    # remove missing data and run Tibshirani method
    # counting outliers
    if (method == 'Tibshirani') {
        for (i in 1:nType) {
            dataMat <- dataList[[i]]
            noData <- is.na(dataMat[1,])
            dataMat <- dataMat[,!noData]
            phenoMat <- phenotype[!noData]
            temp <- outCount(dataMat,phenoMat,tail=tails[i],corr=corr)
            outCts <- outCts + temp
        }
    } else if (method == 'Rank') {
        for (i in 1:nType) {
            dataMat <- dataList[[i]]
            noData <- is.na(dataMat[1,])
            dataMat <- dataMat[,!noData]
            phenoMat <- phenotype[!noData]
            dataList[[i]] <- list(exprs=dataMat,classlab=phenoMat)
        }
        outCts <- outRank(dataList,thres = thres, tail=tails,corr=corr,offsets=offsets)
    } else {
        print('Unknown Outlier Counting Method; Returning Zeros')
    }

    return(outCts)

}