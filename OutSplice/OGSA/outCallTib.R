#  Outlier Gene Set Analysis
#
#  cite with PMID: 24277953
#  Ochs et al IEEE/ACM Trans Comput Biol Bioinform 2013
#
#   outCallTib - counts outliers by the Tibshirani and Hastie method from
#       Biostatistics, 8, 2, 2006
#   and generates a list object with all outliers noted
#
#
# Inputs: dataList - set of data which includes phenotype information
#         tail - vector equal to number of matrices with
#                  values 'left' or 'right' for where to find
#                  outliers
#         corr = whether to correct for normal outliers
#             ONLY for compatibility, since method does not allow
#             determining specific changes in cases, it will just
#             print message if corr = TRUE
#         names - vector equal to number of matrices to name
#                 molecular type of data (e.g., 'CNV')
#
#  Output: a list with all specific outlier calls for each molecular
#          type in each case sample
#
#  v 1.0 - Michael Ochs, TCNJ, 15 Jan 2014


outCallTib <- function(dataList, tail='right', corr=FALSE, names=NULL) {

    if (is.null(names)) {
        names <- vector(length=length(dataList),mode='character')
        for (d in 1:length(dataList)) {
            names[d] <- paste('Data',d)
        }
        
    }
    
    if (corr) {
        print("No Correction Applicable in Tibshirani Indicators")
    }
    temp <- dataList[[1]]
    temp <- temp[[1]]
    nG <- dim(temp)[1]
    nT <- dim(temp)[2]
    outList <- list()

    for (d in 1:length(dataList)) {
        data <- dataList[[d]]
        phenotype <- data[[2]]
        data <- data[[1]]
        nS <- length(phenotype)
        thisTail <- tail[d]
        
        nData <- data[,phenotype==0]
        
        # generate Tibshirani eqn 2.2
        # putting genes on same scale
        medG <- apply(data, 1, median)
        madGN <- apply(nData, 1, mad)

        for (i in 1:nG) {
            temp <- (data[i,] - medG[i])/madGN[i]
            data[i,] <- temp
        }
        tData <- data[,phenotype==1]
        nData <- data[,phenotype==0]
        nT <- dim(tData)[2]

        # count outliers relative to scaled data
        iqrG <- apply(data, 1, IQR)
        if (thisTail == 'right') {
            quantG <- apply(data, 1, quantile, probs = 0.75)
        } else if (thisTail == 'left') {
            quantG <- apply(data, 1, quantile, probs = 0.25)
        }
        
        # generate indicators and place in output matrix
        empirP <- matrix(nrow=nG, ncol=nT)
        if (thisTail == 'right') {
            for (i in 1:nG) {
                Ind <- tData[i,] > (quantG[i] + iqrG[i])
                empirP[i,] <- Ind
            }
        } else if (thisTail == 'left') {
            for (i in 1:nG) {
                Ind <- tData[i,] < (quantG[i] - iqrG[i])
                empirP[i,] <- Ind
            }
        }
        rownames(empirP) <- rownames(data)
        colnames(empirP) <- colnames(tData)
        outList[[d]] <- empirP
    }

    names(outList) <- names
	return(outList)
}		
	

