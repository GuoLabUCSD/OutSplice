#  Outlier Gene Set Analysis
#
#    this routine is generally called by copaInt not used directly
#
#  cite with PMID: 24277953
#  Ochs et al IEEE/ACM Trans Comput Biol Bioinform 2013
#
#  outCount - counts outliers by the Tibshirani and Hastie method from
#         Biostatistics, 8, 2, 2006
#     adds ability to subtract for outliers in the normals
#     using corr = TRUE
#
# Inputs: dataList - set of data matrices
#         phenotype - vector of 1 for case, 0 for control
#         tail - vector equal to number of matrices with
#                  values 'left' or 'right' for where to find
#                  outliers
#         corr = whether to correct for normal outliers
#
#  Output: a vector with outlier counts by gene
#
#  v 1.0 - Michael Ochs, TCNJ, 15 Jan 2014


outCount <- function(data, phenotype, tail='right', corr=FALSE) {

	nG <- dim(data)[1]
	nS <- length(phenotype)

	outCount <- vector(length = nG)

	nData <- data[,phenotype==0]

    # generate Tibshirani eqn 2.2
    # putting genes on same scale
	medG <- apply(data, 1, median)
	madGN <- apply(nData, 1, mad)

    for (i in 1:nG) {
		temp <- (data[i,] - medG[i])/madGN[i]
        data[i,] <- temp
	}

    # count outliers relative to scaled data
	tData <- data[,phenotype==1]
	nData <- data[,phenotype==0]
    iqrG <- apply(data, 1, IQR)
    if (tail == 'right') {
        quantG <- apply(data, 1, quantile, probs = 0.75)
    } else if (tail == 'left') {
        quantG <- apply(data, 1, quantile, probs = 0.25)
    }
    
    if (tail == 'right') {
        for (i in 1:nG) {
            Ind <- tData[i,] > quantG[i] + iqrG[i]
            outCount[i] <- sum(Ind)
        }
    } else if (tail == 'left') {
        for (i in 1:nG) {
            Ind <- tData[i,] < quantG[i] - iqrG[i]
            outCount[i] <- sum(Ind)
        }
    }
    
    # if correcting, subtract number of normals
    # that are called outliers
    if (corr) {
        if (tail == 'right') {
            for (i in 1:nG) {
                Ind <- nData[i,] > quantG[i] + iqrG[i]
                outCount[i] <- outCount[i] - sum(Ind)
            }
        } else if (tail == 'left') {
            for (i in 1:nG) {
                Ind <- nData[i,] < quantG[i] - iqrG[i]
                outCount[i] <- outCount[i] - sum(Ind)
            }
        }
    }


    names(outCount) <- rownames(data)
	return(outCount)
}		
	

