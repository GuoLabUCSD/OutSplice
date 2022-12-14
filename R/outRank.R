#  Outlier Gene Set Analysis
#
#    this routine is generally called by copaInt not used directly
#
#  cite with PMID: 24277953
#  Ochs et al IEEE/ACM Trans Comput Biol Bioinform 2013
#
#  outCount - counts outliers by the Ghosh method from
#         J Biopharm Stat 20, 193, 2010
#      from Y Wei original method
#
# Inputs: dataList - set of data which includes phenotype information
#         thres - alpha value
#         tail - vector equal to number of matrices with
#                  values 'left' or 'right' for where to find
#                  outliers
#         corr = whether to correct for normal outliers
#         offsets = vector equal to number of matrices which
#                    sets minimum value relative to normal to
#                    call outlier (corrected rank only)
#
#  Output: a vector with outlier counts by gene
#
#  v 1.0 - Michael Ochs, TCNJ, 15 Jan 2014

outRank <- function(dataList, thres= 0.05, tail='right', corr=FALSE, offsets=NULL) {
    
    temp <- dataList[[1]]
    temp <- temp[[1]]
    nG <- dim(temp)[1]

    if (!corr) {
        offsets <- rep(0.0,dim(temp)[2])
    } else if (is.null(offsets)) {
        print('No Offsets Set with Correction Requested')
        return()
    }

    #    nT <- 0
    #    for (i in 1:length(dataList)) {
    #        temp <- dim(dataList[[1]])
    #        nT <- nT + length(temp[[2]]==1)
    #    }
    outP <- matrix(nrow=nG, ncol=2)
    outCount <- rep(0,nG)

    for (d in 1:length(dataList)) {
        data <- dataList[[d]]
        phenotype <- data[[2]]
        data <- data[[1]]
        nS <- length(phenotype)
        thisTail <- tail[d]
        
        adjust <- offsets[d]
        nData <- data[,phenotype==0]
        tData <- data[,phenotype==1]
        nT <- dim(tData)[2]
        
        # generate empicical pValues as the
        # number of sum(normals{<,>}tumor)/nN
        empirP <- matrix(nrow=nG, ncol=nT)
        if (thisTail == 'right') {
            for (i in 1:nG) {
                tumor <- tData[i,]
                baseline <- nData[i,]
                result <- sapply(1:length(tumor),function(j)
                    sum((baseline+adjust)>tumor[j]))
                empirP[i,] <- result/length(baseline)
            }
        } else if (thisTail == 'left') {
            for (i in 1:nG) {
                tumor <- tData[i,]
                baseline <- nData[i,]
                result <- sapply(1:length(tumor),function(j)
                    sum((baseline-adjust)<tumor[j]))
                empirP[i,] <- result/length(baseline)
            }
        }
        outP <- cbind(outP,empirP)
    }
    yP <- dim(outP)[2]
    outP <- outP[,3:yP]
    
    # count the number of genes that have outliers based
    # on desired threshold
    for (i in 1:nG) {
        outCount[i] <- sum(outP[i,] < thres)
    }

    names(outCount) <- rownames(data)
	return(outCount)
}		
	

