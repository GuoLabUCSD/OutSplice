#  Outlier Gene Set Analysis
#
#  cite with PMID: 24277953
#  Ochs et al IEEE/ACM Trans Comput Biol Bioinform 2013
#
#  outCallRank - counts outliers by the Ghosh method from
#         J Biopharm Stat 20, 193, 2010
#      and generates list object with all outliers noted,
#      modified from Yingying Wei's original
#
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
#         names - vector equal to number of matrices to name
#                 molecular type of data (e.g., 'CNV')
#
#  Output: a list with all specific outlier calls for each molecular
#          type in each case sample
#
#  v 1.0 - Michael Ochs, TCNJ, 15 Jan 2014
# v 1.5 - Theresa Guo, 7.29.15 add additional "filter" term, to filter out unwanted samples


outCallRankFilter <- function(dataList, filter=NULL, thres= 0.05, tail='right', corr=FALSE, offsets=NULL,names=NULL) {
    
    if (is.null(names)) {
        names <- vector(length=length(dataList),mode='character')
        for (d in 1:length(dataList)) {
            names[d] <- paste('Data',d)
        }
    
    }
    
    temp <- dataList[[1]]
    temp <- temp[[1]]
    nG <- dim(temp)[1]
    outList <- list()

    ## if no filter is specified, make a matrix of all ones (i.e. everything included)
    if (is.null(filter)){
      filter<-temp
      filter<-filter*0+1}
    
    if (!corr) {
        offsets <- rep(0.0,dim(temp)[2])
    } else if (is.null(offsets)) {
        print('No Offsets Set with Correction Requested')
        return()
    }

    outP <- matrix(nrow=nG, ncol=2)
    outCount <- rep(0,nG)

    for (d in 1:length(dataList)) {
        data <- dataList[[d]]
        phenotype <- data[[2]]
        data <- data[[1]]
        ## number of samples
        nS <- length(phenotype)
        thisTail <- tail[d]
        
        adjust <- offsets[d]
        nData <- data[,phenotype==0]
        tData <- data[,phenotype==1]
        nT <- dim(tData)[2] #number of tumors
        
        ## generate filtering matrix for normals and tumors
        nfilter<-filter[,phenotype==0]
        tfilter<-filter[,phenotype==1]
        
        # generate empicical pValues as the
        # number of sum(normals{<,>}tumor)/nN
        empirP <- matrix(nrow=nG, ncol=nT)
        rownames(empirP) <- rownames(data)
        colnames(empirP) <- colnames(tData)
        if (thisTail == 'right') {
            for (i in 1:nG) {
                tumor <- tData[i,]
                ## filter out normals that don't meet criteria
                baseline <- nData[i,which(nfilter[i,]==1)]
                result <- sapply(1:length(tumor),function(j)
                  sum((baseline+adjust)>tumor[j]))
                empirP[i,] <- result/length(baseline)
                ## outliers called here
                empirP[i,] <- empirP[i,] < thres
                ## eliminate tumors at the end (with 0 - cannot be outlier), if they are to be filtered
                empirP[i,tfilter[i,]==0]<-NA

            }
        } else if (thisTail == 'left') {
            for (i in 1:nG) {
                tumor <- tData[i,]
                baseline <- nData[i,which(nfilter[i,]==1)]
                result <- sapply(1:length(tumor),function(j)
                    sum((baseline-adjust)<tumor[j]))
                empirP[i,] <- result/length(baseline)
                ## outliers called here
                empirP[i,] <- empirP[i,] < thres
                empirP[i,tfilter[i,]==0]<-NA
            }
        }

        outList[[d]] <- empirP
    }
    names(outList) <- names
	return(outList)
}		
	

