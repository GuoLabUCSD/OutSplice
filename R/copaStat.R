
# outlier stats by the Tibshirani and Hastie method from
# Biostatistics, 8, 2, 2006
# data is matrix of nGene x nSample
# phenotype is vector of {0,1} of length nSample, 1 = case
# tail indicates whether up (right) or down (left) outliers
# perms is number of permuations
# permType is by all on array or by gene, if by gene increase
# perms significantly and plan on lots of time; in theory
# array should be fine as genes are rescaled

copaStat <- function(data, phenotype, tail='right', perms=100, permType='array') {

	nG <- dim(data)[1]
	nS <- length(phenotype)
    nT <- sum(phenotype)
    
	nData <- data[,phenotype==0]

    dataC <- matrix(nrow=nG,ncol=nS)
    
    copaValue <- vector(length=nG)
    copaStat <- vector(length=nG)
    copaPerm <- matrix(nrow=nG,ncol=perms)
    

    # generate permutation COPA values
	for (j in 1:perms) {
		phenoRand <- rep(0, nS)
		phenoRand[sample(1:nS,nT)] <- 1
        
 
        # generate Tibshirani eqn 2.2
        # putting genes on same scale
		nData <- data[,phenoRand==0]
        medG <- apply(data, 1, median)
        madGN <- apply(nData, 1, mad)
        
        for (i in 1:nG) {
            temp <- (data[i,] - medG[i])/madGN[i]
            dataC[i,] <- temp
        }

        # calc outlier stat on scaled data
        tData <- dataC[,phenoRand==1]
		nData <- dataC[,phenoRand==0]

        iqrG <- apply(dataC, 1, IQR)
        if (tail == 'right') {
            quantG <- apply(dataC, 1, quantile, probs = 0.75)
        } else if (tail == 'left') {
            quantG <- apply(dataC, 1, quantile, probs = 0.25)
        }
        
        if (tail == 'right') {
            for (i in 1:nG) {
                Ind <- tData[i,] > quantG[i] + iqrG[i]
                copaValue[i] <- sum(tData[i,] * Ind)
            }
        } else if (tail == 'left') {
            for (i in 1:nG) {
                Ind <- tData[i,] < quantG[i] - iqrG[i]
                copaValue[i] <- sum(tData[i,] * Ind)
            }
        }
        copaPerm[,j] <- copaValue
	}
	nPerm <- dim(copaPerm)[1] * dim(copaPerm)[2]

    # Now for real data
    # generate Tibshirani eqn 2.2
    # putting genes on same scale
	medG <- apply(data, 1, median)
	madGN <- apply(nData, 1, mad)
    
    for (i in 1:nG) {
		temp <- (data[i,] - medG[i])/madGN[i]
        data[i,] <- temp
	}
    
    # generate stat on scaled data
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
            copaValue[i] <- sum(tData[i,] * Ind)
        }
    } else if (tail == 'left') {
        for (i in 1:nG) {
            Ind <- tData[i,] < quantG[i] - iqrG[i]
            copaValue[i] <- sum(tData[i,] * Ind)
        }
    }

	# generate comparison of COPA values to permutations
	if(permType == 'array') {
		nPerm <- dim(copaPerm)[1] * dim(copaPerm)[2]
        if (tail == 'right') {
            for (i in 1:nG) {
                copaStat[i] <- (sum(copaPerm > copaValue[i])+1)/(nPerm+1)
            }
        } else if (tail == 'left') {
            for (i in 1:nG) {
                copaStat[i] <- (sum(copaPerm < copaValue[i])+1)/(nPerm+1)
            }
		}
	} else if (permType == 'gene') {
		nPerm <- dim(copaPerm)[2]
        if (tail == 'right') {
            for (i in 1:nG) {
                copaStat[i] <- (sum(copaPerm[i,] > copaValue[i])+1)/(nPerm+1)
            }
        } else if (tail == 'left') {
            for (i in 1:nG) {
                copaStat[i] <- (sum(copaPerm < copaValue[i])+1)/(nPerm+1)
            }
		}
    } else {
        print("Unknown Permutation Comparison Type")
    }

	names(copaStat) <- rownames(data)

	return(copaStat)
}		
	

