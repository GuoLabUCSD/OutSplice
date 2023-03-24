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
#  v 1.1 - Joseph Bendik, UCSD, 07 Mar 2023 changed sapply to vapply, changed 1:... to seq_len(n) or seq_along(x)


outCallRank <- function(dataList, thres = 0.05, tail = "right", corr = FALSE, offsets = NULL, names = NULL) {
  if (is.null(names)) {
    names <- vector(length = length(dataList), mode = "character")
    for (d in seq_along(dataList)) {
      names[d] <- paste("Data", d)
    }
  }
  temp <- dataList[[1]]
  temp <- temp[[1]]
  nG <- dim(temp)[1]
  outList <- list()

  if (!corr) {
    offsets <- rep(0.0, dim(temp)[2])
  } else if (is.null(offsets)) {
    print("No Offsets Set with Correction Requested")
    return()
  }

  outP <- matrix(nrow = nG, ncol = 2)
  outCount <- rep(0, nG)

  for (d in seq_along(dataList)) {
    data <- dataList[[d]]
    phenotype <- data[[2]]
    data <- data[[1]]
    nS <- length(phenotype)
    thisTail <- tail[d]

    adjust <- offsets[d]
    nData <- data[, phenotype == 0]
    tData <- data[, phenotype == 1]
    nT <- dim(tData)[2]

    # generate empicical pValues as the
    # number of sum(normals{<,>}tumor)/nN
    empirP <- matrix(nrow = nG, ncol = nT)
    if (thisTail == "right") {
      for (i in seq_len(nG)) {
        tumor <- tData[i, ]
        baseline <- nData[i, ]
        result <- vapply(seq_along(tumor), function(j) {
          sum((baseline + adjust) > tumor[j])
        }, integer(1))
        empirP[i, ] <- result / length(baseline)
      }
    } else if (thisTail == "left") {
      for (i in seq_len(nG)) {
        tumor <- tData[i, ]
        baseline <- nData[i, ]
        result <- vapply(seq_along(tumor), function(j) {
          sum((baseline - adjust) < tumor[j])
        }, integer(1))
        empirP[i, ] <- result / length(baseline)
      }
    }
    empirP <- empirP < thres
    rownames(empirP) <- rownames(data)
    colnames(empirP) <- colnames(tData)
    outList[[d]] <- empirP
  }
  names(outList) <- names
  return(outList)
}
