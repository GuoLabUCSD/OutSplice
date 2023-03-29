#' OutSplice package overview
#'
#' An easy to use tool that can compare splicing events in tumor and normal tissue
#' samples using either a user generated matrix, or data from The Cancer Genome Atlas (TCGA).
#' This package generates a matrix of splicing outliers that are significantly
#' over or underexpressed in tumors samples compared to normal denoted by chromosome location.
#' The package also will calculate the splicing burden in each tumor, characterize
#' the types of splicing events that occur, and allows the user to create waterfall
#' plots of event expression.
#'
#' Below are the available functions provided by OutSplice:
#' \itemize{
#'     \item \code{\link{outspliceAnalysis}}
#'     \item \code{\link{outspliceTCGA}}
#'     \item \code{\link{plotJunctionData}}
#' }
#'
#' @name OutSplice
#'
#' @import GenomicRanges
#' @importFrom MatrixGenerics rowSums colSums rowMeans colMeans
#' @importFrom GenomicFeatures genes exons
#' @importFrom IRanges IRanges overlapsAny
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom Repitools annoGR2DF
#' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene TxDb.Hsapiens.UCSC.hg38.knownGene
#' @importFrom S4Vectors Rle subjectHits queryHits
#' @importFrom AnnotationDbi select
NULL
