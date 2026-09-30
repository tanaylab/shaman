#' shaman.
#'
#' @import misha
#' @name shaman
#' @docType package
#' @useDynLib shaman
#' @importFrom Rcpp sourceCpp
#' @importFrom graphics abline image par
#' @importFrom grDevices colorRampPalette dev.off png rgb
#' @importFrom stats setNames
#' @importFrom utils read.table
NULL

# column names used in ggplot2::aes() and plyr calls
utils::globalVariables(c("chrom1", "score", "start", "start1", "start2", "value"))
