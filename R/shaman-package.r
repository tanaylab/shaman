#' @keywords internal
#' @import misha
#' @useDynLib shaman
#' @importFrom Rcpp sourceCpp
#' @importFrom graphics abline image par
#' @importFrom grDevices colorRampPalette dev.off png rgb
#' @importFrom stats setNames
#' @importFrom utils read.table
"_PACKAGE"

# column names used in ggplot2::aes() and plyr calls
utils::globalVariables(c("score", "start1", "start2", "value"))
