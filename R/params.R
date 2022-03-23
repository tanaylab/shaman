#' Read params from files
init_params <- function(fn, prev_params = NULL) {
    t <- read.table(fn, sep = "=", fill = T, strip.white = T, stringsAsFactors = FALSE, quote = "")
    params <- as.character(t[, 2])
    names(params) <- t[, 1]

    if (any(duplicated(names(params)))) {
        warning("duplicated params in conf file! check ", fn)
        return(NA)
    }

    if (!is.null(prev_params)) {
        if (any(names(prev_params)) %in% names(params)) {
            stop(sprintf("duplicaterd parameters: %s", names(prev_parames)[names(prev_params) %in% names(params)]))
        }
    }
    params <- as.list(params)
    # interpret some values as R exprssions
    expr_idx <- grep("\\@R$", names(params))
    # Strip names from R
    params[expr_idx] <- lapply(params[expr_idx], function(x) eval(parse(text = x)))
    names(params)[expr_idx] <- sub("\\@R$", "", names(params[expr_idx]))

    return(params)
}
#' Get params from saved var
get_param <- function(nm, params) {
    if (nm %in% names(params)) {
        return(params[nm])
    } else {
        message("missing params ", nm)
        return(NA)
    }
}
#' Get params from saved var
get_param_list <- function(nm, params) {
    if (nm %in% names(params)) {
        return(strsplit(as.character(params[nm]), ","))
    } else {
        assign("shaman_miss_conf_err", TRUE, envir = .GlobalEnv)
        message("missing params ", nm)
        return(NA)
    }
}

############################################
#' returns test misha db
#'
#' \code{shaman_get_test_track_db}
#' Sets the current misha DB to the example database provided with shaman.
#' In the example misha database provided in this package we have created a low-footprint
#' matrix to examplify the shaman workflow. We included 4.6 million contacts from
#' ELA K562 dataset covering the hoxd locus (chr2:175e06-178e06) and convergent CTCF regions.
#' Processing the complete matrix from this study requires downloading the full contact list
#' and regenerating the reshuffled matrix.
#' @export


shaman_get_test_track_db <- function() {
    return(sprintf("%s/test", system.file("trackdb", package = "shaman")))
}
