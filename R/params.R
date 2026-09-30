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
#' Returns the path of the example misha database provided with shaman.
#' On first use the database is extracted into the user cache directory
#' (\code{tools::R_user_dir("shaman", "cache")}); later calls reuse it. If the installed
#' package holds only a git-lfs pointer instead of the tarball, the tarball is downloaded
#' from GitHub first (about 100MB).
#' In the example misha database provided in this package we have created a low-footprint
#' matrix to examplify the shaman workflow. We included 4.6 million contacts from
#' ELA K562 dataset covering the hoxd locus (chr2:175e06-178e06) and convergent CTCF regions.
#' Processing the complete matrix from this study requires downloading the full contact list
#' and regenerating the reshuffled matrix.
#' @export


shaman_get_test_track_db <- function() {
    cache_dir <- tools::R_user_dir("shaman", "cache")
    track_db <- file.path(cache_dir, "trackdb", "test")
    if (!dir.exists(track_db)) {
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
        tarball <- system.file("trackdb.tar.gz", package = "shaman")
        # an install from a checkout without git-lfs has a text pointer instead of the tarball
        if (tarball == "" || !identical(readBin(tarball, "raw", 2), as.raw(c(0x1f, 0x8b)))) {
            tarball <- tempfile(tmpdir = cache_dir, fileext = ".tar.gz")
            on.exit(unlink(tarball), add = TRUE)
            old_opts <- options(timeout = max(600, getOption("timeout")))
            on.exit(options(old_opts), add = TRUE)
            message("downloading the shaman example database (about 100MB)")
            utils::download.file("https://media.githubusercontent.com/media/tanaylab/shaman/master/inst/trackdb.tar.gz",
                tarball,
                mode = "wb"
            )
        }
        message("extracting the shaman example database to ", cache_dir)
        # extract next to the cache and move into place, so an interrupted run leaves no partial db
        tmp_dir <- tempfile(tmpdir = cache_dir)
        on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
        if (utils::untar(tarball, exdir = tmp_dir) != 0) {
            stop("failed to extract the shaman example database from ", tarball)
        }
        # The example db has no sequence. Current misha (5.11) refuses a root without a seq directory,
        # and it adds the "chr" prefix to the names in chrom_sizes.txt ("1", "2", ...) only when
        # seq/chr*.seq files exist; without it the names would not match the track files (chr1-chr1, ...).
        # Empty placeholder files are enough. Older misha (4.x) always adds the prefix and ignores them.
        db_dir <- file.path(tmp_dir, "trackdb", "test")
        dir.create(file.path(db_dir, "seq"))
        chroms <- sub("\t.*", "", readLines(file.path(db_dir, "chrom_sizes.txt")))
        file.create(file.path(db_dir, "seq", paste0("chr", chroms, ".seq")))
        # another process may have put the database in place in the meantime
        moved <- suppressWarnings(file.rename(file.path(tmp_dir, "trackdb"), file.path(cache_dir, "trackdb")))
        if (!moved && !dir.exists(track_db)) {
            stop("failed to move the shaman example database into ", cache_dir)
        }
    }
    return(track_db)
}
