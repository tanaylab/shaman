#' Read params from files
#' @param fn Parameter file, one \code{name=value} per line.
#' @param prev_params Previously read parameters; names found in both are an error.
#' @return A named list of the parameters.
init_params <- function(fn, prev_params = NULL) {
    t <- read.table(fn, sep = "=", fill = TRUE, strip.white = TRUE, stringsAsFactors = FALSE, quote = "")
    params <- as.character(t[, 2])
    names(params) <- t[, 1]

    if (any(duplicated(names(params)))) {
        warning("duplicated params in conf file! check ", fn)
        return(NA)
    }

    if (!is.null(prev_params)) {
        if (any(names(prev_params) %in% names(params))) {
            stop(sprintf("duplicated parameters: %s", paste(names(prev_params)[names(prev_params) %in% names(params)], collapse = ", ")))
        }
    }
    params <- as.list(params)
    # interpret some values as R expressions
    expr_idx <- grep("\\@R$", names(params))
    # Strip names from R
    params[expr_idx] <- lapply(params[expr_idx], function(x) eval(parse(text = x)))
    names(params)[expr_idx] <- sub("\\@R$", "", names(params[expr_idx]))

    return(params)
}
#' Get params from saved var
#' @param nm Parameter name.
#' @param params Named list of parameters.
#' @return The parameter, or NA if it is missing.
get_param <- function(nm, params) {
    if (nm %in% names(params)) {
        return(params[nm])
    } else {
        message("missing params ", nm)
        return(NA)
    }
}
#' Get params from saved var
#' @param nm Parameter name.
#' @param params Named list of parameters.
#' @return The comma-separated parameter split into a list, or NA if it is missing.
get_param_list <- function(nm, params) {
    if (nm %in% names(params)) {
        return(strsplit(as.character(params[nm]), ","))
    } else {
        message("missing params ", nm)
        return(NA)
    }
}

############################################
#' Returns the example misha database
#'
#' \code{shaman_get_test_track_db}
#' Returns the path of an example misha database with Hi-C contacts from the ELA K562 dataset.
#'
#' By default this is a small database that is built in the session's temporary directory on first
#' use; later calls in the session reuse it. It has one chromosome (hg19 chr2) and the observed
#' (hic_obs), expected (hic_exp) and score (hic_score) tracks of chr2:176.5e06-177e06, around the
#' hoxd locus (46,735 observed contacts, each stored in both orientations).
#'
#' With \code{full = TRUE} it returns the full example database of earlier shaman versions: 4.6 million
#' contacts from the ELA K562 dataset covering the hoxd locus (chr2:175e06-178e06) and convergent CTCF
#' regions, with scores for chr2:175e06-178e06. It is downloaded (about 100MB, from the lab's public S3 bucket) and extracted into the
#' user cache directory (\code{tools::R_user_dir("shaman", "cache")}, about 480MB) on first use; later
#' calls reuse it.
#' Processing the complete matrix from this study requires downloading the full contact list
#' and regenerating the reshuffled matrix.
#' @param full Whether to return the full example database (downloaded on first use) instead of the small one.
#' @return The path of the database.
#' @examples
#' library(misha)
#' gsetroot(shaman_get_test_track_db())
#' gtrack.ls()
#' @export
shaman_get_test_track_db <- function(full = FALSE) {
    if (full) {
        return(.shaman_get_full_test_track_db())
    }
    track_db <- file.path(tempdir(), "shaman_test_db")
    if (dir.exists(track_db)) {
        return(track_db)
    }
    # built next to its final place and moved there, so a failed build leaves no partial db
    tmp_db <- tempfile("shaman_test_db_")
    on.exit(unlink(tmp_db, recursive = TRUE), add = TRUE)
    dir.create(file.path(tmp_db, "tracks"), recursive = TRUE)
    # hg19 chr2, with an empty seq/chr2.seq (see .shaman_get_full_test_track_db)
    writeLines("2\t243199373", file.path(tmp_db, "chrom_sizes.txt"))
    dir.create(file.path(tmp_db, "seq"))
    file.create(file.path(tmp_db, "seq", "chr2.seq"))
    # the contacts (start1 < start2) of the full example db in chr2:176.5e06-177e06
    contacts <- utils::read.delim(system.file("extdata", "hoxd.tsv.xz", package = "shaman"))
    # gdb.info() needs misha >= 5.3; with older misha the example db stays the current root
    old_root <- tryCatch(gdb.info()$path, error = function(e) NULL)
    if (!is.null(old_root)) {
        on.exit(gsetroot(old_root), add = TRUE)
    }
    gsetroot(tmp_db)
    for (track in unique(contacts$track)) {
        x <- contacts[contacts$track == track, ]
        fn <- tempfile(fileext = ".txt")
        # both orientations, as in the full example db
        data.table::fwrite(data.frame(
            chrom1 = "chr2", start1 = c(x$start1, x$start2), end1 = c(x$start1, x$start2) + 1L,
            chrom2 = "chr2", start2 = c(x$start2, x$start1), end2 = c(x$start2, x$start1) + 1L,
            value = c(x$value, x$value)
        ), fn, sep = "\t")
        gtrack.2d.import(track, paste(track, "of the shaman example database, chr2:176.5e06-177e06"), fn)
        unlink(fn)
    }
    if (!file.rename(tmp_db, track_db) && !dir.exists(track_db)) {
        stop("failed to create the shaman example database in ", track_db)
    }
    return(track_db)
}

.shaman_get_full_test_track_db <- function() {
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
            url <- "https://misha-genome.s3.eu-west-1.amazonaws.com/shaman/trackdb.tar.gz"
            tryCatch(utils::download.file(url, tarball, mode = "wb"), error = function(e) {
                stop("could not download the shaman example database from ", url, ": ", conditionMessage(e), call. = FALSE)
            })
            if (unname(tools::md5sum(tarball)) != "88552541e7bf346ef50187f1e42fc25b") {
                stop("the shaman example database downloaded from ", url, " is not the expected file (md5 mismatch)")
            }
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
