# Regenerates the usage article, vignettes/shaman-package.Rmd, and its figures (vignettes/shaman-package-*.png)
# by running vignettes/shaman-package.Rmd.orig on the full example database, then checks the results (below).
# Edit the .orig, not the .Rmd.
#
# Run from the package root, with this version of shaman installed:
#   Rscript vignettes/precompute.R
# About 13 minutes on 4 cores, 3GB of memory and 2.5GB of disk (measured on 2026-10-01).
#
# CI runs it in .github/workflows/full-example.yaml, called from pkgdown.yaml on pushes to master (the site
# deployed from master shows the regenerated article) and cran-readiness, and also weekly and by hand. The package
# (CRAN, R CMD check) ships the committed .Rmd and figures, which have no code to run: commit them again when they
# change.
#
# The database is the one of shaman_get_test_track_db(full = TRUE), in tools::R_user_dir("shaman", "cache"). If
# R_USER_CACHE_DIR is not set, a new temporary directory is used, so a local run downloads the database and computes
# everything. The tracks the article computes (hic_obs_shuffle, hic_score_new) are reused if they are already in the
# database: CI restores them from its cache when src/, R/, the article and the toolchain are unchanged.

if (!file.exists("vignettes/shaman-package.Rmd.orig")) {
    stop("run from the package root: Rscript vignettes/precompute.R")
}
if (Sys.getenv("R_USER_CACHE_DIR") == "") {
    Sys.setenv(R_USER_CACHE_DIR = tempfile("shaman_cache_"))
}
library(misha)
library(shaman)
track_db <- shaman_get_test_track_db(full = TRUE)
# rescan, for tracks put into the database by a cache restore
gsetroot(track_db, rescan = TRUE)

owd <- setwd("vignettes")
unlink(Sys.glob("shaman-package-*.png"))
knitr::knit("shaman-package.Rmd.orig", "shaman-package.Rmd", envir = new.env())
setwd(owd)

gsetroot(track_db)
failed <- character(0)
check <- function(ok, what) {
    message(if (ok) "ok: " else "FAILED: ", what)
    if (!ok) failed <<- c(failed, what)
}

# 1. Scoring hic_obs against hic_exp reproduces the database's hic_score track (computed in 2017), on its 464,708
#    points in chr2:175e06-178e06: r = 0.9996 on 2026-10-01 (1.3% of the points differ by more than 0.1).
reg <- gintervals.2d(2, 175e06, 178e06, 2, 175e06, 178e06)
old <- gextract("hic_score", reg, colnames = "old")
old <- old[old$start1 < old$start2, c("start1", "start2", "old")]
cor_with_hic_score <- function(points) {
    m <- merge(old, points[points$start1 < points$start2, c("start1", "start2", "score")])
    if (nrow(m) != nrow(old)) {
        stop(nrow(old) - nrow(m), " of the ", nrow(old), " points of hic_score have no new score")
    }
    cor(m$old, m$score)
}
r_exp <- cor_with_hic_score(shaman_score_hic_mat("hic_obs", "hic_exp", reg,
    gintervals.2d(2, 173e06, 182e06, 2, 173e06, 182e06))$points)
check(r_exp >= 0.999, sprintf("scores against hic_exp vs hic_score: r = %.4f (>= 0.999)", r_exp))
# Reported only: the article's hic_score_new (expected shuffled from this database) vs hic_score, r = 0.378 on
# 2026-10-01. Low because the database is a cut-out (see the note in the article).
message(sprintf("hic_score_new vs hic_score: r = %.4f (reported, no threshold)",
    cor_with_hic_score(gextract("hic_score_new", reg, colnames = "score"))))

# 2. The shuffle is reproducible with a seed: chr2 shuffled twice more, with the seed recorded in hic_obs_shuffle and
#    the defaults of shaman_shuffle_hic_track() (which the article uses), gives the same contacts and counts both
#    times; and the same as hic_obs_shuffle itself, unless that track came from the CI cache
#    (SHAMAN_TRACKS_CACHED=true), which may have been computed on another runner image (compiler, libm).
seeds <- do.call(rbind, strsplit(strsplit(gtrack.attr.get("hic_obs_shuffle", "seed"), " ")[[1]], ":"))
seeds <- stats::setNames(suppressWarnings(as.integer(seeds[, 2])), seeds[, 1])
chroms <- gintervals.all()
chr2_end <- chroms$end[chroms$chrom == "chr2"]
shuffle_chr2 <- function() {
    work_dir <- tempfile("shaman_rerun_")
    dir.create(work_dir)
    on.exit(unlink(work_dir, recursive = TRUE))
    seed <- shaman_shuffle_hic_mat_for_track(track_db, "hic_obs", work_dir, "chr2", 0, chr2_end, 0, chr2_end,
        seed = seeds[["chr2"]], sort_uniq = TRUE
    )
    stopifnot(isTRUE(as.integer(seed) == seeds[["chr2"]]))
    x <- as.data.frame(data.table::fread(file.path(work_dir, "hic_obs_chr2_0_0.shuffled.uniq")))
    x[order(x$start1, x$start2), c("start1", "start2", "obs")]
}
same_contacts <- function(x, y) {
    nrow(x) == nrow(y) && all(x$start1 == y$start1) && all(x$start2 == y$start2) && all(x$obs == y$obs)
}
rerun1 <- shuffle_chr2()
rerun2 <- shuffle_chr2()
check(same_contacts(rerun1, rerun2),
    sprintf("chr2 shuffled twice with seed %d: the same %d contacts both times", seeds[["chr2"]], nrow(rerun1)))
if (Sys.getenv("SHAMAN_TRACKS_CACHED") != "true") {
    shuffled <- gextract("hic_obs_shuffle", gintervals.2d(2, 0, chr2_end, 2, 0, chr2_end), colnames = "obs")
    shuffled <- shuffled[order(shuffled$start1, shuffled$start2), c("start1", "start2", "obs")]
    check(same_contacts(rerun1, shuffled),
        sprintf("chr2 shuffled again with seed %d: the same %d contacts as hic_obs_shuffle", seeds[["chr2"]], nrow(shuffled)))
} else {
    message("hic_obs_shuffle came from the CI cache: not compared with the reruns")
}
rm(rerun1, rerun2)

# 3. The shuffle keeps the marginal coverage: on every chromosome, the contacts of hic_obs_shuffle at each position
#    are twice the contact ends of the observed contacts it shuffled (each end kept, each contact stored in both
#    orientations), or equal to them on a chromosome left unshuffled (seed NA: too few contacts).
max_dist <- max(chroms$end)
coverage_kept <- vapply(seq_len(nrow(chroms)), function(i) {
    chrom <- gintervals.2d(chroms$chrom[i], 0, chroms$end[i], chroms$chrom[i], 0, chroms$end[i])
    o <- gextract("hic_obs", chrom, band = c(-max_dist, -1023), colnames = "n")
    e <- gextract("hic_obs_shuffle", chrom, colnames = "n")
    if (is.null(o) || is.null(e)) {
        return(is.null(o) == is.null(e))
    }
    o <- data.table::data.table(pos = c(o$start1, o$start2), n = c(o$n, o$n))[, list(n = sum(n)), by = "pos"]
    e <- data.table::data.table(pos = e$start1, n = e$n)[, list(n = sum(n)), by = "pos"]
    data.table::setorder(o, pos)
    data.table::setorder(e, pos)
    times <- if (is.na(seeds[[chroms$chrom[i]]])) 1 else 2
    nrow(o) == nrow(e) && all(o$pos == e$pos) && all(times * o$n == e$n)
}, TRUE)
check(all(coverage_kept), sprintf("marginal coverage of hic_obs_shuffle kept on %d of %d chromosomes",
    sum(coverage_kept), length(coverage_kept)))

if (length(failed) > 0) {
    stop(length(failed), " check(s) failed: ", paste(failed, collapse = "; "))
}
