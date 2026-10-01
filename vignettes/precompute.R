# Regenerates the usage article, vignettes/shaman-package.Rmd, and its figures (vignettes/shaman-package-*.png)
# by running vignettes/shaman-package.Rmd.orig on the full example database, then checks the new score track
# against the hic_score track shipped in the database. Edit the .orig, not the .Rmd.
#
# Run from the package root, with this version of shaman installed:
#   Rscript vignettes/precompute.R
# About 10 minutes on 4 cores, 2.5GB of memory and 2.5GB of disk (measured on 2026-10-01).
#
# CI runs it in .github/workflows/full-example.yaml, called from pkgdown.yaml on every push to master (the site
# deployed from master shows the regenerated article), and also weekly and by hand. The package (CRAN, R CMD check)
# ships the committed .Rmd and figures, which have no code to run: commit them again when they change.
#
# The database is the one of shaman_get_test_track_db(full = TRUE), in tools::R_user_dir("shaman", "cache"). If
# R_USER_CACHE_DIR is not set, a new temporary directory is used, so a local run downloads the database and computes
# everything. The tracks the article computes (hic_obs_shuffle, hic_score_new) are reused if they are already in the
# database: CI restores them from its cache when src/, R/ and the article are unchanged.

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

# Checks, on the 464,708 observed points of chr2:175e06-178e06 that have a score in the database's hic_score track
# (computed in 2017):
# 1. the article's hic_score_new (expected shuffled with seed 1) against hic_score. On 2026-10-01: r = 0.378, and
#    0.377 with seed 2 (the scores of the two seeds correlate at 0.979). It is low because the database is a
#    cut-out (see the note in the article): its hic_exp, which hic_score was computed with, is not a shuffle of its
#    hic_obs. The threshold, 0.35, allows for a drift ~30 times the difference between the two seeds.
# 2. today's scoring of hic_obs against hic_exp, against hic_score: r = 0.9996 on 2026-10-01 (1.3% of the points
#    differ by more than 0.1). This checks the scoring alone, with the threshold 0.999.
gsetroot(track_db)
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
r_new <- cor_with_hic_score(gextract("hic_score_new", reg, colnames = "score"))
r_exp <- cor_with_hic_score(shaman_score_hic_mat("hic_obs", "hic_exp", reg,
    gintervals.2d(2, 173e06, 182e06, 2, 173e06, 182e06))$points)
message(sprintf("correlation with hic_score: hic_score_new %.4f (threshold 0.35), scored against hic_exp %.4f (threshold 0.999)",
    r_new, r_exp))
if (r_new < 0.35 || r_exp < 0.999) {
    stop("the scores no longer agree with the database's hic_score track")
}
