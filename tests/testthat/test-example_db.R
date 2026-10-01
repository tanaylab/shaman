# Shuffle and score on the small example database (shaman_get_test_track_db()): a few seconds, at most
# 2 processes. The shuffler writes its progress to R's error stream; quiet() hides it.
library(misha)

db <- shaman_get_test_track_db()
quiet <- function(expr) {
    utils::capture.output(res <- suppressMessages(expr), type = "message")
    res
}

test_that("the example database has the hoxd contacts", {
    gsetroot(db)
    expect_true(all(c("hic_obs", "hic_exp", "hic_score") %in% gtrack.ls()))
    obs <- gextract("hic_obs", gintervals.2d(2, 176.5e6, 177e6, 2, 176.5e6, 177e6), colnames = "contacts")
    expect_equal(sum(obs$start1 < obs$start2), 46735)
    # both orientations of each contact
    expect_equal(sum(obs$start1 > obs$start2), 46735)
})

test_that("a shuffle keeps the marginal coverage, is reproducible with a seed and restores the options", {
    gsetroot(db)
    max_dist <- max(gintervals.all()$end)
    # the contacts shaman_shuffle_hic_mat_for_track() shuffles, each repeated by its count
    obs <- gextract("hic_obs", gintervals.2d.all(), band = c(-max_dist, -1023), colnames = "contacts")
    x <- rep(obs$start1, obs$contacts)
    y <- rep(obs$start2, obs$contacts)
    shuffled <- function(seed) {
        work_dir <- tempfile()
        dir.create(work_dir)
        on.exit(unlink(work_dir, recursive = TRUE))
        ret <- quiet(shaman_shuffle_hic_mat_for_track(db, "hic_obs", work_dir, "chr2", 0, max_dist, 0, max_dist,
            shuffle = 1, grid_step_iter = 1, seed = seed
        ))
        expect_equal(ret, seed)
        as.data.frame(data.table::fread(file.path(work_dir, "hic_obs_chr2_0_0.shuffled")))
    }
    old_opts <- options(gmultitasking = TRUE, gmax.data.size = 12345678)
    on.exit(options(old_opts))
    s1 <- shuffled(5)
    expect_identical(getOption("gmultitasking"), TRUE)
    expect_identical(getOption("gmax.data.size"), 12345678)
    # the shuffled file has each contact in both orientations, and swapping partners keeps every end
    expect_equal(sort(s1$start1), sort(c(x, y, x, y)))
    expect_equal(sort(s1$start2), sort(c(x, y, x, y)))
    # but the contacts moved
    expect_lt(mean(paste(s1$start1, s1$start2) %in% paste(c(x, y), c(y, x))), 0.9)
    expect_identical(shuffled(5), s1)
})

test_that("scoring gives a score in [-100, 100] for each point in the focus interval", {
    gsetroot(db)
    focus <- gintervals.2d(2, 176.7e6, 176.8e6, 2, 176.7e6, 176.8e6)
    regional <- gintervals.2d(2, 176.5e6, 177e6, 2, 176.5e6, 177e6)
    res <- quiet(shaman_score_hic_mat("hic_obs", "hic_exp", focus, regional, k = 20))
    p <- res$points
    expect_true(all(c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "score") %in% names(p)))
    expect_gt(nrow(p), 1000)
    expect_false(anyNA(p$score))
    expect_true(all(p$score >= -100 & p$score <= 100))
    expect_true(all(p$start1 >= 176.7e6 & p$start1 < 176.8e6 & p$start2 >= 176.7e6 & p$start2 < 176.8e6))
})

test_that("multi-core mode leaves the foreach backend as it was", {
    skip_if_not_installed("foreach")
    foreach::registerDoSEQ()
    backend <- foreach::getDoParName()
    old_opts <- options(shaman.mc_support = 1)
    on.exit(options(old_opts))
    work_dir <- tempfile()
    dir.create(work_dir)
    on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
    quiet(shaman_score_hic_track(db, work_dir, "hic_score_test", "hic_obs", "hic_exp", near_cis = 1e9, k = 20, max_jobs = 2))
    quiet(shaman_shuffle_hic_track(db, "hic_obs", work_dir, "hic_obs_shuffle_test",
        max_jobs = 2, shuffle = 1, grid_step_iter = 1, seed = 1
    ))
    gdb.reload()
    expect_true(all(c("hic_score_test", "hic_obs_shuffle_test") %in% gtrack.ls()))
    expect_identical(gtrack.attr.get("hic_obs_shuffle_test", "seed"), "chr2:1")
    expect_identical(foreach::getDoParName(), backend)
})
