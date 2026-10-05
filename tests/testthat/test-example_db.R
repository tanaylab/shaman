# Shuffle and score on the small example database (shaman_get_test_track_db()): a few seconds, at most
# 2 processes. The shuffler writes its progress to R's error stream; quiet() hides it.
library(misha)

db <- shaman_get_test_track_db()
quiet <- function(expr) {
    utils::capture.output(res <- suppressMessages(expr), type = "message")
    res
}
digest_of <- function(d) {
    fn <- tempfile()
    on.exit(unlink(fn))
    data.table::fwrite(d, fn, sep = "\t")
    unname(tools::md5sum(fn))
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
    # another seed, another shuffle (a seed that set only the return value would pass the line above)
    expect_false(identical(shuffled(6), s1))
    # the exact shuffle (made with g++ 13 on glibc 2.28): a change in the shuffler's output
    # shows here; another libm or architecture could change it too, so not on CRAN
    skip_on_cran()
    skip_if_not(Sys.info()[["sysname"]] == "Linux" && R.version$arch == "x86_64")
    expect_identical(digest_of(s1), "cd79437a5026d823dc76287418dddaa8")
})

test_that("the shuffler stops with an error when it cannot write its output", {
    skip_if_not(file.exists("/dev/full"))
    out <- file.path(tempfile(), "x.shuffled")
    dir.create(dirname(out))
    on.exit(unlink(dirname(out), recursive = TRUE))
    # a full disk: writes to /dev/full fail with ENOSPC
    file.symlink("/dev/full", out)
    x <- 1e6L + seq_len(5000) * 100L
    shuffle <- function(fn) {
        quiet(shaman_hic_matrix_shuffler_cpp(rbind(x, x + 5000L), fn, 0, 1, 0.5, 2, 1, 5, 0.25, 2e7, 1024, 1, 5e5, 1e6, 5e5, 1, 0, 1, 1))
    }
    expect_error(shuffle(out), "could not write .*x.shuffled: No space left on device")
    # the partial file (here the link) is removed, so a rerun does not take it as done
    expect_false(basename(out) %in% list.files(dirname(out)))
    expect_error(shuffle(file.path(dirname(out), "no_dir", "y.shuffled")), "could not open output file .*y.shuffled")
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
