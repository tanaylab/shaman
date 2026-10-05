# Inputs that used to crash the session or never finish. The shuffler runs in a forked child, so
# that a crash fails the test instead of killing the test process.
library(misha)

in_child <- function(expr) parallel::mccollect(parallel::mcparallel(expr))[[1]]

# chr1: 4,000 contacts (8,000 in the shuffler), shuffled because the chromosome has hg38's chr1
# size; chr2: 3,000 contacts (6,000 in both orientations) for local mode
db <- make_track(tempfile(), c(chr1 = 248956422, chr2 = 4e6), c(chr1 = 4000, chr2 = 3000))

test_that("a chromosome with fewer than 10,000 contacts is shuffled instead of dividing by zero", {
    gsetroot(db)
    work_dir <- tempfile()
    dir.create(work_dir)
    on.exit(unlink(work_dir, recursive = TRUE))
    ret <- in_child(shaman_shuffle_hic_mat_for_track(db, "hic_obs", work_dir, "chr1", 0, 248956422, 0, 248956422,
        shuffle = 2, grid_step_iter = 1
    ))
    expect_false(is.null(ret))
    obs <- gextract("hic_obs", gintervals.2d(1), band = c(-248956422, -1023), colnames = "contacts")
    s <- as.data.frame(data.table::fread(file.path(work_dir, "hic_obs_chr1_0_0.shuffled")))
    # each contact in both orientations, every contact end kept
    expect_equal(sort(s$start1), sort(rep(c(obs$start1, obs$start2), 2)))
})

test_that("local mode on fewer than 10,000 contacts runs instead of dividing by zero", {
    gsetroot(db)
    work_dir <- tempfile()
    dir.create(work_dir)
    on.exit(unlink(work_dir, recursive = TRUE))
    res <- in_child(shaman_shuffle_and_score_hic_mat("hic_obs", gintervals.2d(2, 0, 4e6, 2, 0, 4e6), work_dir,
        shuffle = 2, grid_step_iter = 1
    ))
    expect_false(is.null(res))
    expect_equal(nrow(res$points), 3000)
    expect_true(all(res$points$score >= -100 & res$points$score <= 100))
})

test_that("too few expected contacts for k_exp give no score instead of an error", {
    # 1,200 observed and 150 expected contacts (both orientations)
    db <- make_track(tempfile(), c(chr1 = 5e6), c(chr1 = 600))
    make_track(db, c(chr1 = 5e6), c(chr1 = 75), "hic_exp")
    all <- gintervals.2d(1, 0, 5e6, 1, 0, 5e6)
    expect_null(shaman_score_hic_mat("hic_obs", "hic_exp", all, all, k = 100, k_exp = 200))
    expect_false(is.null(shaman_score_hic_mat("hic_obs", "hic_exp", all, all, k = 100, k_exp = 150)))
})

test_that("the score step stops after 3 rounds and names the matrices that got no score", {
    db <- make_track(tempfile(), c(chr1 = 5e6), c(chr1 = 1000))
    make_track(db, c(chr1 = 5e6), c(chr1 = 1000), "hic_exp")
    calls <- tempfile()
    on.exit(unlink(calls))
    # every matrix fails without writing its score file, as when knn does not complete
    local_mocked_bindings(shaman_score_hic_mat_for_track = function(track_db, work_dir, obs_track_nms, exp_track_nms,
                                                                    points_track_nms, chrom, start1, end1, start2, end2, ...) {
        cat(chrom, start1, start2, "\n", file = calls, append = TRUE)
        # without a cap the rounds never end
        if (length(readLines(calls)) > 100) stop("still retrying")
        -1
    })
    old_opts <- options(shaman.mc_support = 1, shaman.sge_support = 0)
    on.exit(options(old_opts), add = TRUE)
    work_dir <- tempfile()
    dir.create(work_dir)
    on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
    expect_error(
        shaman_score_hic_track(db, work_dir, "hic_score_test", "hic_obs", "hic_exp", near_cis = 2e6, max_jobs = 2),
        "6 matrices have no score after 3 rounds: chr1:0-2000000 x 0-2000000, "
    )
    expect_equal(as.vector(table(readLines(calls))), rep(3, 6))
    expect_false(gtrack.exists("hic_score_test"))
})
