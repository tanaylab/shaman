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
