# The shortcuts of the score step give the same results as the direct path: matrices of one row
# scored in one call (one extraction of each track for the row) and one by one, and the focus points
# taken from the regional extraction and from their own extraction.
library(misha)

db <- make_track(tempfile(), c(chr1 = 5e6), c(chr1 = 30000))
make_track(db, c(chr1 = 5e6), c(chr1 = 60000), "hic_exp", seed = 2)

test_that("matrices of one row scored in one call give the same score files as one by one", {
    gsetroot(db)
    score_dir <- function() {
        d <- tempfile()
        dir.create(d)
        d
    }
    together <- score_dir()
    alone <- score_dir()
    on.exit(unlink(c(together, alone), recursive = TRUE))
    start2 <- c(1e6, 2e6, 3e6)
    ret <- shaman_score_hic_mat_for_track(db, together, "hic_obs", "hic_exp", "hic_obs", "chr1", 1e6, 2e6, start2, start2 + 1e6,
        expand = 5e5, k = 20
    )
    # the third matrix has too few points, so its file has only the header
    expect_equal(ret, c(1, 1, 0))
    for (s in start2) {
        shaman_score_hic_mat_for_track(db, alone, "hic_obs", "hic_exp", "hic_obs", "chr1", 1e6, 2e6, s, s + 1e6, expand = 5e5, k = 20)
    }
    files <- list.files(alone)
    expect_length(files, 3)
    expect_identical(list.files(together), files)
    expect_identical(unname(tools::md5sum(file.path(together, files))), unname(tools::md5sum(file.path(alone, files))))
    # the files have scores, not only a header
    expect_gt(nrow(data.table::fread(file.path(alone, files[1]))), 1000)
})

test_that("focus points from the regional extraction give the same scores as their own extraction", {
    gsetroot(db)
    focus <- gintervals.2d(1, 2e6, 2.5e6, 1, 2e6, 3e6)
    regional <- gintervals.2d(1, 1.5e6, 3e6, 1, 1.5e6, 3.5e6)
    res <- shaman_score_hic_mat("hic_obs", "hic_exp", focus, regional, k = 20)
    # the points as shaman_score_hic_mat() extracted them before the shortcut
    p <- gextract("hic_obs", focus, colnames = "contacts")
    p <- unique(p[abs(p$start1 - p$start2) > 1024, c("chrom1", "start1", "end1", "chrom2", "start2", "end2")])
    ref <- shaman_score_hic_points("hic_obs", "hic_exp", p, regional, k = 20)
    expect_gt(nrow(ref$points), 1000)
    expect_identical(res$points$score, ref$points$score)
    expect_identical(res, ref)
})
