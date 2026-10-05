# The in-memory count of identical shuffled contacts (sort_uniq = TRUE) against the
# sort | uniq -c | awk pipeline it replaced.
test_that("sort_uniq gives the lines of sort | uniq -c, ordered by start1, start2", {
    skip_if(any(Sys.which(c("sort", "uniq", "awk", "grep")) == ""), "needs sort, uniq, awk and grep")
    set.seed(3)
    work_dir <- tempfile()
    dir.create(work_dir)
    on.exit(unlink(work_dir, recursive = TRUE))
    # a shuffled file in both orientations, with many repeated contacts and coordinates of 1-9 digits
    x <- sample(c(1:50, 2e8L + 1:50), 20000, replace = TRUE)
    y <- sample(c(1:50, 2e8L + 1:50), 20000, replace = TRUE)
    shuf_fn <- file.path(work_dir, "obs_chrX_0_0.shuffled")
    data.table::fwrite(data.frame(start1 = c(x, y), start2 = c(y, x)), shuf_fn, sep = "\t")
    # the shuffled file exists, so this only counts
    shaman_shuffle_hic_mat_for_track(NULL, "obs", work_dir, "chrX", 0, 1, 0, 1, sort_uniq = TRUE)
    new <- readLines(paste0(shuf_fn, ".uniq"))
    old_fn <- file.path(work_dir, "old.uniq")
    system(sprintf("echo 'chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tobs' > %s", old_fn))
    system(sprintf(
        "cat %s | grep -v start | sort -T %s | uniq -c | awk '{ print \"%s\" \"\t\" $2 \"\t\" ($2+1) \"\t\" \"%s\" \"\t\" $3 \"\t\" ($3+1) \"\t\" $1}' >> %s",
        shuf_fn, work_dir, "chrX", "chrX", old_fn
    ))
    old <- readLines(old_fn)
    expect_gt(length(old), 5000)
    expect_identical(new[1], old[1])
    expect_identical(sort(new[-1]), sort(old[-1]))
    u <- data.table::fread(paste0(shuf_fn, ".uniq"))
    expect_identical(order(u$start1, u$start2), seq_len(nrow(u)))
})

test_that("a shuffle with sort_uniq = TRUE returns the seed it used", {
    db <- make_track(tempfile(), c(chr1 = 5e6), c(chr1 = 12000))
    work_dir <- tempfile()
    dir.create(work_dir)
    on.exit(unlink(c(work_dir, db), recursive = TRUE))
    ret <- shaman_shuffle_hic_mat_for_track(db, "hic_obs", work_dir, "chr1", 0, 5e6, 0, 5e6, seed = 7, sort_uniq = TRUE)
    expect_identical(as.integer(ret), 7L)
    expect_true(file.exists(file.path(work_dir, "hic_obs_chr1_0_0.shuffled.uniq")))
})
