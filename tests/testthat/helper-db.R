# Adds a 2D points track of n random contacts per chromosome (start1 < start2, 2kb-1Mb apart, each
# stored in both orientations as in the example database) to the misha database in dir, creating
# the database with chromosomes of the given sizes (> 1.1Mb) if it does not exist. Returns dir.
make_track <- function(dir, sizes, n, track = "hic_obs", seed = 1) {
    set.seed(seed)
    if (!dir.exists(dir)) {
        dir.create(file.path(dir, "tracks"), recursive = TRUE)
        dir.create(file.path(dir, "seq"))
        writeLines(paste0(sub("^chr", "", names(sizes)), "\t", sizes), file.path(dir, "chrom_sizes.txt"))
        file.create(file.path(dir, "seq", paste0(names(sizes), ".seq")))
    }
    misha::gsetroot(dir)
    d <- do.call(rbind, lapply(names(n)[n > 0], function(chrom) {
        s1 <- sample.int(sizes[[chrom]] - 2^20 - 1, 2 * n[[chrom]], replace = TRUE)
        s2 <- s1 + round(2^stats::runif(2 * n[[chrom]], 11, 20))
        u <- unique(data.frame(s1, s2))[seq_len(n[[chrom]]), ]
        data.frame(chrom = chrom, s1 = u$s1, s2 = u$s2)
    }))
    fn <- tempfile(fileext = ".txt")
    on.exit(unlink(fn))
    utils::write.table(data.frame(
        chrom1 = d$chrom, start1 = c(d$s1, d$s2), end1 = c(d$s1, d$s2) + 1,
        chrom2 = d$chrom, start2 = c(d$s2, d$s1), end2 = c(d$s2, d$s1) + 1, value = 1
    ), fn, sep = "\t", quote = FALSE, row.names = FALSE)
    misha::gtrack.2d.import(track, "random contacts", fn)
    dir
}
