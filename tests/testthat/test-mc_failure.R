# Multi-core shuffle where one chromosome's process fails: no track is imported, the shuffled
# chromosomes stay in work_dir, and a rerun completes the track. At most 2 processes.
library(misha)

test_that("a failed chromosome stops the multi-core shuffle before the import, and a rerun completes it", {
    # fewer than 2,000 contacts per chromosome: written unshuffled, which is fast
    db <- make_track(tempfile(), c(chr1 = 2e6, chr2 = 2e6, chr3 = 2e6), c(chr1 = 1500, chr2 = 1500, chr3 = 1500))
    old_opts <- options(shaman.mc_support = 1, shaman.sge_support = 0, shaman.debug = FALSE)
    on.exit(options(old_opts))
    real <- shaman_shuffle_hic_mat_for_track
    for (failure in c("killed", "error")) {
        work_dir <- tempfile()
        dir.create(work_dir)
        failing <- function(track_db, track, work_dir, chrom, ...) {
            if (chrom == "chr2") {
                if (failure == "killed") tools::pskill(Sys.getpid(), tools::SIGKILL) else stop("out of memory")
            }
            real(track_db, track, work_dir, chrom, ...)
        }
        expect_error(
            suppressWarnings(with_mocked_bindings(
                shaman_shuffle_hic_track(db, "hic_obs", work_dir, max_jobs = 2, shuffle = 1, grid_step_iter = 1),
                shaman_shuffle_hic_mat_for_track = failing
            )),
            if (failure == "killed") "the shuffle of chr2 failed \\(no result" else "the shuffle of chr2 failed \\(first error: out of memory\\)"
        )
        gdb.reload()
        expect_false(gtrack.exists("hic_obs_shuffle"))
        done <- file.path(work_dir, paste0("hic_obs_", c("chr1", "chr3"), "_0_0.full_chrom_shuffled"))
        expect_true(all(file.exists(done, paste0(done, ".uniq"))))
        expect_false(file.exists(file.path(work_dir, "hic_obs_chr2_0_0.full_chrom_shuffled")))
        mtimes <- file.mtime(done)
        Sys.sleep(1.1)
        # keep the files this time, to check that chr1 and chr3 are not shuffled again
        options(shaman.debug = TRUE)
        shaman_shuffle_hic_track(db, "hic_obs", work_dir, max_jobs = 2, shuffle = 1, grid_step_iter = 1)
        options(shaman.debug = FALSE)
        gdb.reload()
        expect_setequal(as.character(unique(gextract("hic_obs_shuffle", gintervals.2d.all())$chrom1)), c("chr1", "chr2", "chr3"))
        expect_identical(file.mtime(done), mtimes)
        gtrack.rm("hic_obs_shuffle", force = TRUE)
        unlink(work_dir, recursive = TRUE)
    }
})
