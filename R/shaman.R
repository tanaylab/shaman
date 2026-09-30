#'  generate an expected hic track based on observed hic data
#'
#' \code{shaman_shuffle_hic_track}
#'
#' This function generates an expected 2D hic track based on observed hic data.
#' Each chromosome is shuffled seperately, to generate an expected shuffled contact matrix
#' Note that this function requires sge (qsub) or multicore to be enabled.
#' Parameter can be set via shaman.sge_support or shaman.mc_support in shaman.conf file.
#' Reshuffling of an entire dataset will require 7 hours per 1 billion reads on a machine
#' with one core per chromosome.
#'
#' Each step creates temporary files of the shuffled matrices which are then joined to a track.
#' Temporary files are deleted upon track creation.
#' @param track_db Directory of the misha database.
#' @param obs_track_nm Name of observed 2D genomic track for the hic data.
#' @param work_dir Centralized directory to store temporary files.
#' @param exp_track_nm Name of expected 2D genomic track.
#' @param max_jobs Maximal number of qsub or local jobs - for optimal performance provide the number of chromosomes.
#' @param shuffle Average number of shuffling transitions for each observed point in the chromosomal contact matrix.
#' @param grid_small Initial size of maximum distance between contact pairs consdered for switching
#' @param grid_high Final size of maximum distance between contact pairs consdered for switching
#' @param grid_step_iter Number of iterations in each grid size
#' @param dist_resolution Number of bins in each log2 distance unit. If NA, value is determined
#' based on observed data (recommended).
#' @param smooth Number of bins to use for smoothing the MCMC target function: the decay curve.
#' If NA, value is determined based on observed data (recommended).
#' @param seed Seed for the shuffler: NULL (default) seeds each chromosome from the time (in seconds) its
#' shuffle starts. With a whole number, chromosome i of \code{gintervals.all()} is shuffled with seed + i - 1
#' (modulo 2^31), so the same seed gives the same track. The seed each chromosome was shuffled with is
#' stored in the track attribute \code{seed} (NA for chromosomes this call did not shuffle).
#' A chromosome's seed depends on its position in \code{gintervals.all()}, so the same seed reproduces a
#' track only in a database with the same chromosomes; the \code{seed} attribute reproduces each chromosome.
#' @return No return value, called for side effects: creates the 2D track exp_track_nm in track_db, with the
#' seeds in its attribute \code{seed}.
#'
#' @examples
#'
#' # The example below runs on the test misha db provided with shaman.
#' # Note that this is a toy db sampled from K562 ela data - shuffling the observed track
#' # will not produce the expected track.
#' library(misha)
#' track_db <- shaman_get_test_track_db()
#' gsetroot(track_db)
#' # options(shaman.sge_support=1) #configuring sge engine mode - preferred
#' old_opts <- options(shaman.mc_support = 1) # configuring multi-core mode
#' if (gtrack.exists("hic_obs_shuffle")) {
#'     gtrack.rm("hic_obs_shuffle", force = TRUE)
#'     gdb.reload()
#' }
#' \donttest{
#' ret <- shaman_shuffle_hic_track(track_db,
#'     obs_track_nm = "hic_obs",
#'     # work_dir can be tempdir() only in multi-core mode.
#'     # For sge mode, work_dir must be accessible by all jobs.
#'     work_dir = tempdir(),
#'     shuffle = 1, # default is set to 80
#'     grid_step_iter = 1, # default is set to 40
#'     max_jobs = 2
#' ) # optimally set to number of chromosomes
#' gdb.reload()
#' gtrack.ls("hic_obs_shuffle") # new shuffled track that was created
#' }
#' options(old_opts)
#' @export
##########################################################################################################
shaman_shuffle_hic_track <- function(track_db, obs_track_nm, work_dir,
                                     exp_track_nm = paste0(obs_track_nm, "_shuffle"), max_jobs = 25,
                                     shuffle = 80, grid_small = 500000, grid_high = 1000000, grid_step_iter = 40,
                                     dist_resolution = NA, smooth = NA, seed = NULL) {
    .shaman_check_seed(seed)
    gsetroot(track_db)
    # check tracks
    if (!gtrack.exists(obs_track_nm)) {
        stop(paste("Missing obs_track_nm (", obs_track_nm, ") in track db"))
    }
    if (gtrack.exists(exp_track_nm)) {
        stop(paste("exp_track_nm (", exp_track_nm, ") already exists in track db"))
    }
    # check sge
    .shaman_check_config("shaman.sge_support")
    sge_support <- getOption("shaman.sge_support")
    mc_support <- getOption("shaman.mc_support")
    if (!sge_support & !mc_support) {
        stop(paste("shuffle_hic_track function requires SGE or multicore support. If available, set configuration parameter shaman.sge_support or shaman.mc_support to 1"))
    }
    sge_flags <- getOption("shaman.sge_flags")

    # check work_dir
    if (substr(work_dir, nchar(work_dir), nchar(work_dir)) != "/") {
        work_dir <- paste0(work_dir, "/")
    }
    if (!dir.exists(work_dir)) {
        stop(paste("work_dir (", work_dir, ") does not exists"))
    }

    intervals <- gintervals.all()
    seeds <- if (!is.null(seed)) (seed + seq_len(nrow(intervals)) - 1) %% 2^31

    if (sge_support) {
        commands <- paste0(
            "{library(shaman); shaman_shuffle_hic_mat_for_track(\"", track_db, "\",\"", obs_track_nm, "\",\"",
            work_dir, "\", \"", intervals$chrom, "\", ", intervals$start, ", ",
            intervals$end, ", ", intervals$start, ", ", intervals$end,
            ", min_dist=1024, dist_resolution=", dist_resolution, ", decay_smooth=",
            smooth, ", shuffle=", shuffle, ", grid_small=", grid_small, ", grid_high=", grid_high,
            ", grid_step_iter=", grid_step_iter, ", seed=", if (is.null(seeds)) "NULL" else seeds,
            ", raw_ext=\"full_chrom_raw\", shuffled_ext=\"full_chrom_shuffled\", sort_uniq=TRUE)}"
        )
        res <- .gcluster.run2(command.list = commands, opt.flags = sge_flags, max.jobs = max_jobs)
        # a failed job's retv is the error message
        used_seeds <- sapply(res, function(r) if (is.numeric(r$retv)) r$retv[1] else NA)
    } else {
        doMC::registerDoMC(cores = max_jobs)
        res <- plyr::ddply(intervals, plyr::.(chrom, start), function(x) {
            shaman_shuffle_hic_mat_for_track(track_db, obs_track_nm, work_dir, x$chrom[1],
                x$start[1], x$end[1], x$start[1], x$end[1],
                min_dist = 1024,
                dist_resolution = dist_resolution, decay_smooth = smooth, shuffle = shuffle,
                grid_small = grid_small, grid_high = grid_high, grid_step_iter = grid_step_iter,
                raw_ext = "full_chrom_raw", shuffled_ext = "full_chrom_shuffled", sort_uniq = TRUE,
                seed = seeds[match(x$chrom[1], intervals$chrom)]
            )
        }, .parallel = TRUE)
        used_seeds <- res$V1[match(intervals$chrom, res$chrom)]
    }
    exp_shuf_files <- paste0(obs_track_nm, "_", intervals$chrom, "_0_0.full_chrom_shuffled.uniq")
    obs_shuf_files <- list.files(work_dir, pattern = paste0(obs_track_nm, ".*full_chrom_shuffled.uniq"))
    missing_files <- exp_shuf_files[!exp_shuf_files %in% obs_shuf_files]
    if (length(missing_files) > 0) {
        warning(paste(length(missing_files), "full chrom files were not shuffled:\n", paste(missing_files, collapse = ",")))
    }

    # reaching this point means that all matrices should be in work_dir --> creating track this section
    # should be replaced after misha generates a proper import for overlapping contacts
    files <- c()
    for (chrom in intervals$chrom) {
        fn <- paste0(work_dir, obs_track_nm, "_", chrom, "_0_0.full_chrom_shuffled.uniq")
        message(chrom)
        if (file.exists(fn)) {
            files <- c(files, fn)
        } else {
            message("missing file")
        }
    }
    gtrack.2d.import(exp_track_nm, paste(
        "shuffled 2d track with shuffle factor =",
        shuffle, ", based on", obs_track_nm
    ), files)
    gtrack.attr.set(exp_track_nm, "seed", paste0(intervals$chrom, ":", used_seeds, collapse = " "))

    # cleanup work dir from all temporary files
    debug <- getOption("shaman.debug")
    if (!debug) {
        try(system(sprintf("rm %s%s*", work_dir, obs_track_nm)))
    }
}



####################################################################################################
#' Generate an expected matrix from observed data as a process for generating an expected track
#'
#' \code{shuffle_hic_mat_for_track}
#'
#' This function generates an expected 2D hic matrix from observed hic data. Should not be called externally,
#  but rather as jobs when creating a 2D shuffled track.
#' The observed data is a combination of observed contacts in scope plus already shuffled
#' near-cis contacts (stored in work_dir) which we sample from to maintain the decay probability curve.
#'
#' Each step creates temporary files of the shuffled matrices which are then joined to a track.
#' Temporary files are deleted upon track creation.
#' @param track_db Directory of the misha database.
#' @param track Name of observed 2D genomic track for the hic data.
#' @param work_dir Centralized directory to store temporary files.
#' @param chrom The chormosome of the matrix.
#' @param start1 The start coordinate of the first dimension.
#' @param end1 The end coordinate of the first dimension.
#' @param start2 The start coordinate of the second dimension.
#' @param end2 The end coordinate of the second dimension.
#' @param min_dist The minimum distance between contact end points.
#' @param max_dist The maximum distance between contact end points.
#' @param dist_resolution Number of bins in each log2 distance unit. If NA, value is determined
#' based on observed data (recommended).
#' @param decay_smooth Number of bins to use for smoothing the MCMC target function: the decay curve.
#' If NA, value is determined based on observed data (recommended).
#' @param proposal_iterations Number of MCMC sampling iterations between proposal corrections.
#' @param shuffle Number of shuffling rounds for each observed point.
#' @param hic_mcmc_max_resolution Maximum number of bins for each log2 unit
#' @param raw_ext File extension of the observed data.
#' @param shuffled_ext File extension of the shuffled data.
#' @param grid_small Initial size of maximum distance between contact pairs consdered for switching
#' @param grid_high Final size of maximum distance between contact pairs consdered for switching
#' @param grid_increase Grid increase size
#' @param grid_step_iter Number of iterations in each grid size
#' @param sort_uniq Binary flag, indicating whether the shuffled matrix file should be sorted and
#' contacts combined. This is required prior to importing the track to misha, and should be applied
#' to full chromosomes only.
#' @param seed Seed for the shuffler: NULL (default) seeds it from the current time (in seconds), or a
#' whole number between 0 and 2^31 - 1. The shuffler prints the seed it uses.
#'
#' @return The seed the shuffler used, or NA if the matrix was not shuffled here (no or too few
#' contacts, or the shuffled file already existed).
#'
#' @export
##########################################################################################################
shaman_shuffle_hic_mat_for_track <- function(track_db, track, work_dir, chrom, start1, end1, start2, end2,
                                             min_dist = 1024, max_dist = max(gintervals.all()$end), dist_resolution = NA,
                                             decay_smooth = NA, proposal_iterations = 1e+07, shuffle = 80, hic_mcmc_max_resolution = 400,
                                             raw_ext = "raw", shuffled_ext = "shuffled", grid_small = 500000, grid_high = 1000000, grid_increase = 500000,
                                             grid_step_iter = 40, sort_uniq = FALSE, seed = NULL) {
    seed <- .shaman_check_seed(seed)
    raw_fn <- paste0(work_dir, "/", track, "_", chrom, "_", start1, "_", start2, ".", raw_ext)
    shuf_fn <- paste0(work_dir, "/", track, "_", chrom, "_", start1, "_", start2, ".", shuffled_ext)
    ret <- NA
    if (!file.exists(shuf_fn)) {
        x <- sample(1:10, 1)
        system(paste("sleep", x))
        old_opts <- options(gmultitasking = FALSE, gmax.data.size = 1e+09)
        on.exit(options(old_opts), add = TRUE)
        gsetroot(track_db)
        scope <- gintervals.2d(chrom, start1, end1, chrom, start2, end2)
        a <- gextract(track, scope, band = c(-max_dist, -min_dist + 1), colnames = "contact")
        if (is.null(nrow(a))) {
            message("not shuffling, no data")
            system(paste("echo 'start1\tstart2' > ", shuf_fn))
            return(ret)
        }
        # multiply each line according to the number of observed counts
        a <- plyr::ddply(a, c("contact"), function(x) {
            return(x[rep(seq_len(nrow(x)), each = x$contact[1]), c("start1", "start2")])
        })[, 2:3]

        # automatic decision on resolution, smoothing and samples per correction
        if (is.na(dist_resolution)) {
            dist_resolution <- min(floor((nrow(a) / 200) / (log2(max_dist) - log2(min_dist))), hic_mcmc_max_resolution)
        }

        message(paste("writing", nrow(a), "raw data, dist resolution=", dist_resolution))
        if (nrow(a) < 2000 | dist_resolution == 0) {
            message("not shuffling, leaving raw")
            data.table::fwrite(format(rbind(a, setNames(rev(a), names(a))), scientific = FALSE), shuf_fn,
                quote = FALSE, row.names = FALSE,
                sep = "\t"
            )
        } else {
            if (is.na(decay_smooth)) {
                decay_smooth <- min(floor(dist_resolution / 10), 20)
            }

            ret <- shaman_hic_matrix_shuffler_cpp(
                t(a[, c("start1", "start2")]),
                shuf_fn, shuffle, 1, 0.5, dist_resolution, decay_smooth, 5, 0.25,
                max_dist, 1024, 1, grid_small, grid_high, grid_increase, grid_step_iter, 0, 1, seed
            )
        }
    }
    if (sort_uniq) {
        # count identical contacts in memory (the lines sort | uniq -c gave, ordered by start1, start2),
        # in one thread like the shuffle itself (more threads gain little here)
        threads <- data.table::setDTthreads(1)
        on.exit(data.table::setDTthreads(threads), add = TRUE)
        a <- data.table::fread(shuf_fn, sep = "\t", header = TRUE, colClasses = "integer")
        data.table::setorderv(a)
        id <- data.table::rleidv(a)
        obs <- tabulate(id, max(0L, id)) # max(0L, ...): no contacts give no bins, not one empty bin
        rm(id)
        first <- cumsum(obs) - obs + 1L
        start1 <- a$start1[first]
        start2 <- a$start2[first]
        rm(a)
        chroms <- rep(chrom, length(obs))
        data.table::fwrite(list(
            chrom1 = chroms, start1 = start1, end1 = start1 + 1L,
            chrom2 = chroms, start2 = start2, end2 = start2 + 1L, obs = obs
        ), paste0(shuf_fn, ".uniq"), sep = "\t")
    }
    return(ret)
}

# Checks a shuffler seed (NULL or a whole number in [0, 2^31 - 1]) and returns it as the integer
# the C++ shuffler takes: -1 for NULL, which seeds from the current time.
.shaman_check_seed <- function(seed) {
    if (is.null(seed)) {
        return(-1L)
    }
    if (!is.numeric(seed) || length(seed) != 1 || is.na(seed) || seed != round(seed) ||
        seed < 0 || seed > .Machine$integer.max) {
        stop("seed must be NULL or a whole number between 0 and ", .Machine$integer.max)
    }
    return(as.integer(seed))
}

##########################################################################################################
#'  generate a score hic track based on observed and expected (shuffled) hic data
#'
#' \code{shaman_score_hic_track}
#'
#' This function generates a 2D score track based on observed and expected hic data.
#' The score is computed by generating a grid of small matrices spanning all chromosomes
#' and computing the score of each matrix independantly.
#' The model for computing the score relies on the KS D statistic computed for each observed point,
#' over the distances of the k-nearest neighbors in the observed compared to the expected.
#' High scores represent contact enrichment while low scores depict insulation.
#' Note that this function requires either sge (qsub) or multicore to compute in a timely manner.
#' Parameters can be set via shaman.sge_support or shaman.mc_support in shaman.conf file.
#' Score computation on 1 billion reads on a distributed system may take 4-10 hours (with default parameters),
#' depending on the number of cores available.
#' \code{options(shaman.score.threads = N)} computes the kNN distances and scores of each matrix on N
#' threads (default 1); in multi-core mode each of the max_jobs processes uses N threads, and in SGE mode
#' each job does, so shaman.sge_flags should ask for N slots.
#' Each extraction of a 2D track can keep about the whole chromosome's track file in memory;
#' misha's \code{options(gtrack.num.chunks = 1000)} bounds that (2-3x lower peak memory on dense
#' matrices in our runs) without changing the scores.
#'
#' Each step creates temporary files of the matrix scores which are then joined to a track.
#' Temporary files are deleted upon track creation.
#' Matrices that still have no score after 3 rounds of jobs stop the run with an error that names them.
#' The finished score files stay in work_dir, so a rerun with the same work_dir computes only the missing ones.
#' @param track_db Directory of the misha database.
#' @param work_dir Centralized directory to store temporary files.
#' @param score_track_nm Score track that will be created.
#' @param obs_track_nms Names of observed 2D genomic tracks for the hic data. Pooling of multiple
#' observed tracks is supported.
#' @param exp_track_nms Names of expected (shuffled) 2D genomic tracks. Pooling of multiple expected
#' tracks is supported.
#' @param points_track_nms Names of 2D genomic tracks that contain points on which to compute
#' normalized score. Pooling points from multiple tracks is supported.
#' @param near_cis Size of matrix in grid.
#' @param expand Size of expansion, points to include outside the matrix for accurate computing of the score.
#' Note that for each observed point, its k-nearest neighbors must be included in the expanded matrix.
#' @param k The number of neighbor distances used for the score. For higher resolution maps, increase k. For
#' lower resolution maps, decrease k.
#' @param max_jobs Maximal number of qsub jobs.
#' @return No return value, called for side effects: creates the 2D score track score_track_nm in track_db.
#'
#' @examples
#'
#' # The example below runs on the test misha db provided with shaman.
#' # Note that this is a toy db sampled from K562 ela data -
#' # scoring based on the observed and expected tracks will not produce the score track,
#' # as most of the genome is missing.
#' library(misha)
#' track_db <- shaman_get_test_track_db()
#' gsetroot(track_db)
#' # options(shaman.sge_support=1) #configuring sge engine mode - preferred
#' old_opts <- options(shaman.mc_support = 1) # configuring multi-core mode
#' if (gtrack.exists("hic_score_new")) {
#'     gtrack.rm("hic_score_new", force = TRUE)
#'     gdb.reload()
#' }
#' ret <- shaman_score_hic_track(track_db,
#'     # work_dir can be tempdir() only in multi-core mode.
#'     # For sge mode, work_dir must be accessible by all jobs.
#'     work_dir = tempdir(),
#'     score_track_nm = "hic_score_new",
#'     obs_track_nms = "hic_obs",
#'     exp_track_nms = "hic_exp",
#'     near_cis = 1e09, # this test db contains very little data, can increase the size of each job
#'     k = 20, # default is set to 100
#'     max_jobs = 2
#' ) # increase number of jobs for optimal runtime when running in sge mode
#' gdb.reload()
#' gtrack.ls("hic_score_new") # new score track that was created
#' options(old_opts)
#' @export
##########################################################################################################
shaman_score_hic_track <- function(track_db, work_dir, score_track_nm, obs_track_nms,
                                   exp_track_nms = paste0(obs_track_nms, "_shuffle"), points_track_nms = obs_track_nms,
                                   near_cis = 5e06, expand = 2e06, k = 100, max_jobs = 100) {
    gsetroot(track_db)
    # check tracks
    if (sum(gtrack.exists(obs_track_nms)) < length(obs_track_nms)) {
        stop(paste("Missing obs_track_nm (", obs_track_nms[!gtrack.exists(obs_track_nms)], ") in track db"))
    }
    if (sum(gtrack.exists(exp_track_nms)) < length(exp_track_nms)) {
        stop(paste("Missing exp_track_nm (", exp_track_nms[!gtrack.exists(exp_track_nms)], ") in track db"))
    }
    if (gtrack.exists(score_track_nm)) {
        stop(paste("score_track_nm (", score_track_nm, ") already exists in track db"))
    }
    # check sge
    .shaman_check_config("shaman.sge_support")
    sge_support <- getOption("shaman.sge_support")
    mc_support <- getOption("shaman.mc_support")
    if (!sge_support & !mc_support) {
        stop(paste("score_hic_track function requires SGE or multicore support. If available, set configuration parameter shaman.sge_support or shaman.mc_support to 1"))
    }
    sge_flags <- getOption("shaman.sge_flags")

    # split all chromosomes to rectangles of size near_cis*near_cis
    band <- c(-max(gintervals.all()$end), 0)
    near_cis_intervals <- gintervals.force_range(plyr::ddply(gintervals.all(), c("chrom"), function(x) {
        data.frame(
            chrom = rep(x$chrom[1], nrow(x)),
            start = seq(0, floor(x$end / near_cis) * near_cis, by = near_cis),
            end = seq(near_cis, ceiling(x$end / near_cis) * near_cis, by = near_cis)
        )
    }))
    g <- expand.grid(1:nrow(near_cis_intervals), 1:nrow(near_cis_intervals))
    g <- g[near_cis_intervals$chrom[g$Var1] == near_cis_intervals$chrom[g$Var2], ]
    near_cis_2d <- gintervals.2d(
        near_cis_intervals$chrom[g$Var1], near_cis_intervals$start[g$Var1], near_cis_intervals$end[g$Var1],
        near_cis_intervals$chrom[g$Var2], near_cis_intervals$start[g$Var2], near_cis_intervals$end[g$Var2]
    )

    # selecting only the upper half of the matrix
    near_cis_2d_upper_mat <- gintervals.2d.band_intersect(near_cis_2d, band)
    expected_files <- paste0(
        work_dir, "/", paste0(obs_track_nms, collapse = "."), ".",
        near_cis_2d_upper_mat$chrom1, ".", near_cis_2d_upper_mat$start1, ".", near_cis_2d_upper_mat$start2, ".score"
    )

    existing_files <- file.exists(expected_files)
    missing_files <- expected_files[!existing_files]
    message(paste("missing", length(missing_files), "score files"))
    near_cis_2d_upper_mat <- near_cis_2d_upper_mat[!existing_files, ]
    expected_files <- paste0(
        work_dir, "/", paste0(obs_track_nms, collapse = "."), ".",
        near_cis_2d_upper_mat$chrom1, ".", near_cis_2d_upper_mat$start1, ".", near_cis_2d_upper_mat$start2, ".score"
    )
    # ponytail: a fixed number of rounds, not an option; a matrix that failed this often (e.g. a job
    # that is always killed for memory) would fail again
    max_rounds <- 3
    rounds <- 0
    while (nrow(near_cis_2d_upper_mat) > 0) {
        if (rounds == max_rounds) {
            m <- near_cis_2d_upper_mat
            failed <- sprintf("%s:%.0f-%.0f x %.0f-%.0f", m$chrom1, m$start1, m$end1, m$start2, m$end2)
            stop(sprintf(
                "%d matrices have no score after %d rounds: %s%s. The scores of the other matrices are kept in %s; a rerun with the same work_dir computes only the missing ones.",
                length(failed), max_rounds, paste(failed[seq_len(min(length(failed), 20))], collapse = ", "),
                if (length(failed) > 20) ", ..." else "", work_dir
            ))
        }
        rounds <- rounds + 1
        # compute scores for each of the small matrices
        if (sge_support) {
            # gcluster.run jobs do not inherit the misha root, so each job sets it
            commands <- paste0(
                "{library(shaman); gsetroot(track_db); shaman_score_hic_mat_for_track(track_db, work_dir, obs_track_nms, exp_track_nms, points_track_nms, \"",
                near_cis_2d_upper_mat$chrom1, "\", ", near_cis_2d_upper_mat$start1, ", ",
                near_cis_2d_upper_mat$end1, ",", near_cis_2d_upper_mat$start2, ", ",
                near_cis_2d_upper_mat$end2, ", ", expand, ", ", k, ")}"
            )
            # commands <- paste(commands, collapse=",")

            res <- .gcluster.run2(command.list = commands, opt.flags = sge_flags, max.jobs = max_jobs)
        } else {
            doMC::registerDoMC(cores = max_jobs)
            res <- plyr::ddply(near_cis_2d_upper_mat, plyr::.(chrom1, start1, start2), function(x) {
                shaman_score_hic_mat_for_track(
                    track_db, work_dir, obs_track_nms, exp_track_nms, points_track_nms,
                    x$chrom1[1], x$start1[1], x$end1[1], x$start2[1], x$end2[1], expand, k
                )
            },
            .parallel = TRUE
            )
        }
        # res <- eval(parse(text=paste("gcluster.run(", commands, ",opt.flags=\"", sge_flags,  "\" ,max.jobs=", max_jobs, ")")))
        # check to see if there are any missing files
        existing_files <- file.exists(expected_files)
        missing_files <- expected_files[!existing_files]
        message(paste("missing", length(missing_files), "score files"))
        near_cis_2d_upper_mat <- near_cis_2d_upper_mat[!existing_files, ]
        expected_files <- paste0(
            work_dir, "/", paste0(obs_track_nms, collapse = "."), ".",
            near_cis_2d_upper_mat$chrom1, ".", near_cis_2d_upper_mat$start1, ".", near_cis_2d_upper_mat$start2, ".score"
        )
    }
    # reaching this point means that all points have a score - need to create track
    score_files <- list.files(work_dir, paste0(obs_track_nms, ".*.score$"), full.names = TRUE)

    gtrack.2d.import_contacts(score_track_nm, paste(
        "normalized score of 2d track with k =", k, "based on",
        paste(obs_track_nms, collapse = ", "), "and",
        paste(exp_track_nms, collapse = ", ")
    ), score_files, allow.duplicates = FALSE)

    # cleanup work dir from all temporary files
    for (track in obs_track_nms) {
        try(system(sprintf("rm %s/%s*", work_dir, track)))
    }
}

##########################################################################################################
#'  generate a score matrix for observed data based on the expected
#'
#' \code{shaman_score_hic_mat_for_track}
#'
#' This function extracts observed data and expected data in an expanded matrix and computes
#  the score for each observed point.
#' The score for a point is the KS D-statistic of the distances to the points k-nearest-neighbors
#  in the observed data compared the the expected data.
#'
#' @param track_db Directory of the misha database.
#' @param work_dir Centralized directory to store temporary files.
#' @param obs_track_nms Names of observed 2D genomic tracks for the hic data. Pooling of multiple
#' observed tracks is supported.
#' @param exp_track_nms Names of expected (shuffled) 2D genomic tracks. Pooling of multiple expected
#' tracks is supported.
#' @param points_track_nms Names of 2D genomic tracks that contain points on which to compute
#' normalized score. Pooling points from multiple tracks is supported.
#' @param chrom The chormosome of the matrix.
#' @param start1 The start coordinate of the first dimension.
#' @param end1 The end coordinate of the first dimension.
#' @param start2 The start coordinate of the second dimension.
#' @param end2 The end coordinate of the second dimension.
#' @param expand Size of expansion, points to include outside the matrix for accurate computing of the score.
#' Note that for each observed point, its k-nearest neighbors must be included in the expanded matrix.
#' @param k The number of neighbor distances used for the score. For higher resolution maps, increase k. For
#' lower resolution maps, decrease k.
#' @param min_dist The minimum distance between points.
#'
#' @details chrom, start1, end1, start2 and end2 can be vectors, to score several matrices in one call.
#' Matrices with the same chrom, start1 and end1 then share one extraction of each track over the union
#' of their expanded intervals, instead of one extraction per matrix (each extraction of a 2D track reads
#' the whole chromosome pair). The output files are the same as scoring the matrices one by one.
#'
#' \code{options(shaman.score.threads = N)} computes the kNN distances and scores on N threads (default 1).
#'
#' @return 0, 1 or -1 per matrix.
#' @export
##########################################################################################################
shaman_score_hic_mat_for_track <- function(track_db, work_dir, obs_track_nms, exp_track_nms, points_track_nms,
                                           chrom, start1, end1, start2, end2, expand = 2e06, k = 100, min_dist = 1024) {
    if (max(length(chrom), length(start1), length(end1), length(start2), length(end2)) > 1) {
        rects <- data.frame(chrom = as.character(chrom), start1 = start1, end1 = end1, start2 = start2, end2 = end2, stringsAsFactors = FALSE)
        fns <- paste0(work_dir, "/", paste0(obs_track_nms, collapse = "."), ".", rects$chrom, ".", rects$start1, ".", rects$start2, ".score")
        row <- paste(rects$chrom, rects$start1, rects$end1)
        todo <- !file.exists(fns)
        ret <- rep(0, nrow(rects))
        on.exit(rm(list = ls(.shaman_extract_cache), envir = .shaman_extract_cache))
        old_opts <- options(gmax.data.size = 1e09)
        on.exit(options(old_opts), add = TRUE)
        for (r in unique(row[todo])) {
            i <- which(row == r & todo)
            if (length(i) > 1) {
                .shaman_cache_extractions(unique(c(obs_track_nms, exp_track_nms, points_track_nms)), gintervals.force_range(data.frame(
                    chrom1 = rects$chrom[i[1]], start1 = rects$start1[i[1]] - expand, end1 = rects$end1[i[1]] + expand,
                    chrom2 = rects$chrom[i[1]], start2 = min(rects$start2[i]) - expand, end2 = max(rects$end2[i]) + expand
                )))
            }
            for (j in i) {
                ret[j] <- shaman_score_hic_mat_for_track(track_db, work_dir, obs_track_nms, exp_track_nms, points_track_nms,
                    rects$chrom[j], rects$start1[j], rects$end1[j], rects$start2[j], rects$end2[j], expand = expand, k = k, min_dist = min_dist)
            }
            rm(list = ls(.shaman_extract_cache), envir = .shaman_extract_cache)
        }
        return(ret)
    }
    fn <- paste0(work_dir, "/", paste0(obs_track_nms, collapse = "."), ".", chrom, ".", start1, ".", start2, ".score")
    if (file.exists(fn)) {
        return(0)
    }
    old_opts <- options(gmax.data.size = 1e09)
    on.exit(options(old_opts), add = TRUE)
    regional_interval <- gintervals.force_range(data.frame(
        chrom1 = chrom, start1 = start1 - expand, end1 = end1 + expand,
        chrom2 = chrom, start2 = start2 - expand, end2 = end2 + expand
    ))
    focus_interval <- gintervals.2d(chrom, start1, end1, chrom, start2, end2)
    n <- shaman_score_hic_mat(obs_track_nms, exp_track_nms, focus_interval, regional_interval, points_track_nms = points_track_nms, min_dist = min_dist, k = k)
    if (is.null(n)) {
        system(paste("echo 'chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tscore' > ", fn))
        # insufficient data in region - not writing region file
        return(0)
    }
    if (length(n) < 1) {
        # knn did not complete
        return(-1)
    }
    p <- n$points
    data.table::fwrite(format(p[
        p$start1 <= p$start2,
        c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "score")
    ],
    scientific = FALSE
    ),
    fn,
    row.names = FALSE, quote = FALSE, sep = "\t"
    )
    return(1)
}

##########################################################################################################
#'  generate a score matrix for observed data based on the expected for a given 2D focus interval
#'
#' \code{shaman_score_hic_mat}
#'
#' This function extracts observed data and expected data in an expanded matrix and computes
#  the score for each observed point.
#' The score for a point is the KS D-statistic of the distances to the points k-nearest-neighbors
#  in the observed data compared the the expected data.
#'
#' @param obs_track_nms Names of observed 2D genomic tracks for the hic data. Pooling of multiple
#' observed tracks is supported.
#' @param exp_track_nms Names of expected (shuffled) 2D genomic tracks. Pooling of multiple expected
#' tracks is supported.
#' @param focus_interval 2D interval on which to compute the scores.
#' @param regional_interval An expansion of the focus interval, inclusing  points outside the focus matrix
#' for accurate computing of the score. Note that for each observed point, its k-nearest neighbors must be
#' included in the expanded matrix.
#' @param points_track_nms Names of 2D genomic tracks that contain points on which to compute
#' normalized score. Pooling points from multiple tracks is supported.
#' @param min_dist The minimum distance between points.
#' @param k The number of neighbor distances used for the score. For higher resolution maps, increase k. For
#' lower resolution maps, decrease k.
#' @param k_exp The number of neighbor distances used for the score on the expected tracks. Note, that when
#' comparing expected generated by shuffling the observed, k_exp should be 2*k as the number of contacts in
#' the expected track will always be twice the observed. However, if comparing between two datasets that are
#' independent, k_exp should be set to NA and the will be determined by the ratio between the total number
#' of contacts in this region.
#'
#' @return NULL if insufficient observed data, otherwise resturns a list containing 3 elements:
#' 1) points - start1, start2 and score for all observed points.
#' 2) obs - the observed points.
#' 3) exp - the expected points.
#'
#' @examples
#'
#' # Set misha db to test
#' library(misha)
#' gsetroot(shaman_get_test_track_db())
#' mat_score <- shaman_score_hic_mat(
#'     obs_track_nms = "hic_obs", exp_track_nms = "hic_exp",
#'     focus_interval = gintervals.2d(2, 176.7e06, 176.8e06, 2, 176.7e06, 176.8e06),
#'     regional_interval = gintervals.2d(2, 176.5e06, 177e06, 2, 176.5e06, 177e06)
#' )
#' shaman_gplot_map_score(mat_score$points)
#' @export
##########################################################################################################
shaman_score_hic_mat <- function(obs_track_nms, exp_track_nms, focus_interval, regional_interval,
                                 points_track_nms = obs_track_nms, min_dist = 1024, k = 100, k_exp = 2 * k) {
    if (identical(points_track_nms, obs_track_nms) && .shaman_interval_within(focus_interval, regional_interval) &&
        all(vapply(obs_track_nms, function(x) gtrack.info(x)$type == "points", TRUE))) {
        # The points are then the observed contacts inside the focus interval: take them from the
        # observed contacts of the regional interval (see .shaman_points_in) instead of another pass
        # over the track.
        obs <- .shaman_combine_points_multi_tracks(obs_track_nms, regional_interval, min_dist)
        if (NROW(obs) < 1000) {
            # then so are the points in the focus interval
            message("number of points in focus interval < 1000")
            return(NULL)
        }
        points <- .shaman_points_in(obs, focus_interval)
    } else {
        obs <- NULL
        points <- .shaman_combine_points_multi_tracks(points_track_nms, focus_interval, min_dist)
    }
    if (is.null(points)) {
        message("number of points in focus interval = 0")
        return(NULL)
    }
    if (nrow(points) < 1000) {
        message("number of points in focus interval < 1000")
        return(NULL)
    }

    points <- unique(points[, c("chrom1", "start1", "end1", "chrom2", "start2", "end2")])
    message(paste0("kk norm on ", nrow(points), " points"))
    return(.shaman_score_hic_points(obs_track_nms, exp_track_nms, points, regional_interval, min_dist, k = k, k_exp = k_exp, obs = obs))
}

##########################################################################################################
#'  generate a score matrix for observed data based on the expected for a given set of points
#'
#' \code{shaman_score_hic_points}
#'
#' This function extracts observed data and expected data in an expanded matrix and computes
#  the score for each observed point.
#' The score for a point is the KS D-statistic of the distances to the points k-nearest-neighbors
#  in the observed data compared the the expected data.
#'
#' @param obs_track_nms Names of observed 2D genomic tracks for the hic data. Pooling of multiple
#' observed tracks is supported.
#' @param exp_track_nms Names of expected (shuffled) 2D genomic tracks. Pooling of multiple expected
#' tracks is supported.
#' @param points A score will be computed for each of the points.
#' @param regional_interval An expansion of the focus interval, inclusing  points outside the focus matrix
#' for accurate computing of the score. Note that for each observed point, its k-nearest neighbors must be
#' included in the expanded matrix.
#' @param min_dist The minimum distance between points.
#' @param k The number of neighbor distances used for the score. For higher resolution maps, increase k. For
#' lower resolution maps, decrease k.
#' @param k_exp The number of neighbor distances used for the score on the expected tracks (see
#' \code{shaman_score_hic_mat}).
#'
#' @return NULL if insufficient observed data, otherwise resturns a list containing 3 elements:
#' 1) points - start1, start2 and score for all observed points.
#' 2) obs - the observed points.
#' 3) exp - the expected points.
#'
#' @examples
#'
#' # Set misha db to test
#' library(misha)
#' gsetroot(shaman_get_test_track_db())
#' focus <- gintervals.2d(2, 176.7e06, 176.8e06, 2, 176.7e06, 176.8e06)
#' points <- gextract("hic_obs", focus, band = c(-5e05, -1024))
#' mat_score <- shaman_score_hic_points(
#'     obs_track_nms = "hic_obs", exp_track_nms = "hic_exp",
#'     points = points, regional_interval = gintervals.2d(2, 176.5e06, 177e06, 2, 176.5e06, 177e06)
#' )
#' shaman_gplot_map_score(mat_score$points)
#' @export
##########################################################################################################
shaman_score_hic_points <- function(obs_track_nms, exp_track_nms, points, regional_interval, min_dist = 1024, k = 100, k_exp = 2 * k) {
    .shaman_score_hic_points(obs_track_nms, exp_track_nms, points, regional_interval, min_dist, k = k, k_exp = k_exp)
}

# obs: the observed contacts in regional_interval as .shaman_combine_points_multi_tracks() returns
# them, if already extracted
.shaman_score_hic_points <- function(obs_track_nms, exp_track_nms, points, regional_interval, min_dist, k, k_exp, obs = NULL) {
    message(paste("obs = ", paste(obs_track_nms, collapse = ",")))
    if (is.null(obs)) {
        obs <- .shaman_combine_points_multi_tracks(obs_track_nms, regional_interval, min_dist)
    }
    if (is.null(obs) | nrow(points) == 0) {
        message(paste("0 data found in intervals, focus interval=", nrow(points)))
        return(NULL)
    }
    # repeat each point by its number of contacts, in the row order plyr::ddply(obs, "contacts", ...) gave
    o <- order(obs$contacts)
    obs <- as.data.frame(lapply(obs, `[`, rep(o, obs$contacts[o])))
    if (nrow(obs) < k) {
        message(paste("insufficient data found in intervals: obs=", nrow(obs)))
        return(NULL)
    }
    n_obs <- nrow(obs)

    exp <- .shaman_combine_points_multi_tracks(exp_track_nms, regional_interval, min_dist)
    if (is.null(exp)) {
        message(paste("0 data found in intervals: exp"))
        return(NULL)
    }
    # repeat each point by its number of contacts, in the row order plyr::ddply(exp, "contacts", ...) gave
    o <- order(exp$contacts)
    exp <- as.data.frame(lapply(exp, `[`, rep(o, exp$contacts[o])))
    if (nrow(exp) < k) {
        message(paste("insufficient data found in intervals: exp=", nrow(exp)))
        return(NULL)
    }
    n_exp <- nrow(exp)
    if (is.na(k_exp)) {
        # computing k_exp by number of points
        k_exp <- round(k * n_exp / n_obs)
    }
    if (n_exp < k_exp) {
        message(paste("insufficient data found in intervals: exp=", n_exp, "< k_exp =", k_exp))
        return(NULL)
    }
    message(paste0("n_obs = ", n_obs, ", n_exp = ", n_exp, ", k_exp = ", k_exp))
    # same values as shaman_merge_ks_cpp(round(RANN::nn2(obs, points, k)$nn.dist),
    # round(RANN::nn2(exp, points, k_exp)$nn.dist)), without the n x k matrices
    s_ks <- shaman_knn_ks_cpp(
        obs$start1, obs$start2, exp$start1, exp$start2, points$start1, points$start2,
        k, k_exp, getOption("shaman.score.threads", 1)
    )
    rm(obs, exp)

    points$score <- 100 * ifelse(-s_ks$V1 < s_ks$V2, s_ks$V2, s_ks$V1)

    return(list(points = points))
}

#########################################################################################################
#'  Inline function for generating an expected matrix and computing the score for a given interval
#'
#' \code{shaman_shuffle_and_score_hic_mat}
#'
#' This function generates an expected 2D hic matrix based on observed hic data, and computes its score.
#' @param obs_track_nms Name of observed 2D genomic tracks for the hic data.
#' @param interval 2D interval on which to compute the scores.
#' @param work_dir Centralized directory to store temporary files.
#' @param expand Size of expansion, points to include outside the matrix for accurate computing of the score.
#' Note that for each observed point, its k-nearest neighbors must be included in the expanded matrix.
#' @param min_dist The minimum distance between points.
#' @param k The number of neighbor distances used for the score. For higher resolution maps, increase k. For
#' @param dist_resolution Number of bins in each log2 distance unit. If NA, value is determined
#' based on observed data (recommended).
#' @param decay_smooth Number of bins to use for smoothing the MCMC target function: the decay curve.
#' If NA, value is determined based on observed data (recommended).
#' @param hic_mcmc_max_resolution Maximum number of bins for each log2 unit.
#' @param shuffle Number of shuffling rounds for each observed point.
#' @param grid_small Initial size of maximum distance between contact pairs consdered for switching
#' @param grid_high Final size of maximum distance between contact pairs consdered for switching
#' @param grid_increase Grid increase size
#' @param grid_step_iter Number of iterations in each grid size
#' @param seed Seed for the shuffler: NULL (default) seeds it from the current time (in seconds), or a
#' whole number between 0 and 2^31 - 1. The shuffler prints the seed it uses.

#' @return NULL if insufficient observed data, otherwise resturns a list containing:
#' 1) points - start1, start2 and score for all observed points.
#' 2) obs - the observed points.
#' 3) exp - the expected points.
#' 4) exp_fn - the name of the expected (shuffled) data file
#' 5) seed - the seed the shuffler used
#'
#' @examples
#'
#' # Set misha db to test
#' library(misha)
#' gsetroot(shaman_get_test_track_db())
#' mat_score <- shaman_shuffle_and_score_hic_mat(
#'     obs_track_nms = "hic_obs",
#'     interval = gintervals.2d(2, 176.6e06, 176.9e06, 2, 176.6e06, 176.9e06),
#'     expand = 1e05,
#'     work_dir = tempdir(),
#'     shuffle = 2, # default is set to 80
#'     grid_step_iter = 1 # default is set to 40
#' )
#' shaman_gplot_map_score(mat_score$points)
#' @export
##########################################################################################################

shaman_shuffle_and_score_hic_mat <- function(obs_track_nms, interval, work_dir, expand = 1e06, min_dist = 1024, k = 100,
                                             dist_resolution = NA, decay_smooth = NA, hic_mcmc_max_resolution = 400, shuffle = 80,
                                             grid_small = 500000, grid_high = 1000000, grid_increase = 500000,
                                             grid_step_iter = 40, seed = NULL) {
    seed <- .shaman_check_seed(seed)
    if (interval$chrom1 != interval$chrom2) {
        stop("Only cis intervals supported")
    }
    points <- .shaman_combine_points_multi_tracks(obs_track_nms, interval, min_dist)
    points <- points[points$start2 > points$start1, ]
    if (is.null(points)) {
        message("number of points in interval = 0")
        return(NULL)
    }
    message(paste("working on ", nrow(points), "points"))
    if (is.na(dist_resolution)) {
        dist_resolution <- min(floor((nrow(points) / 200) / (log2(max(points$start2 - points$start1)) -
            log2(1024))), hic_mcmc_max_resolution)
    }
    # adjust expand such that first and last dist bins are complete
    message(paste("adjusting expand"))
    max_d <- interval$end2 - interval$start1
    expand <- (2**(ceiling(log2(max_d + 2 * expand) * dist_resolution) / dist_resolution) - max_d) / 2
    regional_interval <- gintervals.force_range(data.frame(
        chrom1 = interval$chrom1,
        start1 = interval$start1 - expand, end1 = interval$end1 + expand,
        chrom2 = interval$chrom2, start2 = interval$start2 - expand,
        end2 = interval$end2 + expand
    ))
    message(paste("combining multi tracks"))
    obs <- .shaman_combine_points_multi_tracks(obs_track_nms, regional_interval, min_dist)
    message(paste("replicating multi-contacts"))
    obs <- plyr::ddply(obs, c("contacts"), function(x) {
        return(x[rep(seq_len(nrow(x)), each = x$contact[1]), ])
    })
    if (nrow(obs) < 1000) {
        message(paste("insufficient data found in intervals: obs=", nrow(obs)))
        return(NULL)
    }

    # shuffle contacts in expanded interval
    message(paste0("shuffling ", nrow(obs), " observed points"))
    shuf_fn <- paste0(work_dir, "/", interval$chrom1, "_", interval$start1, "_", interval$start2, "_", expand, ".shuffled")

    if (is.na(decay_smooth)) {
        decay_smooth <- min(floor(dist_resolution / 10), 20)
    }
    samples_per_proposal_correction <- floor(nrow(obs) / 20)
    used_seed <- shaman_hic_matrix_shuffler_cpp(
        t(obs[, c("start1", "start2")]), shuf_fn, shuffle, 1, 0.5, dist_resolution, decay_smooth, 5, 0.25,
        max(obs$start2 - obs$start1), 1024, 1, grid_small, grid_high, grid_increase, grid_step_iter, 1, 1, seed
    )


    exp <- as.data.frame(data.table::fread(shuf_fn, header = TRUE))

    ret <- shaman_kk_norm(obs, exp, points, k = k, k_exp = 2 * k)
    ret$exp_fn <- shuf_fn
    ret$seed <- used_seed
    return(ret)
}

##########################################################################################################
#'  Compute a score matrix for observed data based on the expected for a given set of points
#'
#' \code{shaman_kk_norm}
#'
#' This function receives observed and expected data and compute the score on a given set of points.
#' The score for a point is the KS D-statistic of the distances to the points k-nearest-neighbors
#  in the observed data compared the the expected data.
#'
#' @param obs Dataframe containing the observed points
#' @param exp Dataframe containing the expected (shuffled) points.
#' @param points A score will be computed for each of the points.
#' @param k The number of neighbor distances used for the score on observed data.
#' For higher resolution maps, increase k. For lower resolution maps, decrease k.
#' @param k_exp The number of neighbor distances used for the score on expected data. Should reflect
#' the ratio between the total number of observed and expected over the entire chromosome.
#'
#' @return NULL if insufficient observed data, otherwise resturns a list containing 3 elements:
#' 1) points - start1, start2 and score for all observed points.
#' 2) obs - the observed points.
#' 3) exp - the expected points.
#'
#' @examples
#'
#' # Set misha db to test
#' library(misha)
#' gsetroot(shaman_get_test_track_db())
#' focus <- gintervals.2d(2, 176.6e06, 176.9e06, 2, 176.6e06, 176.9e06)
#' regional <- gintervals.2d(2, 176.5e06, 177e06, 2, 176.5e06, 177e06)
#' points <- gextract("hic_obs", focus, band = c(-5e05, -1024))
#' obs <- gextract("hic_obs", regional, band = c(-5e05, -1024))
#' exp <- gextract("hic_exp", regional, band = c(-5e05, -1024))
#' mat_score <- shaman_kk_norm(obs, exp, points, k = 100, k_exp = 200)
#' shaman_gplot_map_score(mat_score$points)
#' @export
##########################################################################################################
shaman_kk_norm <- function(obs, exp, points, k = 100, k_exp = 100) {
    message(paste0("going into knn witn ", nrow(obs), " observed and ", nrow(exp), " expected"))
    o_knn <- RANN::nn2(obs[, c("start1", "start2")], points[, c("start1", "start2")], k = k)
    message("going into shuffled knn")
    e_knn <- RANN::nn2(exp[, c("start1", "start2")], points[, c("start1", "start2")], k = k_exp)

    s_ks <- shaman_merge_ks_cpp(round(o_knn$nn.dist), round(e_knn$nn.dist))

    points$score <- 100 * ifelse(-s_ks$V1 < s_ks$V2, s_ks$V2, s_ks$V1)

    return(list(points = points, obs = obs, exp = exp))
}

##########################################################################################################
#  .shaman_combine_points_multi_tracks
#
##########################################################################################################
.shaman_combine_points_multi_tracks <- function(tracks, interval, min_dist) {
    points <- plyr::adply(tracks, 1, function(x) {
        p <- .shaman_cached_gextract(x, interval)
        if (isFALSE(p)) {
            p <- gextract(x, interval, colnames = c("contacts"))
        }
        return(p[abs(p$start1 - p$start2) > min_dist, ])
    })
    return(points[, -1])
}

# Extractions of whole rows of matrices, kept while shaman_score_hic_mat_for_track() scores the row
.shaman_extract_cache <- new.env(parent = emptyenv())

# TRUE when interval and outer are single 2D intervals and interval has integer coordinates and lies
# inside outer
.shaman_interval_within <- function(interval, outer) {
    if (NROW(interval) != 1 || NROW(outer) != 1) {
        return(FALSE)
    }
    co <- c(interval$start1, interval$end1, interval$start2, interval$end2)
    all(co == round(co)) && as.character(interval$chrom1) == as.character(outer$chrom1) &&
        as.character(interval$chrom2) == as.character(outer$chrom2) &&
        interval$start1 >= outer$start1 && interval$end1 <= outer$end1 &&
        interval$start2 >= outer$start2 && interval$end2 <= outer$end2
}

# The contacts of a points track in interval, from its contacts in a larger interval (p, rows as
# gextract() returns them: [x, x + 1) x [y, y + 1)). gextract() of a 2D track visits every object of
# the chromosome pair in a fixed order and returns those inside a single scope interval, so these are
# the same rows, in the same order, as gextract() of interval gives.
.shaman_points_in <- function(p, interval) {
    p <- p[p$start1 >= interval$start1 & p$start1 < interval$end1 & p$start2 >= interval$start2 & p$start2 < interval$end2, ]
    rownames(p) <- NULL
    p
}

.shaman_cache_extractions <- function(tracks, interval) {
    for (x in tracks) {
        if (gtrack.info(x)$type == "points") {
            p <- gextract(x, interval, colnames = c("contacts"))
            if (is.null(p)) {
                p <- data.frame(start1 = numeric(0), end1 = numeric(0), start2 = numeric(0), end2 = numeric(0))
            }
            assign(x, list(interval = interval, points = p), envir = .shaman_extract_cache)
        }
    }
}

# gextract(track, interval, colnames = "contacts") taken from a cached extraction of a larger
# interval (see .shaman_points_in), or FALSE if there is none
.shaman_cached_gextract <- function(track, interval) {
    ce <- .shaman_extract_cache[[track]]
    if (is.null(ce) || !.shaman_interval_within(interval, ce$interval)) {
        return(FALSE)
    }
    p <- .shaman_points_in(ce$points, interval)
    # gextract() returns NULL, not an empty frame, when there are no rows
    if (nrow(p) == 0) {
        return(NULL)
    }
    p
}

.shaman_compute_marginal_multi_tracks <- function(tracks, interval, min_dist) {
    total <- plyr::adply(tracks, 1, function(x) {
        gvtrack.create("v_sum", x, "weighted.sum")
        return(sum(gextract("v_sum", interval, iterator = interval, band = c(-max(gintervals.all()$end), -min_dist))$v_sum))
    })[, -1]
    return(sum(total))
}





# Runs commands on SGE via misha::gcluster.run. Commands are given either as expressions in '...'
# or as strings in 'command.list'. The call is evaluated in the caller's frame, so the jobs get
# the caller's variables (e.g. track_db, work_dir), as with a direct gcluster.run call.
.gcluster.run2 <- function(..., command.list = NULL, opt.flags = "", max.jobs = 400, debug = FALSE, R = "R") {
    if (!is.null(command.list)) {
        commands <- lapply(command.list, str2lang)
    } else {
        commands <- as.list(substitute(list(...))[-1L])
    }
    do.call(misha::gcluster.run,
        c(commands, list(opt.flags = opt.flags, max.jobs = max.jobs, debug = debug, R = R)),
        envir = parent.frame()
    )
}
