#'  generate an expected hic track based on observed hic data
#'
#' \code{shaman_shuffle_hic_track}
#'
#' This function generates an expected 2D hic track based on observed hic data.
#' Each chromosome is shuffled seperately, to generate an expected shuffled contact matrix
#' Note that this function requires sge (qsub) or multicore to be enables.
#' Parameter can be set via shaman.sge_support or shaman.mc_support in shaman.conf file.
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
#'
#' @export
##########################################################################################################
shaman_shuffle_hic_track <- function(track_db, obs_track_nm, work_dir,
    exp_track_nm=paste0(obs_track_nm, "_shuffle"), max_jobs=25,
    shuffle=80, grid_small = 500000, grid_high=1000000, grid_step_iter = 40,
    dist_resolution=NA, smooth=NA)
{
    gsetroot(track_db)
    # check tracks
    if (!gtrack.exists(obs_track_nm)) {
	stop(paste("Missing obs_track_nm (", obs_track_nm, ") in track db"))
    }
   if (gtrack.exists(exp_track_nm)) {
	stop(paste("exp_track_nm (", exp_track_nm, ") already exists in track db"))
   }
   #check sge
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

    if (sge_support) {
    commands <- paste0("{library(shaman); shaman_shuffle_hic_mat_for_track(\"", track_db, "\",\"", obs_track_nm, "\",\"",
	      work_dir, "\", \"", intervals$chrom, "\", ", intervals$start, ", ",
        intervals$end, ", ", intervals$start, ", ", intervals$end,  
        ", min_dist=1024, dist_resolution=", dist_resolution, ", decay_smooth=",
        smooth, ", shuffle=", shuffle, ", grid_small=", grid_small, ", grid_high=", grid_high,
        ", grid_step_iter=", grid_step_iter,
        ", raw_ext=\"full_chrom_raw\", shuffled_ext=\"full_chrom_shuffled\", sort_uniq=TRUE)}")
     res <- .gcluster.run2(command.list = commands, opt.flags=sge_flags, max.jobs=max_jobs)
    } else {
	doMC::registerDoMC(cores=max_jobs)
	res=plyr::ddply(intervals, .(chrom, start), function(x) {
		shaman_shuffle_hic_mat_for_track(track_db, obs_track_nm, work_dir, x$chrom[1],
			x$start[1], x$end[1], x$start[1], x$end[1], min_dist=1024,
			dist_resolution=dist_resolution, decay_smooth=smooth, shuffle=shuffle,
			grid_small=grid_small, grid_high=grid_high, grid_step_iter=grid_step_iter,
			raw_ext="full_chrom_raw", shuffled_ext="full_chrom_shuffled", sort_uniq=TRUE)
	}, .parallel=TRUE)
    }
    exp_shuf_files <- paste0(obs_track_nm, "_", intervals$chrom, "_0_0.full_chrom_shuffled.uniq")
    obs_shuf_files <- list.files(work_dir, pattern=paste0(obs_track_nm, ".*full_chrom_shuffled.uniq"))
    missing_files <- exp_shuf_files[!exp_shuf_files %in% obs_shuf_files]
    if (length(missing_files) > 0) {
        warning(paste(length(missing_files), "full chrom files were not shuffled:\n", paste(missing_files, collapse=",")))
    }

    # reaching this point means that all matrices should be in work_dir --> creating track this section
    # should be replaced after misha generates a proper import for overlapping contacts
    files <-  c()
    for (chrom in intervals$chrom) {
        fn <- paste0(work_dir, obs_track_nm, "_", chrom, "_0_0.full_chrom_shuffled.uniq")
        message(chrom)
        if (file.exists(fn)) {
            files <- c(files, fn)
        } else {
            message("missing file")
        }
    }
    gtrack.2d.import(exp_track_nm, paste("shuffled 2d track with shuffle factor =",
      shuffle, ", based on", obs_track_nm), files)

    # cleanup work dir from all temporary files
    try(system(sprintf("rm %s%s*", work_dir, obs_track_nm)))
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
#'
#' @export
##########################################################################################################
shaman_shuffle_hic_mat_for_track <- function(track_db, track, work_dir, chrom, start1, end1, start2, end2,
  min_dist=1024, max_dist=max(gintervals.all()$end), dist_resolution=NA,
  decay_smooth=NA, proposal_iterations=1e+07, shuffle=80, hic_mcmc_max_resolution=400,
  raw_ext="raw", shuffled_ext="shuffled", grid_small = 500000, grid_high=1000000, grid_increase=500000,
  grid_step_iter = 40, sort_uniq=FALSE) {

  raw_fn <- paste0(work_dir, "/", track, "_", chrom, "_", start1, "_", start2, ".", raw_ext)
  shuf_fn <- paste0(work_dir, "/", track, "_", chrom, "_", start1, "_", start2, ".", shuffled_ext)
  if (!file.exists(shuf_fn)) {
      x <- sample(1:10, 1)
      system(paste("sleep", x))
      options(gmultitasking=FALSE)
      options(gmax.data.size=1e+09)
      gsetroot(track_db)
      scope <- gintervals.2d(chrom, start1, end1, chrom, start2, end2)
      a <- gextract(track, scope, band=c(-max_dist, -min_dist + 1), colnames="contact")
      if (is.null(nrow(a))) {
          message("not shuffling, no data")
          system(paste("perl -e'print\"start1\tstart2\n\";' > ", shuf_fn))
          return(0)
      }
      # multiply each line according to the number of observed counts
      a <- plyr::ddply(a, c("contact"), function(x) {
          return(x[rep(seq_len(nrow(x)), each=x$contact[1]), c("start1", "start2")])
      })[, 2:3]

      # automatic decision on resolution, smoothing and samples per correction
      if (is.na(dist_resolution)) {
          dist_resolution=min(floor((nrow(a)/200)/(log2(max_dist) - log2(min_dist))), hic_mcmc_max_resolution)
      }

      message(paste("writing", nrow(a), "raw data, dist resolution=", dist_resolution))
      if (nrow(a) < 2000 | dist_resolution == 0) {
          message("not shuffling, leaving raw")
          data.table::fwrite(format(rbind(a,setNames(rev(a), names(a))), scientific=FALSE), shuf_fn, quote=FALSE, row.names=F,
              sep="\t")
      } else {
          if (is.na(decay_smooth)) {
              decay_smooth=min(floor(dist_resolution/10), 20)
          }

          ret = shaman_hic_matrix_shuffler_cpp(t(a[, c("start1", "start2")]), 
		shuf_fn, shuffle, 1, 0.5, dist_resolution, decay_smooth, 5, 0.25,
		max_dist, 1024, 1, grid_small, grid_high, grid_increase, grid_step_iter, 0, 1)
      }
  }
  if (sort_uniq) {
      ret=1
      system(sprintf("echo -e \"chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tobs\" > %s.uniq", shuf_fn))
      system(sprintf("cat %s | grep -v start | sort | uniq -c | awk '{ print \"%s\" \"\t\" $2 \"\t\" ($2+1) \"\t\" \"%s\" \"\t\" $3 \"\t\" ($3+1) \"\t\" $1}' >> %s.uniq",
          shuf_fn, chrom, chrom, shuf_fn))
  }
  return(ret)
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
#'
#' Each step creates temporary files of the matrix scores which are then joined to a track.
#' Temporary files are deleted upon track creation.
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
#'
#' @export
##########################################################################################################
shaman_score_hic_track <- function(track_db, work_dir, score_track_nm, obs_track_nms,
  exp_track_nms=paste0(obs_track_nms, "_shuffle"), points_track_nms=obs_track_nms,
  near_cis=5e06, expand=2e06, k=100, max_jobs=100)
{
    gsetroot(track_db)
    # check tracks
    if (sum(gtrack.exists(obs_track_nms)) < length(obs_track_nms)) {
	stop(paste("Missing obs_track_nm (", obs_track_nm[!gtrack.exists(obs_track_nms)], ") in track db"))
    }
    if (sum(gtrack.exists(exp_track_nms)) < length(exp_track_nms)) {
	stop(paste("Missing exp_track_nm (", exp_track_nm[!gtrack.exists(exp_track_nms)], ") in track db"))
    }
   if (gtrack.exists(score_track_nm)) {
	stop(paste("score_track_nm (", score_track_nm, ") already exists in track db"))
   }
   #check sge
   .shaman_check_config("shaman.sge_support")
   sge_support <- getOption("shaman.sge_support")
   mc_support <- getOption("shaman.mc_support")
   if (!sge_support & !mc_support) {
     stop(paste("score_hic_track function requires SGE or multicore support. If available, set configuration parameter shaman.sge_support or shaman.mc_support to 1"))
   }
   sge_flags <- getOption("shaman.sge_flags")

  #split all chromosomes to rectangles of size near_cis*near_cis
  band <- c(-max(gintervals.all()$end), 0)
  near_cis_intervals <- gintervals.force_range(plyr::ddply(gintervals.all(), c("chrom"), function(x) {
    data.frame(chrom=rep(x$chrom[1], nrow(x)),
    start=seq(0, floor(x$end/near_cis)*near_cis, by=near_cis),
    end=seq(near_cis, ceiling(x$end/near_cis)*near_cis, by=near_cis))}))
  g = expand.grid(1:nrow(near_cis_intervals), 1:nrow(near_cis_intervals))
  g = g[near_cis_intervals$chrom[g$Var1] == near_cis_intervals$chrom[g$Var2],]
  near_cis_2d = gintervals.2d(
        near_cis_intervals$chrom[g$Var1], near_cis_intervals$start[g$Var1], near_cis_intervals$end[g$Var1],
        near_cis_intervals$chrom[g$Var2], near_cis_intervals$start[g$Var2], near_cis_intervals$end[g$Var2])

  #selecting only the upper half of the matrix
  near_cis_2d_upper_mat = gintervals.2d.band_intersect(near_cis_2d, band)
  expected_files = paste0(work_dir , "/", paste0(obs_track_nms, collapse=".") , ".",
        near_cis_2d_upper_mat$chrom1, ".", near_cis_2d_upper_mat$start1, ".", near_cis_2d_upper_mat$start2, ".score")

   existing_files = file.exists(expected_files)
   missing_files = expected_files[!existing_files]
   message(paste("missing", length(missing_files), "score files"))
   near_cis_2d_upper_mat = near_cis_2d_upper_mat[!existing_files,]
   expected_files = paste0(work_dir , "/", paste0(obs_track_nms, collapse=".") , ".",
        near_cis_2d_upper_mat$chrom1, ".", near_cis_2d_upper_mat$start1, ".", near_cis_2d_upper_mat$start2, ".score")
  while (nrow(near_cis_2d_upper_mat)>0) {
    #compute scores for each of the small matrices
    if (sge_support) {
      commands = paste0("{library(shaman); shaman_score_hic_mat_for_track(track_db, work_dir, obs_track_nms, exp_track_nms, points_track_nms, \"",
        near_cis_2d_upper_mat$chrom1, "\", ", near_cis_2d_upper_mat$start1, ", ",
        near_cis_2d_upper_mat$end1, ",", near_cis_2d_upper_mat$start2, ", ",
        near_cis_2d_upper_mat$end2, ", ", expand, ", ", k, ")}")
    #commands <- paste(commands, collapse=",")

     res <- .gcluster.run2(command.list = commands, opt.flags=sge_flags, max.jobs=max_jobs)
   } else {
     doMC::registerDoMC(cores=max_jobs)
     res=plyr::ddply(near_cis_2d_upper_mat, .(chrom1, start1, start2), function(x) {
        shaman_score_hic_mat_for_track(track_db, work_dir, obs_track_nms, exp_track_nms, points_track_nms,
        x$chrom1[1], x$start1[1],x$end1[1], x$start2[1], x$end2[1], expand,  k)},
	.parallel=TRUE)
   }
    #res <- eval(parse(text=paste("gcluster.run(", commands, ",opt.flags=\"", sge_flags,  "\" ,max.jobs=", max_jobs, ")")))
    #check to see if there are any missing files
    existing_files = file.exists(expected_files)
    missing_files = expected_files[!existing_files]
    message(paste("missing", length(missing_files), "score files"))
    near_cis_2d_upper_mat = near_cis_2d_upper_mat[!existing_files,]
    expected_files = paste0(work_dir , "/", paste0(obs_track_nms, collapse=".") , ".",
        near_cis_2d_upper_mat$chrom1, ".", near_cis_2d_upper_mat$start1, ".", near_cis_2d_upper_mat$start2, ".score")
  }
  #reaching this point means that all points have a score - need to create track
  score_files = list.files(work_dir, paste0(obs_track_nms, ".*.score$"), full.names=TRUE)

  gtrack.2d.import_contacts(score_track_nm, paste("normalized score of 2d track with k =", k, "based on",
	paste(obs_track_nms, collapse=", "), "and",
	paste(exp_track_nms, collapse=", ")), score_files, allow.duplicates=FALSE)

  #cleanup work dir from all temporary files
  for (track in obs_track_nms) {
    try(system(sprintf('rm %s%s*', work_dir, track)))
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
#' @export
##########################################################################################################
shaman_score_hic_mat_for_track <- function(track_db, work_dir, obs_track_nms, exp_track_nms, points_track_nms,
  chrom, start1, end1, start2, end2, expand=2e06, k=100, min_dist=1024)
{
  fn <- paste0(work_dir , "/", paste0(obs_track_nms, collapse=".") , ".", chrom, ".", start1, ".", start2, ".score")
  if(file.exists(fn)) {
     return(0)
  }
  options(gmax.data.size=1e09)
  regional_interval <- gintervals.force_range(data.frame(chrom1=chrom, start1 = start1-expand, end1 = end1+expand,
    chrom2=chrom, start2=start2-expand, end2 = end2+expand))
  focus_interval = gintervals.2d(chrom, start1, end1, chrom, start2, end2)
  n = shaman_score_hic_mat(obs_track_nms, exp_track_nms, focus_interval, regional_interval, points_track_nms=points_track_nms, min_dist=min_dist, k=k)
  if (is.null(n)) {
    system(paste("perl -e'print\"chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tscore\";' > ", fn))
    #insufficient data in region - not writing region file
        return(0);
  }
  if (length(n) < 1) {
       #knn did not complete
      return(-1);
  }
  p = n$points
 data.table::fwrite(format(p[p$start1<=p$start2,
 	c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "score")],
       scientific=FALSE),
       fn, row.names=FALSE, quote=FALSE, sep="\t")
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
#' @export
##########################################################################################################
shaman_score_hic_mat <- function(obs_track_nms, exp_track_nms, focus_interval, regional_interval,
  points_track_nms=obs_track_nms, min_dist=1024, k=100, k_exp=2*k)
{
  points <- .shaman_combine_points_multi_tracks(points_track_nms, focus_interval, min_dist)
  if (is.null(points)) {
    message("number of points in focus interval = 0")
    return(NULL)
  }
  if (nrow(points) < 1000) {
    message("number of points in focus interval < 1000")
    return(NULL)
  }

  points = unique(points[,c("chrom1", "start1", "end1", "chrom2", "start2", "end2")])
  message(paste0("kk norm on ", nrow(points), " points"))
  return(shaman_score_hic_points(obs_track_nms, exp_track_nms, points, regional_interval, min_dist, k=k, k_exp=k_exp))
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
#'
#' @return NULL if insufficient observed data, otherwise resturns a list containing 3 elements:
#' 1) points - start1, start2 and score for all observed points.
#' 2) obs - the observed points.
#' 3) exp - the expected points.
#'
#' @export
##########################################################################################################
shaman_score_hic_points <- function(obs_track_nms, exp_track_nms, points, regional_interval, min_dist=1024, k=100, k_exp=2*k)
{
  knn_pl <- sprintf("%s/%s", system.file("perl", package="shaman"), getOption("shaman.ks_pl"))
  o_knn.tmp <- tempfile("knn_o_")
  e_knn.tmp <- gsub("knn_o_", "knn_e_", o_knn.tmp)
  ks_knn.tmp <- gsub("knn_o_", "knn_ks_", o_knn.tmp)
  message(paste("obs = ", paste(obs_track_nms, collapse=",")))
  obs <- .shaman_combine_points_multi_tracks(obs_track_nms, regional_interval, min_dist)
  if (is.null(obs) | nrow(points)==0) {
     message(paste("0 data found in intervals, focus interval=", nrow(points)))
     return(NULL)
  }
  obs <- plyr::ddply(obs, c("contacts"), function(x) { return(x[rep(seq_len(nrow(x)), each=x$contact[1]), ]) })
  if (nrow(obs) < k) {
     message(paste("insufficient data found in intervals: obs=", nrow(obs)))
     return(NULL)
  }
  n_obs = nrow(obs)
  o_knn <- RANN::nn2(obs[,c("start1", "start2")], points[,c("start1", "start2")], k=k)
  message("write tab 1")
  data.table::fwrite(as.data.frame(round(o_knn$nn.dist)), o_knn.tmp, sep="\t", col.names=F, quote=F, row.names=F)
  if (!file.exists(o_knn.tmp)) {
    message(paste0("problem writing ", o_knn.tmp))
    return(0);
  }
  rm(o_knn)
  rm(obs)
  gc()

  exp <- .shaman_combine_points_multi_tracks(exp_track_nms, regional_interval, min_dist)
  if (is.null(exp)) {
      message(paste("0 data found in intervals: exp"))
      return(NULL)
  }
  exp <- plyr::ddply(exp, c("contacts"), function(x) { return(x[rep(seq_len(nrow(x)), each=x$contact[1]), ]) })
  if (nrow(exp) < k) {
      message(paste("insufficient data found in intervals: exp=", nrow(exp)))
      return(NULL)
  }
  n_exp = nrow(exp)
  if (is.na(k_exp)) {
	#computing k_exp by number of points
	k_exp = round(k * n_exp / n_obs)
  }
  message(paste0("n_obs = ", n_obs, ", n_exp = ", n_exp, ", k_exp = ", k_exp))
  e_knn <- RANN::nn2(exp[,c("start1", "start2")], points[,c("start1", "start2")], k=k_exp)
  message("write tab 2")
  data.table::fwrite(as.data.frame(round(e_knn$nn.dist)), e_knn.tmp, sep="\t", col.names=F, quote=F, row.names=F)
  if (!file.exists(e_knn.tmp)) {
    message(paste0("problem writing ", e_knn.tmp))
    return(0);
  }
  rm(e_knn)
  rm(exp)
  gc()
  system(sprintf("perl %s %s %s >%s", knn_pl, o_knn.tmp, e_knn.tmp, ks_knn.tmp) )
  s_ks <- as.data.frame(data.table::fread(ks_knn.tmp))

  points$score <- 100 *  ifelse(-s_ks$V1 < s_ks$V2, s_ks$V2,  s_ks$V1)
  try(system(sprintf('rm %s %s %s', o_knn.tmp, e_knn.tmp, ks_knn.tmp)))

  return(list(points=points))
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

#' @return NULL if insufficient observed data, otherwise resturns a list containing 3 elements:
#' 1) points - start1, start2 and score for all observed points.
#' 2) obs - the observed points.
#' 3) exp - the expected points.
#' 4) obs_fn - the name of the observed data file
#' 4) exp_fn - the name of the expected (shuffled) data file
#'
#' @examples
#' init_shaman_examples...
#' shaman_shuffle_and_score_hic_mat(...)
#' @export
##########################################################################################################

shaman_shuffle_and_score_hic_mat <- function(obs_track_nms, interval, work_dir, expand=1e06, min_dist=1024, k=100,
  dist_resolution=NA, decay_smooth=NA, hic_mcmc_max_resolution=400, shuffle=80,
  grid_small = 500000, grid_high=1000000, grid_increase=500000,
  grid_step_iter = 40)
{

  if (interval$chrom1 != interval$chrom2) {
	stop("Only cis intervals supported")
  }
  points <-  .shaman_combine_points_multi_tracks(obs_track_nms, interval, min_dist)
  points <- points[points$start2 > points$start1,]
  if (is.null(points)) {
    message("number of points in interval = 0")
    return(NULL)
 }
  message(paste("working on ", nrow(points), "points"))
  if (is.na(dist_resolution)) {
        dist_resolution = min(floor((nrow(points)/200) / (log2(max(points$start2-points$start1)) -
	log2(1024))), hic_mcmc_max_resolution);
  }
  # adjust expand such that first and last dist bins are complete
  message(paste("adjusting expand"))
  max_d <- interval$end2 - interval$start1
  expand <- (2 ** (ceiling(log2(max_d + 2 * expand) * dist_resolution) / dist_resolution) - max_d ) / 2
  regional_interval = gintervals.force_range(data.frame(chrom1=interval$chrom1,
        start1=interval$start1 - expand, end1 =interval$end1 + expand,
        chrom2=interval$chrom2, start2=interval$start2 - expand,
        end2=interval$end2 + expand))
  message(paste("combining multi tracks"))
  obs <-  .shaman_combine_points_multi_tracks(obs_track_nms, regional_interval, min_dist)
  message(paste("replicating multi-contacts"))
  obs = plyr::ddply(obs, c("contacts"), function(x) { return(x[rep(seq_len(nrow(x)), each=x$contact[1]), ]) })
  if (nrow(obs) < 1000) {
    message(paste("insufficient data found in intervals: obs=", nrow(obs)))
    return(NULL)
  }

  # shuffle contacts in expanded interval
  message(paste0("shuffling ", nrow(obs), " observed points"))
  shuf_fn = paste0(work_dir, "/",  interval$chrom1, "_", interval$start1, "_", interval$start2, "_", expand, ".shuffled")

  if (is.na(decay_smooth)) {
        decay_smooth = min(floor(dist_resolution / 10),20)
  }
  samples_per_proposal_correction = floor(nrow(obs)/20)
  shaman_hic_matrix_shuffler_cpp(t(obs[, c("start1", "start2")]), shuf_fn, shuffle, 1, 0.5, dist_resolution, decay_smooth, 5, 0.25,
	max(obs$start2-obs$start1), 1024, 1, grid_small, grid_high, grid_increase, grid_step_iter, 1, 1)


  exp = as.data.frame(data.table::fread(shuf_fn, header=T))

  ret = shaman_kk_norm(obs, exp, points, k=k, k_exp=2*k)
  ret$obs_fn = raw_fn
  ret$exp_fn = shuf_fn
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
#' @export
##########################################################################################################
shaman_kk_norm <- function(obs, exp, points, k=100, k_exp=100)
{
  .shaman_check_config("shaman.ks_pl")
  knn_pl <- sprintf("%s/%s", system.file("perl", package="shaman"), getOption("shaman.ks_pl"))
  message(paste0("going into knn witn ", nrow(obs), " observed and ", nrow(exp), " expected"))
  o_knn <- RANN::nn2(obs[,c("start1", "start2")], points[,c("start1", "start2")], k=k)
  message("going into shuffled knn")
  e_knn <- RANN::nn2(exp[,c("start1", "start2")], points[,c("start1", "start2")], k=k_exp)

  o_knn.tmp <- tempfile("knn_o_")
  e_knn.tmp <- gsub("knn_o_", "knn_e_", o_knn.tmp)
  ks_knn.tmp <- gsub("knn_o_", "knn_ks_", o_knn.tmp)
  message("write tab 1")
  data.table::fwrite(as.data.frame(round(o_knn$nn.dist)), o_knn.tmp, sep="\t", col.names=F, quote=F, row.names=F)
  if (!file.exists(o_knn.tmp)) {
    message(paste0("problem writing ", o_knn.tmp))
    return(0);
  }
  message("write tab 2")
  data.table::fwrite(as.data.frame(round(e_knn$nn.dist)), e_knn.tmp, sep="\t", col.names=F, quote=F, row.names=F)
  if (!file.exists(e_knn.tmp)) {
    message(paste0("problem writing ", e_knn.tmp))
    return(0);
  }

  system(sprintf("perl %s %s %s >%s", knn_pl, o_knn.tmp, e_knn.tmp, ks_knn.tmp) )

  s_ks <- as.data.frame(data.table::fread(ks_knn.tmp))

  points$score <- 100 *  ifelse(-s_ks$V1 < s_ks$V2, s_ks$V2,  s_ks$V1)
  try(system(sprintf('rm %s %s %s', o_knn.tmp, e_knn.tmp, ks_knn.tmp)))

  return(list(points=points, obs=obs, exp=exp))
}

##########################################################################################################
#  .shaman_combine_points_multi_tracks
#
##########################################################################################################
.shaman_combine_points_multi_tracks <- function(tracks, interval, min_dist)
{
  points <- plyr::adply(tracks, 1, function(x) {
    p <- gextract(x, interval, colnames=c("contacts"))
     return(p[abs(p$start1-p$start2) > min_dist,]) })
  return(points[,-1])
}

.shaman_compute_marginal_multi_tracks <- function(tracks, interval, min_dist)
{
   total = plyr::adply(tracks, 1, function(x) {
     gvtrack.create("v_sum", x, "weighted.sum")
     return(sum(gextract("v_sum", interval, iterator=interval, band=c(-max(gintervals.all()$end), -min_dist))$v_sum))
     })[,-1]
   return(sum(total))
}





.gcluster.run2 <- function (..., command.list = NULL, opt.flags = "", max.jobs = 400, debug = FALSE,  R = "R")
{
    if (!is.null(command.list)){
        commands <- plyr::llply(command.list, function(x) parse(text=x))
    } else {
        commands <- as.list(substitute(list(...))[-1L])
    }

    if (length(commands) < 1)
        stop("Usage: gculster.run(..., command.list=NULL, opt.flags = \"\" max.jobs = 400, debug = FALSE)",
            call. = F)
    if (!length(system("which qsub", ignore.stderr = T, intern = T)))
        stop("gcluster.run must run on a host that supports Sun Grid Engine (qsub)",
            call. = F)
    .gcheckroot()
    tmp.dirname <- ""
    submitted.jobs <- c()
    tryCatch({
        tmp.dirname <- tempfile(pattern = "", tmpdir = paste(get("GROOT"),
            "/tmp", sep = ""))
        if (!dir.create(tmp.dirname, recursive = T, mode = "0777"))
            stop(sprintf("Failed to create a directory %s", tmp.dirname),
                call. = F)
        cat("Preparing for distribution...\n")
        save(.GLIBDIR, file = paste(tmp.dirname, "libdir", sep = "/"))
        vars <- ls(all.names = TRUE, envir = parent.frame())
        envir <- parent.frame()
        while (!identical(envir, .GlobalEnv)) {
            envir <- parent.env(envir)
            vars <- union(vars, ls(all.names = TRUE, envir = envir))
        }
        save(list = vars, file = paste(tmp.dirname, "envir",
            sep = "/"), envir = parent.frame())
        .GSGECMD <- commands
        save(.GSGECMD, file = paste(tmp.dirname, "commands",
            sep = "/"))
        opts <- options()
        save(opts, file = paste(tmp.dirname, "opts", sep = "/"))
        cat("Running the commands...\n")
        completed.jobs <- c()
        progress <- -1
        repeat {
            num.running.jobs <- length(submitted.jobs) - length(completed.jobs)
            if (length(submitted.jobs) < length(commands) &&
                num.running.jobs < max.jobs) {
                istart <- length(submitted.jobs) + 1
                iend <- min(length(commands), istart + (max.jobs -
                  num.running.jobs) - 1)
                for (i in istart:iend) {
                  out.file <- sprintf("%s/%d.out", tmp.dirname,
                    i)
                  err.file <- sprintf("%s/%d.err", tmp.dirname,
                    i)
                  script <- paste(get(".GLIBDIR"), "exec", "sgjob.sh",
                    sep = "/")
                  command <- sprintf("unset module; qsub -terse -S /bin/bash -o %s -e %s -V %s %s %d '%s' '%s'",
                    out.file, err.file, opt.flags, script, i,
                    tmp.dirname, R)
                  jobid <- system(command, intern = TRUE)
                  if (length(jobid) != 1)
                    stop("Failed to run qsub", call. = FALSE)
                  if (debug)
                    cat(sprintf("\tSubmitted job %d (id: %s)\n",
                      i, jobid))
                  submitted.jobs <- c(submitted.jobs, jobid)
                }
            }
            Sys.sleep(3)
            running.jobs <- .gcluster.running.jobs(submitted.jobs)
            old.completed.jobs <- completed.jobs
            completed.jobs <- setdiff(submitted.jobs, running.jobs)
            if (debug) {
                delta.jobs <- setdiff(completed.jobs, old.completed.jobs)
                if (length(delta.jobs) > 0) {
                  for (jobid in delta.jobs) cat(sprintf("\tJob %d (id: %s) completed\n",
                    match(jobid, submitted.jobs), jobid))
                }
                if (!length(running.jobs) && length(submitted.jobs) ==
                  length(commands))
                  break
                new.progress <- length(completed.jobs)
                if (new.progress != progress) {
                  progress <- new.progress
                  cat(sprintf("\t%d job(s) still in progress\n",
                    length(commands) - progress))
                }
            }
            else {
                if (!length(running.jobs) && length(submitted.jobs) ==
                  length(commands))
                  break
                new.progress <- as.integer(100 * length(completed.jobs)/length(commands))
                if (new.progress != progress) {
                  progress <- new.progress
                  cat(sprintf("%d%%...", progress))
                }
                else cat(".")
            }
        }
        if (!debug && progress != -1 && progress != 100)
            cat("100%\n")
    }, interrupt = function(interrupt) {
        cat("\n")
        stop("Command interrupted!", call. = FALSE)
    }, finally = {
        cleanup.finished <- FALSE
        while (!cleanup.finished) {
            tryCatch({
                if (length(submitted.jobs) > 0) {
                  running.jobs <- .gcluster.running.jobs(submitted.jobs)
                  answer <- c()
                  for (i in 1:length(commands)) {
                    res <- list()
                    res$exit.status <- NA
                    res$retv <- NA
                    res$stdout <- NA
                    res$stderr <- NA
                    if (submitted.jobs[i] %in% running.jobs)
                      res$exit.status <- "interrupted"
                    else {
                      fname <- sprintf("%s/%d.retv", tmp.dirname,
                        i)
                      if (file.exists(fname)) {
                        load(fname)
                        res$exit.status <- "success"
                        res$retv <- retv
                      }
                      else res$exit.status <- "failure"
                    }
                    out.file <- sprintf("%s/%d.out", tmp.dirname,
                      i)
                    if (file.exists(out.file)) {
                      f <- file(out.file, "rc")
                      res$stdout <- readChar(f, 1000000)
                      close(f)
                    }
                    err.file <- sprintf("%s/%d.err", tmp.dirname,
                      i)
                    if (file.exists(err.file)) {
                      f <- file(err.file, "rc")
                      res$stderr <- readChar(f, 1000000)
                      close(f)
                    }
                    answer[[i]] <- res
                  }
                  for (job in running.jobs) system(sprintf("qdel %s",
                    job), ignore.stderr = T, intern = T)
                  unlink(tmp.dirname, recursive = TRUE)
                  return(answer)
                }
                unlink(tmp.dirname, recursive = TRUE)
                cleanup.finished <- TRUE
            }, interrupt = function(interrupt) {
            })
        }
    })
}


