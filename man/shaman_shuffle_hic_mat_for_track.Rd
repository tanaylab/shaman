% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shaman.R
\name{shaman_shuffle_hic_mat_for_track}
\alias{shaman_shuffle_hic_mat_for_track}
\title{Generate an expected matrix from observed data as a process for generating an expected track}
\usage{
shaman_shuffle_hic_mat_for_track(
  track_db,
  track,
  work_dir,
  chrom,
  start1,
  end1,
  start2,
  end2,
  min_dist = 1024,
  max_dist = max(gintervals.all()$end),
  dist_resolution = NA,
  decay_smooth = NA,
  proposal_iterations = 10000000,
  shuffle = 80,
  hic_mcmc_max_resolution = 400,
  raw_ext = "raw",
  shuffled_ext = "shuffled",
  grid_small = 500000,
  grid_high = 1000000,
  grid_increase = 500000,
  grid_step_iter = 40,
  sort_uniq = FALSE
)
}
\arguments{
\item{track_db}{Directory of the misha database.}

\item{track}{Name of observed 2D genomic track for the hic data.}

\item{work_dir}{Centralized directory to store temporary files.}

\item{chrom}{The chormosome of the matrix.}

\item{start1}{The start coordinate of the first dimension.}

\item{end1}{The end coordinate of the first dimension.}

\item{start2}{The start coordinate of the second dimension.}

\item{end2}{The end coordinate of the second dimension.}

\item{min_dist}{The minimum distance between contact end points.}

\item{max_dist}{The maximum distance between contact end points.}

\item{dist_resolution}{Number of bins in each log2 distance unit. If NA, value is determined
based on observed data (recommended).}

\item{decay_smooth}{Number of bins to use for smoothing the MCMC target function: the decay curve.
If NA, value is determined based on observed data (recommended).}

\item{proposal_iterations}{Number of MCMC sampling iterations between proposal corrections.}

\item{shuffle}{Number of shuffling rounds for each observed point.}

\item{hic_mcmc_max_resolution}{Maximum number of bins for each log2 unit}

\item{raw_ext}{File extension of the observed data.}

\item{shuffled_ext}{File extension of the shuffled data.}

\item{grid_small}{Initial size of maximum distance between contact pairs consdered for switching}

\item{grid_high}{Final size of maximum distance between contact pairs consdered for switching}

\item{grid_increase}{Grid increase size}

\item{grid_step_iter}{Number of iterations in each grid size}

\item{sort_uniq}{Binary flag, indicating whether the shuffled matrix file should be sorted and
contacts combined. This is required prior to importing the track to misha, and should be applied
to full chromosomes only.}
}
\description{
\code{shuffle_hic_mat_for_track}
}
\details{
This function generates an expected 2D hic matrix from observed hic data. Should not be called externally,
The observed data is a combination of observed contacts in scope plus already shuffled
near-cis contacts (stored in work_dir) which we sample from to maintain the decay probability curve.

Each step creates temporary files of the shuffled matrices which are then joined to a track.
Temporary files are deleted upon track creation.
}
