% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shaman_feature_grid.R
\name{shaman_plot_feature_grid}
\alias{shaman_plot_feature_grid}
\title{Plot spatial enrichment of contacts between two features}
\usage{
shaman_plot_feature_grid(
  grids,
  range,
  grid_resolution,
  plot_resolution,
  type = "enrichment",
  fig_fn = "",
  fig_width = 400,
  fig_height = 400,
  pal = colorRampPalette(c("blue", "#87FFFF", "white", "#FF413D", "black")),
  zlim = c(-1, 1)
)
}
\arguments{
\item{grids}{List of grids generated for the same two feature sets from replicate observed / expected hic tracks.}

\item{range}{Window size around the grid point that was used to generate the grid.}

\item{grid_resolution}{Size of bins that was used to generate the grid.}

\item{plot_resolution}{Size of bins to display.}

\item{type}{The type of data to plot, can be "observed" for log10(obs/sum(obs)), "expected" for log10(exp/sum(exp)) or
"enrichment" for log2(obs/sum(obs) / exp/sum(exp))}

\item{fig_fn}{Name of png file to output to. Empty string will cause the figure to be plotted to the current device.}

\item{fig_width}{Width in pixels of output png figure.}

\item{fig_height}{Height in pixels of output png figure.}

\item{pal}{Color palette to use for image}

\item{zlim}{The minimum and maximum values for which colors should be plotted. Suggested zlim values by type:
enrichment: (-1,1), obs,exp: (-4.5, -3)}
}
\description{
\code{shaman_plot_feature_grid}
}
\details{
After contact feature grids are generated via shaman_generate_feature_grid, use this function to visualize them.
}
\examples{

# Set misha db to test
gsetroot(shaman_get_test_track_db())
grid <- shaman_generate_feature_grid(shaman::ctcf_forward, shaman::ctcf_reverse, "hic_obs",
    exp_track_nm = "hic_exp", score_track_nm = "hic_score"
)
shaman_plot_feature_grid(list(grid), 25000, 500, 500)
shaman_plot_feature_grid(list(grid), 25000, 500, 1000)
shaman_plot_feature_grid(list(grid), 25000, 500, 2000)
}
