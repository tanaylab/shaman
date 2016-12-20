#' plot an observed or expected hic map
#'
#' \code{shaman_gplot_map}
#'
#' Plots hic contact matrix.
#'
#' @param points A dataframe containing the points (start1, start2).
#' @param rotate Binary flag, indicating if the plot should be rotated by 45 degrees.
#' @param point_size Cex size of the points in the plot.
#' @param add_axis Binary flag, indicating if axis should be added to plot.
#' @return gplot containing the map
#' @export
##########################################################################################################
shaman_gplot_map <- function(points, interval_range=NA, rotate=TRUE, point_size=0.1, add_axis=TRUE) {
  if (!all(c("start1", "start2") %in% colnames(points))) {
    stop("points_score data frame must contain the following columns: start1, start2")
  }
  if (is.na(interval_range)[1]) {
        interval_range=data.frame(start=min(points$start1), end=max(points$start1))
  }
  if (rotate) {
    points = points[points$start2 > points$start1,]
    map_gplot <- ggplot2::ggplot(points,
      ggplot2::aes(x=(start1+start2)/2, y=(start2-start1)/2)) +
      ggplot2::coord_cartesian(xlim=c(interval_range$start, interval_range$end),
        ylim=c(0, (interval_range$end-interval_range$start)/2)) +
      ggplot2::theme (panel.grid.major = ggplot2::element_blank(),
                      panel.grid.minor = ggplot2::element_blank())
  } else {
    map_gplot <- ggplot2::ggplot(points,
      ggplot2::aes(x=start1, y=start2)) +
      ggplot2::scale_x_continuous(position = "top") +
      ggplot2::theme(axis.line.y = ggplot2::element_line(size=0.2))
  }
  map_gplot <- map_gplot +
    ggplot2::geom_point(size=point_size, alpha=0.1) +
    ggplot2::scale_x_continuous(expand=c(0,0)) +
    ggplot2::scale_y_continuous(expand=c(0,0)) +
    ggplot2::theme_bw() + ggplot2::theme( legend.position="none",
         panel.border = ggplot2::element_blank(),
         axis.title.x = ggplot2::element_blank(),
         axis.title.y = ggplot2::element_blank(),
	 axis.line.x = ggplot2::element_line(size=0.2))
  if (add_axis == FALSE) {
	map_gplot <- map_gplot + ggplot2::theme(
         axis.ticks = ggplot2::element_blank(),
         axis.text.x = ggplot2::element_blank(),
         axis.text.y = ggplot2::element_blank(),
	 axis.line.x = ggplot2::element_blank(),
	 axis.line.y = ggplot2::element_blank(),
         axis.ticks.length = ggplot2::unit(0,"null") )
  }
  return(map_gplot)
}

#' plot a normlized hic map
#'
#' \code{shaman_gplot_map_score}
#'
#' Plots observerved hic contact matrix color-coded by normalized scores.
#' Data can be either extracted directly from score track or computed via the functions:
#' score_hic_mat, shuffle_and_score_hic_mat
#'
#' @param points_score A dataframe containing the points (start1, start2) and their normalized score.
#' @param rotate Binary flag, indicating if the plot should be rotated by 45 degrees.
#' @param point_size Cex size of the points in the plot.
#' @param add_axis Binary flag, indicating if axis should be added to plot.
#' @return gplot containing the map
#' @export
##########################################################################################################
shaman_gplot_map_score <- function(points_score, interval_range=NA, rotate=TRUE, point_size=0.1, add_axis=TRUE) {
  if (!all(c("start1", "start2", "score") %in% colnames(points_score))) {
    stop("points_score data frame must contain the following columns: start1, start2, score")
  }
  if (is.na(interval_range)[1]) {
	interval_range=data.frame(start=min(points_score$start1), end=max(points_score$start1))
  }
  col.scores = shaman_score_pal()
  if (rotate) {
    points_score = points_score[points_score$start2 > points_score$start1,]
    map_gplot <- ggplot2::ggplot(points_score[order(points_score$score),],
      ggplot2::aes(x=(start1+start2)/2, y=(start2-start1)/2, color=factor(floor(score)))) +
      ggplot2::coord_cartesian(xlim=c(interval_range$start, interval_range$end),
	ylim=c(0, (interval_range$end-interval_range$start)/2)) +
      ggplot2::theme (panel.grid.major = ggplot2::element_blank(),
                      panel.grid.minor = ggplot2::element_blank())
  } else {
    map_gplot <- ggplot2::ggplot(points_score[order(points_score$score),],
      ggplot2::aes(x=start1, y=start2, color=factor(floor(score)))) +
      ggplot2::scale_x_continuous(position = "top") + 
      ggplot2::theme(axis.line.y = ggplot2::element_line(size=0.2))
  }
  map_gplot <- map_gplot +
    ggplot2::geom_point(size=point_size) +
    ggplot2::scale_colour_manual(values=col.scores[101 + sort(unique(floor(points_score$score)))]) +
    ggplot2::scale_x_continuous(expand=c(0,0)) +
    ggplot2::scale_y_continuous(expand=c(0,0)) +
    ggplot2::theme_bw() + ggplot2::theme( legend.position="none",
         panel.border = ggplot2::element_blank(),
         axis.title.x = ggplot2::element_blank(),
         axis.title.y = ggplot2::element_blank(),
	 axis.line.x = ggplot2::element_line(size=0.2))
  if (add_axis == FALSE) {
	map_gplot <- map_gplot + ggplot2::theme(
         axis.ticks = ggplot2::element_blank(),
         axis.text.x = ggplot2::element_blank(),
         axis.text.y = ggplot2::element_blank(),
	 axis.line.x = ggplot2::element_blank(),
	 axis.line.y = ggplot2::element_blank(),
         axis.ticks.length = ggplot2::unit(0,"null") )
  }
  return(map_gplot)
}


##########################################################################################################
#' plots misha 1D tracks and annotations
#'
#' \code{shaman_gplot_tracks_and_annotations}
#'
#' Plots 1D tracks in UCSC format. Along with gene annotations and chromosomal axis, this
#' function plots one dimensional misha tracks (e.g. rna-seq, chip-seq data) and marks
#' annotated intervals in the plotting region.
#'
#'
#' @param genome Name of reference genome (e.g. "hg19", "mm10")
#' @param interval_range 1D interval (chrom, start, end) specifying the region to plot
#' @param misha_tracks List of 1D track expressions (virtual tracks also supported) which can be extracted from.
#' @param mt_colors Array of colors, one for each misha_track, which will be used to plot each 1d track.
#' @param my_ylims Y-axis limits for each misha track. If not provided, generated automatically based on the data
#' to be displayed in the region.
#' @param annotations List of gintervals highlighting annotated regions.
#' @param a_colors Array of colors, one for each annotation set.
#' @param add_genes Boolean flag indicating whether to show gene annotations.
#' @param add_ideogram Boolean flag indicating whether to add an chromosomal ideogram.
#' @param add_axis Boolean flag indicating whether to show chromosmal axis.
#' @param gene_stacking Describes the viewing option of the genes track. Can be either "squish", or "dense"
#' @param gene_size Size of gene annotation view
#' @param track_size Size of track view
#' @param annot_size Size of annotation view
#'
#' @export
##########################################################################################################
shaman_plot_tracks_and_annotations <- function(genome, interval_range,
  misha_tracks=list(), mt_colors=getOption("shaman.track_colors"), mt_ylims=NULL,
  annotations=list(), a_colors=getOption("shaman.annotation_colors"),
  add_genes=T, add_ideogram=T, add_axis=T, gene_stacking="squish", gene_size=0.7,
  track_size=0.8, annotation_size=0.7)
{
  tracks = list()
  if (add_axis) {
    tracks = list(Gviz::GenomeAxisTrack())
  }
  if (add_genes) {
	gene_track <- .shaman_get_ucsc(genome, interval_range, stacking=gene_stacking)
	Gviz::displayPars(gene_track) <- list(size=gene_size)
	tracks[[length(tracks)+1]] <- gene_track
  }
  if (length(misha_tracks) > 0) {
    for (t in 1:length(misha_tracks)) {
      data_track = .shaman_get_data(genome, interval_range, misha_tracks[[t]], col=mt_colors[t], ylim=mt_ylims[[t]])
      Gviz::displayPars(data_track) <- list(size=track_size)
      tracks[[length(tracks)+1]] <- data_track
    }
  }
  if (length(annotations) > 0) {
    for (t in 1:length(annotations)) {
      annot_track = .shaman_get_annotation(genome, interval_range, annotations[[t]],
	fill_col=a_colors[t], border_col=a_colors[t])
      Gviz::displayPars(annot_track) <- list(size=annotation_size)
      tracks[[length(tracks)+1]] <- annot_track
    }
  }
  if (add_ideogram) {
	ideo_track <- Gviz::IdeogramTrack(genome=genome, chromosome=as.character(interval_range$chrom))
	tracks[[length(tracks)+1]] <- ideo_track
  }
  Gviz::plotTracks(tracks, from=interval_range$start, to=interval_range$end, panel.only=T, labelPos="below")
}

##########################################################################################################
#' plots hic normalized map with annotations
#'
#' \code{shaman_plot_map_score_with_annotations}
#'
#' Plots observerved hic contact matrix color-coded by normalized scores, aligned with
#' linear data such as gene annotations, rna, chromatin, transcription factor binding, etc').
#'
#'
#' @param genome Name of reference genome (e.g. "hg19", "mm10")
#' @param points_score A dataframe containing the points (start1, start2) and their normalized score.
#' @param interval_range 1D interval (chrom, start, end) specifying the region to plot
#' @param misha_tracks List of 1D track expressions (virtual tracks also supported) which can be extracted from.
#' @param mt_colors Array of colors, one for each misha_track, which will be used to plot each 1d track.
#' @param my_ylims Y-axis limits for each misha track. If not provided, generated automatically based on the data
#' to be displayed in the region.
#' @param annotations List of gintervals highlighting annotated regions.
#' @param a_colors Array of colors, one for each annotation set.
#' @param add_genes Boolean flag indicating whether to show gene annotations.
#' @param add_ideogram Boolean flag indicating whether to add an chromosomal ideogram.
#' @param add_axis Boolean flag indicating whether to show chromosmal axis.
#' @param gene_stacking Describes the viewing option of the genes track. Can be either "squish", or "dense".
#' @param gene_size Size of gene annotation view.
#' @param track_size Size of track view.
#' @param annot_size Size of annotation view.
#' @param fig_fn Name of png file to output to. Empty string will cause the figure to be plotted to the current device.
#' @param fig_width Width in pixels of output png figure.
#' @param fig_height Height in pixels of output png figure.
#'
#' @export
##########################################################################################################

shaman_plot_map_score_with_annotations <- function(genome, points_score, interval_range, point_size=0.1,
  misha_tracks=list(), mt_colors=getOption("shaman.track_colors"), mt_ylims=NULL,
  annotations=list(), a_colors=getOption("shaman.annotation_colors"),
  add_genes=T, add_ideogram=T, add_axis=T, gene_stacking="squish", gene_size=0.7,
  track_size=0.8, annotation_size=0.7, fig_fn="", fig_width=900, fig_height=5/6*900)
{
  if (!all(c("start1", "start2", "score") %in% colnames(points_score))) {
    stop("points_score data frame must contain the following columns: start1, start2, score")
  }
  map_score <- shaman_gplot_map_score(points_score, interval_range, rotate=TRUE, point_size=point_size, add_axis=FALSE)
  if (fig_fn != "") {
	png(fig_fn, width=fig_width, height=fig_height)
  }
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(5, 1)))
  print(map_score, vp=grid::viewport(layout.pos.row=1:3, layout.pos.col=1))
  grid::pushViewport(grid::viewport(layout.pos.row=4:5, layout.pos.col=1))
  shaman_plot_tracks_and_annotations(genome, interval_range,  misha_tracks=misha_tracks,
	mt_colors=mt_colors, mt_ylims=mt_ylims, annotations=annotations, a_colors=a_colors,
  	add_genes=add_genes, add_ideogram=add_ideogram, add_axis=add_axis, gene_stacking=gene_stacking,
	gene_size=gene_size, track_size=track_size, annotation_size=annotation_size)
  if (fig_fn != "") {
	dev.off()
  }
}

##########################################################################################################
#' Returns the color palette used to display normalized scores 
#'
#' \code{shaman_score_pal}
#'
#' @export
##########################################################################################################

shaman_score_pal <- function() {
  .shaman_check_config(c("shaman.score_pal_0_col", "shaman.score_pal_pos_col",
	"shaman.score_pal_neg_col", "shaman.score_pal_pos_breaks", "shaman.score_pal_neg_breaks"))
  col.neg <- rev(.shaman_palette_breaks(100 , getOption("shaman.score_pal_neg_col"),
    getOption("shaman.score_pal_neg_breaks") ))
  col.pos <- .shaman_palette_breaks(100 , getOption("shaman.score_pal_pos_col"),
    getOption("shaman.score_pal_pos_breaks") )
  col.scores = c(col.neg, getOption("shaman.score_pal_0_col"), col.pos)
  return(col.scores)
}

.shaman_palette_breaks = function(n, colors, breaks)
{
     colspec = colorRampPalette(c(colors[1],colors[1]))(breaks[1])
     for(i in 2:(length(colors)) ){
         colspec = c(colspec, colorRampPalette(c(colors[i-1], colors[i]))(abs(breaks[i]-breaks[i-1])))
     }
     colspec = c( colspec,
          colorRampPalette(c(colors[length(colors)],colors[length(colors)]))(n-breaks[length(colors)])
     )
     return(colspec)
}

.shaman_get_ucsc <- function(genome, interv, stacking="dense") {
  genetrack=Gviz::UcscTrack(track="RefSeq Genes", table="refGene",
	trackType="GeneRegionTrack",chromosome=as.character(interv$chrom), genome=genome,
        rstart="exonStarts", rends="exonEnds", gene="name", symbol="name2",
        transcript="name", strand="strand", name="RefSeq Genes",
        feature="name2", stacking=stacking,
        showId=T, from=interv$start, to=interv$end)
  return(genetrack)

}

.shaman_get_data <- function(genome, interv, track, window=-1, windowSize=500, col="grey", ylim=NULL) {
  d = gextract(track, interv, colnames=c("data"))
  d = d[!is.na(d$data),]
  dtrack = Gviz::DataTrack(start=d$start, end=d$end, data = d$data, chromosome=as.character(interv$chrom),
     genome=genome, type="hist", window=window, windowSize=windowSize, col.histogram=col, fill.histogram=col, ylim=ylim)
  return(dtrack)
}


.shaman_get_annotation <- function(genome, interv, annotations, plot_strand=FALSE,
   fill_col="lightblue", border_col="lightgrey")
{
  if (plot_strand) {
    annot = plyr::ddply(annotations, c("strand"), function(x) {
      y = gintervals.intersect(x, interv)
      y$strand = ifelse(x$strand[1] == 1, "+", "-")
      return(y) })
  } else {
     annot = gintervals.intersect(annotations, interv)
     annot$strand="*"
  }
  if (is.null(annot)) {
      return(null)
  }
  atrack = Gviz::AnnotationTrack(start=annot$start, width=annot$end-annot$start, chromosome=as.character(annot$chrom),
          strand=annot$strand, genome=genome, name="", fill=fill_col, stacking="dense", col=border_col)
  return(atrack)
}

