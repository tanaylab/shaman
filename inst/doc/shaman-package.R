## ---- eval=FALSE---------------------------------------------------------
#  install.packages("http://www.wisdom.weizmann.ac.il/~nettam/shaman/misha_3.4.3.tar.gz", repos=NULL) #Installs misha package from file

## ---- eval=FALSE---------------------------------------------------------
#  source("https://bioconductor.org/biocLite.R")
#  biocLite("Gviz")

## ---- eval=FALSE---------------------------------------------------------
#  devtools::install_bitbucket("tanaylab/shaman")
#  library(shaman)

## ---- eval=FALSE, warning=FALSE------------------------------------------
#  shaman_dump_config(config_dir = "my_shaman_proj/conf")

## ---- eval=FALSE, warning=FALSE------------------------------------------
#  shaman_load_config("my_shaman_proj/conf")

## ---- eval=FALSE, warning=FALSE------------------------------------------
#  shuffle_hic_track(track_db="db", obs_track_nm="obs", work_dir="work_dir")

## ---- eval=FALSE, warning=FALSE------------------------------------------
#  score_hic_track(track_db="db", work_dir="work_dir", score_track_nm="score_track", obs_track_nms=c("obs"))

## ---- eval=FALSE, warning=FALSE------------------------------------------
#  point_score = gextract("score_track", region, colnames="score")
#  plot_map_score_with_annotations("hg19", point_score$points, region, misha_tracks=list("K56.k27ac", "rna-seq"), annotations=list("ctcf_pos", "ctcf_neg"), a_colors=c("#4572A7", "#AA4643"))

## ---- eval=FALSE, warning=FALSE------------------------------------------
#  point_score = shuffle_and_score_hic_mat(obs_track_nms="obs", interval=interval, work_dir="work_dir")

## ---- eval=FALSE, warning=FALSE------------------------------------------
#  plot_map_score_with_annotations("hg19", point_score$points, region, misha_tracks=list("K56.k27ac", "rna-seq"), annotations=list("ctcf_pos", "ctcf_neg"), a_colors=c("#4572A7", "#AA4643"))

## ---- eval=FALSE, warning=FALSE------------------------------------------
#  grid = shaman_generate_feature_grid(feature1, feature2, obs_track, exp_track, range=25000, resolution=500)
#  shaman_plot_feature_grid(grid, range=25000, grid_resolution=500, plot_resolution=1000)

