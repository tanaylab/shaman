# The shaman pacakge - sampling HiC contact matrices for a-parametric normalization#

This package provides a normalization scheme, along with basic analysis and statistical visualization of HiC experimental data. The normalization workflow consists of the following steps:

1. Start with observed HiC data.
1. Shuffle observed data to generate expected data.
1. Compute the normalized score of each observed data point.
1. Visualize normalized contact maps.
1. Quantify contact enrichment acoording to features.

### Installation 
#### Requirements 
- _Perl_
- R packages:
    * _devtools_.
    * _misha_.
    * _RANN_.
    * _data.table_.
    * _ggplot2_.
    * _reshape2_.
    * _grid_.
    * _Gviz (from bioconductor)_.

#### Installing misha package:
```
#!r
install.packages("http://www.wisdom.weizmann.ac.il/~nettam/shaman/misha_3.4.3.tar.gz", repos=NULL) # Download and install misha package
```

#### installing Gviz from bioconductor:
```
#!r
source("https://bioconductor.org/biocLite.R")
biocLite("Gviz")
```

#### Installing shaman package:
Download and install *shaman*: 
```
#!r
devtools::install_bitbucket("tanaylab/shaman", ref='default', vignette = TRUE)
library(shaman)
```

#### Using the package
Please refer to the package vignette for usage and workflow.
```
#!r
browseVignettes('shaman') 
```