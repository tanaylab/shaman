The shaman pacakge - sampling HiC contact matrices for a-parametric normalization
=================================================================================

This package provides a normalization scheme, along with basic analysis and statistical visualization of HiC experimental data. The normalization workflow consists of the following steps:

1.  Start with observed HiC data.
2.  Shuffle observed data to generate expected data.
3.  Compute the normalized score of each observed data point.
4.  Visualize normalized contact maps.
5.  Quantify contact enrichment acoording to features.

### Code

Source code can be found at: <https://bitbucket.org/tanaylab/shaman>

### Installation

#### Requirements

-   *Perl*
-   R packages:
    -   *devtools*.
    -   *misha*.
    -   *RANN*.
    -   *data.table*.
    -   *ggplot2*.
    -   *reshape2*.
    -   *grid*.
    -   *Gviz (from bioconductor)*.

#### Installing misha package:

``` r
install.packages("http://www.wisdom.weizmann.ac.il/~nettam/shaman/misha_3.4.3.tar.gz", repos=NULL) # Download and install misha package
```

#### Installing Gviz from bioconductor:

``` r
source("https://bioconductor.org/biocLite.R")
biocLite("Gviz")
```

#### Installing shaman package:

Download and install *shaman*:

``` r
devtools::install_bitbucket("tanaylab/shaman", ref='default', vignette = TRUE)
library(shaman)
```

#### Using the package

Please refer to the package vignette for usage and workflow.

``` r
browseVignettes('shaman') 
```

Alternatively, see the usage section in the website: <https://tanaylab.bitbucket.io/shaman/articles/shaman-package.html>
