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

#### Requirements

-   *Perl*
-   *devtools* R package
-   *Gviz* package from bioconductor

### Installation

The quickest way to install *shaman* is use the following command:

``` r
devtools::install_bitbucket('tanaylab/shaman@default', 
                            vignette = TRUE, 
                            repos=c(getOption('repos'), 'https://tanaylab.bitbucket.io/repo'))
```

If devtools fails to install all the bioconductor requirments please install *Gviz* manually from bioconductor:

``` r
source("https://bioconductor.org/biocLite.R")
biocLite("Gviz")
```

#### Using the package

Please refer to the package vignette for usage and workflow.

``` r
browseVignettes('shaman') 
```

Alternatively, see the usage section in the website: <https://tanaylab.bitbucket.io/shaman/articles/shaman-package.html>
