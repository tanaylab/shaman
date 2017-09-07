
The shaman pacakge - sampling HiC contact matrices for a-parametric normalization
=================================================================================

This package provides a normalization scheme, along with basic analysis and statistical visualization of HiC experimental data. The normalization workflow consists of the following steps:

1.  [Import](https://tanaylab.bitbucket.io/shaman/articles/import.html) observed HiC data to misha db.
2.  Shuffle observed data to generate expected data.
3.  Compute the normalized score of each observed data point.
4.  Visualize normalized contact maps.
5.  Quantify contact enrichment acoording to genomic features.

### Code

Source code can be found at: <https://bitbucket.org/tanaylab/shaman>

### Requirements

-   *Perl*
-   *devtools* R package (optional, for automatic installation of *bioconductor* dependencies)
-   *System:* multi-core unix / linux based system or SGE (sun grid engine) cluster for distributed computing are required for large HiC datasets.

### Installation

The quickest way to install *shaman* is use the following command:

``` r
devtools::install_bitbucket('tanaylab/shaman@default', 
                            vignette = TRUE, 
                            repos=c(getOption('repos'), 'https://tanaylab.bitbucket.io/repo'))
```

If devtools fails to install all the bioconductor requirments please install *Gviz* and *GenomeInfoDb* manually from bioconductor:

``` r
source("https://bioconductor.org/biocLite.R")
biocLite("Gviz")
biocLite("GenomeInfoDb")
```

#### Installing without the *devtools* package

In order to install without the *devtools* package, please take the following steps:

1.  Install *Gviz* and *GenomeInfoDb* from *bioconductor*:

``` r
source("https://bioconductor.org/biocLite.R")
biocLite("Gviz")
biocLite("GenomeInfoDb")
```

1.  Install *shaman*:

``` r
install.packages('shaman', repos=c(getOption('repos'), 'https://tanaylab.bitbucket.io/repo')) 
```

### Using the package

Please refer to <https://tanaylab.bitbucket.io/shaman/articles/shaman-package.html> for usage and workflow.
