
# The shaman package - sampling HiC contact matrices for a-parametric normalization

Shaman package impelements functions for resampling Hi-C matrices in order to generate expected contact distributions given constraints on marginal coverage and contact-distance probability distributions. The package also provides support for visualizing normalized matrices and statistical analysis of contact distributions around selected landmarks. It is an R package embedding algorithms implemened in C++, which is built over support from the Tanay's lab *Misha* genomic database system that is provided with the package.

The normalization workflow consists of the following steps:

1.  [Import](https://tanaylab.github.io/shaman/articles/import.html) observed HiC data to misha db. In order to utilize full functionality of the shaman package, you need to construct a Misha database for your genome of interest, and import your Hi-C data to Misha binary Quad-tree format.
2.  Shuffle observed data to generate expected data. Reshuffling of an entire dataset will require 7 hours per 1 billion reads on a machine with one core per chromosome. Note that shaman can also run distributed using Sun Grid Engine.
3.  Compute the normalized score of each observed data point. This computation can be conducted inline for a specific genomic region, or pre-processed on the entire dataset using distributed computation. Score computation on 1 billion reads on a distributed system may take 4-10 hours, depending on the number of cores available.
4.  Visualize normalized contact maps. Precomputed normalized scores allow for flexible analysis of contact distributions.
5.  Quantify contact enrichment acoording to genomic features.

In the example misha database provided in this package we have created a low-footprint matrix to examplify the shaman workflow. We included 4.6 million contacts from ELA K562 dataset covering the hoxd locus (chr2:175e06-178e06) and convergent CTCF regions. Processing the complete matrix from this study requires downloading the full contact list and regenerating the reshuffled matrix.

### Code

Source code can be found at: <https://github.com/tanaylab/shaman>

### Requirements

-   *Perl*
-   *remotes* R package (optional, for automatic installation of *bioconductor* dependencies)
-   *System:* multi-core unix / linux based system or SGE (sun grid engine) cluster for distributed computing are required for large HiC datasets.

### Installation

The quickest way to install *shaman* is to use the following command:

``` r
remotes::install_github('tanaylab/shaman')
```

If remotes fails to install all the bioconductor requirments please install *Gviz* and *GenomeInfoDb* manually from bioconductor:

``` r
source("https://bioconductor.org/biocLite.R")
biocLite("Gviz")
biocLite("GenomeInfoDb")
```

#### When all else fails:

In order to install from source, please take the following steps:

``` r
install.packages("http://www.wisdom.weizmann.ac.il/~nettam/shaman/misha_3.5.6.tar.gz", , repos=NULL) # Download and install misha package) 
source("https://bioconductor.org/biocLite.R") #installing Gviz
biocLite("Gviz")
biocLite("GenomeInfoDb")
install.packages("http://www.wisdom.weizmann.ac.il/~nettam/shaman/shaman_2.0.tar.gz", , repos=NULL) # Download and install shaman package)
```

### Using the package

Please refer to <https://tanaylab.github.io/shaman/articles/shaman-package.html> for usage and workflow.
