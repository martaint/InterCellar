
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![install with
bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/bioconductor-intercellar/README.html)

<!-- badges: end -->

# InterCellar <img src="inst/app/www/logo_lowres.png" align="right" alt="" width="150" />

an R/Shiny app for interactive analysis and exploration of cell-cell
communication based on single-cell transcriptomics data

## Description

`InterCellar` allows researchers to interactively analyze the results of
cell-cell communication from scRNA-seq data. Starting from pre-computed
ligand-receptor interactions, `InterCellar` provides filtering options,
annotations and multiple visualizations to explore clusters, genes and
functions. Moreover, based on functional annotation from Gene Ontology
and pathway databases, `InterCellar` implements data-driven analyses to
investigate cell-cell communication in one or multiple conditions.

Every step of the analysis can be performed interactively, thus not
requiring any programming skills. Moreover, `InterCellar` runs on your
local machine, avoiding issues related to data privacy.

## Installation

### Bioconductor

`InterCellar` is distributed as a
[Bioconductor](https://www.bioconductor.org/) package and requires R
(version 4.1) and Bioconductor (version 3.14).

To install `InterCellar` package enter:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("InterCellar")
```

### Bioconda

Alternatively, `InterCellar` can be installed through Bioconda. We
recommend installing `InterCellar` in a fresh environment, such as:

``` r
conda create --name=intercellar_env 
conda activate intercellar_env
conda install bioconductor-intercellar
```

Once the installation is done, you can start R simply by

``` r
R
```

## Launching the app

Once `InterCellar` is successfully installed, it can be loaded as
follow:

``` r
library(InterCellar)
```

In order to start the app, please run the following command:

``` r
InterCellar::run_app( reproducible = TRUE )
```

`InterCellar` should be opening in a browser. If this does not happen
automatically, please open a browser and navigate to the address shown
(for example, `Listening on http://127.0.0.1:6134`). The flag
`reproducible = TRUE` ensures that your results will be reproducible
across R sessions.

## Troubleshooting

### Bioconda

For Bioconda users, running `InterCellar` might fail due to this error:

``` r
Error in utils::browseURL(appUrl) : 
  'browser' must be a non-empty character string
```

Try this solution:

``` r
options(browser="firefox")

# and then as usual
InterCellar::run_app( reproducible = TRUE )
```

## User Guide

First time here? Please have a look at `InterCellar` user guide
[here](http://bioconductor.org/packages/devel/bioc/vignettes/InterCellar/inst/doc/user_guide.html).

## Paper reproducibility

Please have a look at
(InterCellar-reproducibility)\[<https://github.com/martaint/InterCellar-reproducibility>\]
if you are interested in data and results showed in the
(manuscript)\[<https://www.researchsquare.com/article/rs-525466/v1>\].

## Help and Suggestions

If you have any question, problem or suggestion, please feel free to
open an [issue](https://github.com/martaint/InterCellar/issues) or
contact Marta Interlandi at <marta.interlandi@uni-muenster.de>

## Citation

Marta Interlandi, Kornelius Kerl, Martin Dugas et al. Interactive
analysis and exploration of cell-cell communication in single-cell
transcriptomics with InterCellar, 21 June 2021, PREPRINT (Version 1)
available at Research Square
\[<https://doi.org/10.21203/rs.3.rs-525466/v1>\]

## Code of Conduct

Please note that the InterCellar project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
