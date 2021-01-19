
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Canek

<!-- badges: start -->

[![R-CMD-check](https://github.com/MartinLoza/Canek/workflows/R-CMD-check/badge.svg)](https://github.com/MartinLoza/Canek/actions)
<!-- badges: end -->

The goal of Canek is to …

## Installation

<!--
You can install the released version of Canek from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("Canek")
```
-->

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("MartinLoza/Canek")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(Canek)
#> Registered S3 method overwritten by 'spdep':
#>   method   from
#>   plot.mst ape

res <- RunCanek(SimBatches$batches)
res[1:5, 1:5]
#>        Cell701  Cell702  Cell703  Cell704  Cell705
#> Gene1 6.664089 5.711871 6.277570 6.117087 6.720634
#> Gene2 6.962410 6.934708 7.857783 6.586155 7.419143
#> Gene3 6.063297 5.654124 5.545338 6.164500 5.734249
#> Gene4 6.641050 5.922194 5.971686 6.360457 6.731420
#> Gene5 4.503925 3.358179 1.160355 2.163020 3.391893
```

## Citing Canek
If you use Canek, please cit our paper..... 


