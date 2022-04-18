
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!--![Canek_logo](images/logo.png){width="100"}<!-- -->
<p align="center">
<img src="man/figures/README-logo.png" width="50%"  />
</p>

# Canek

<!-- badges: start -->

[![R-CMD-check](https://github.com/MartinLoza/Canek/workflows/R-CMD-check/badge.svg)](https://github.com/MartinLoza/Canek/actions)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/Canek)](https://cran.r-project.org/package=Canek)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/Canek)](https://cran.r-project.org/package=Canek)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/Canek)](https://cran.r-project.org/package=Canek)
<!-- badges: end -->

## Citation

If you use Canek in your research please cite our work using:

Loza M, Teraguchi S, Standley DM, Diez D (2022). “Unbiased integration
of single cell transcriptome replicates.” *NAR Genomics and
Bioinformatics*, *1*(4), lqac022. doi: 10.1093/nargab/lqac022 (URL:
<https://doi.org/10.1093/nargab/lqac022>).

-   [Canek
    manuscript](https://academic.oup.com/nargab/article/4/1/lqac022/6548822?login=true)
-   [PDF](https://academic.oup.com/nargab/article-pdf/4/1/lqac022/42899055/lqac022.pdf)

## Installation

You can install the release version of Canek from
[CRAN](https://CRAN.R-project.org) with:

    install.packages("Canek")

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("MartinLoza/Canek")
```

## Example

This is a basic example:

``` r
library(Canek)

res <- RunCanek(SimBatches$batches)
res[1:5, 1:5]
#>          Cell1    Cell2    Cell3    Cell4    Cell5
#> Gene1 6.145644 6.228590 6.960904 6.665023 6.586427
#> Gene2 6.919552 7.048549 7.328128 7.649456 7.005002
#> Gene3 5.468337 6.115350 6.146068 6.085750 6.142014
#> Gene4 6.810614 7.076788 6.872226 6.359445 6.303407
#> Gene5 2.848188 2.160687 4.454532 2.467829 3.165861
```

For more tutorials using Canek check out our github page and vignettes:

-   [Canek website](https://martinloza.github.io/Canek/index.html)
-   [Run Canek on a toy example
    vignette](https://martinloza.github.io/Canek/articles/toy_example.html)
-   [Run Canek on Seurat objects
    vignette](https://martinloza.github.io/Canek/articles/seurat.html)
