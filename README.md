
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!--![Canek_logo](images/logo.png){width="100"}<!-- -->
<p align="center">
<img src="man/figures/README-logo.png" width="50%"  />
</p>

# Canek

<!-- badges: start -->

[![R-CMD-check](https://github.com/MartinLoza/Canek/workflows/R-CMD-check/badge.svg)](https://github.com/MartinLoza/Canek/actions)
<!-- badges: end -->

## Citation

Please cite Canek if you use it in your research with:

*Loza,M., Teraguchi,S., Standley,D.M. and Diez,D. (2022) Unbiased
integration of single cell transcriptome replicates. Nar Genom
Bioinform, 4, lqac022.*

Check out the manuscript in NAR Genomics and Bioinformatics:

-   [Canek
    manuscript](https://academic.oup.com/nargab/article/4/1/lqac022/6548822?login=true)
-   [PDF](https://academic.oup.com/nargab/article-pdf/4/1/lqac022/42899055/lqac022.pdf)

## Installation

<!-- You can install the released version of Canek from [CRAN](https://CRAN.R-project.org) with: -->
<!-- ``` r -->
<!-- install.packages("Canek") -->
<!-- ``` -->

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
#> Gene1 6.153816 6.216487 6.955690 6.680703 6.599896
#> Gene2 6.925125 7.046951 7.335873 7.623013 7.017121
#> Gene3 5.451978 6.140471 6.149198 6.046806 6.158075
#> Gene4 6.793623 7.044712 6.857325 6.274883 6.287759
#> Gene5 2.791353 2.271425 4.410185 2.571911 3.155063
```

For more tutorials using Canek check out our github page and vignettes:

-   [Canek website](https://martinloza.github.io/Canek/index.html)
-   [Run Canek on a toy example
    vignette](https://martinloza.github.io/Canek/articles/toy_example.html)
-   [Run Canek on Seurat objects
    vignette](https://martinloza.github.io/Canek/articles/seurat.html)
