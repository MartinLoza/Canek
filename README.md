
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

*Canek is an R package to correct batch effects from single-cell RNA-seq
biological replicates*.

### Motivation to develop Canek

As single-cell genomics technologies become mainstream, more
laboratories will perform experiments under different conditions with
biological replicates obtained using a common technology. In this
scenario, integration of datasets with minimal impact on cell phenotype
is essential.

### The workflow

Canek leverages information from mutual nearest neighbor to combine
local linear corrections with cell-specific non-linear corrections
within a fuzzy logic framework.

<p align="center">
<img src="man/figures/README-workflow.png" width="100%"/>
</p>
<style type="text/css">
    ol { list-style-type: upper-alpha; }
</style>

<font size="2">

> A. Canek starts with a reference batch and query batch, assuming a
> predominantly linear batch effect.
>
> B. Cell clusters are defined on the query batch and MNN pairs (arrows)
> are used to define batch effect observations.
>
> C. The MNN pairs from each cluster are used to estimate cluster
> specific correction vectors. These vectors can be used to correct the
> batch effect or, (D) a non-linear correction can be applied by
> calculating cell-specific correction vectors using fuzzy logic.

<font size="3">

### Results

Canek was the highest scored method in tests specifically designed to
assess **over-correction**, where Canek corrected batch effects without
distortion to the structures of cells as compared with a gold standard.

For more information about Canek check out our manuscript in **NAR
Genomics and Bioinformatics**:

-   [Canek
    manuscript](https://academic.oup.com/nargab/article/4/1/lqac022/6548822?login=true)
-   [PDF](https://academic.oup.com/nargab/article-pdf/4/1/lqac022/42899055/lqac022.pdf)

## Usage

You can use Canek directly with *normalized-count matrices*, *Seurat*
objects or *SingleCellExperiment* objects. For more details, check out
our GitHub page and vignettes:

-   [Canek website](https://martinloza.github.io/Canek/index.html)
-   [Run Canek on a toy example
    vignette](https://martinloza.github.io/Canek/articles/toy_example.html)
-   [Run Canek on Seurat objects
    vignette](https://martinloza.github.io/Canek/articles/seurat.html)

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

## Citation

If you use Canek in your research please cite our work using:

Loza M, Teraguchi S, Standley DM, Diez D (2022). “Unbiased integration
of single cell transcriptome replicates.” *NAR Genomics and
Bioinformatics*, *1*(4), lqac022. doi: 10.1093/nargab/lqac022 (URL:
<https://doi.org/10.1093/nargab/lqac022>).
