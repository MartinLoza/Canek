# Canek 0.3.0

## Latests updates

- `RunCanek()` on Seurat objects now defaults to `correctEmbeddings = TRUE`:
  correction happens in PCA-embedding space instead of on gene expression
  directly, and the result is stored as a new `"canek"` dimensionality
  reduction rather than a new `"Canek"` assay. Code that reads
  `assay = "Canek"` after `RunCanek()` needs updating — use
  `Embeddings(x, "canek")` instead, or pass `correctEmbeddings = FALSE` to
  keep the original assay-based behavior.

## Changes

- `pcaDim` is now inferred automatically from an existing `"pca"` reduction
  on the object when not specified, instead of always defaulting to 50;
  falls back to 30 with a warning if no `"pca"` reduction exists.

- When `correctEmbeddings = TRUE` and a `"pca"` reduction already exists,
  `RunCanek()` reuses it directly instead of recomputing one internally
  (new `precomputedEmbeddings` parameter on `CorrectBatches()`/
  `CorrectBatch()`).

- Added SCTransform support: `RunCanek()` auto-detects an existing `"SCT"`
  assay, and for `correctEmbeddings = TRUE` requires an existing `"pca"`
  reduction computed on properly reconciled data (see SCTransform documentation). Errors
  with guidance towards the reconciliation steps
  (`SelectIntegrationFeatures()` + `PrepSCTIntegration()`) are thrown.

- Added a iterative-correction loop for `correctEmbeddings = TRUE`: new
  `maxLoop`/`loopTol` parameters (`RunCanek()` on Seurat objects defaults to
  `maxLoop = 5`) iterate the correction, feeding each
  pass's result into the next, and stop early once further passes stop
  improving the correction.

- Fixed a bug where clustering and the fuzzy logic step would error
  (`subscript out of bounds`) whenever `correctEmbeddings = TRUE` was used
  with `pcaDim` below 10.

- Rewrote the Seurat vignette and added a new SCTransform vignette
  demonstrating the current default workflow for each normalization method.
  Pre-0.3.0 vignettes are archived under
  [Previous versions](https://martinloza.github.io/Canek/articles/legacy.html).

# Canek 0.2.5

## Changes

- Update unit tests for compatibility with changes in igraph ([#24](https://github.com/MartinLoza/Canek/issues/24)).

# Canek 0.2.4

## Changes

- Update Canek for compatibility with Seurat v5 ([#20](https://github.com/MartinLoza/Canek/issues/20))

# Canek 0.2.3

## Changes

- Update Canek for compatibility with Seurat v5 ([#20](https://github.com/MartinLoza/Canek/issues/20))

# Canek 0.2.2

## Changes

- Fix bug when using SCTransform.
([#14](https://github.com/MartinLoza/Canek/issues/14))

- Add BugReports and URL.([#13](https://github.com/MartinLoza/Canek/pull/13))

- Relaxing tests to no consider the names attribute.([#12](https://github.com/MartinLoza/Canek/issues/12))

- Refactor `RunCanek()` Seurat and SingleCellExperiment interfaces to preserve original objects.

- Added integration.name argument to `RunCanek()` to customize the name of the integrated assay.

# Canek 0.2.1

## Changes

- Fix error in `test-Clustering` unit-test for latest `igraph` version. ([#7](https://github.com/MartinLoza/Canek/issues/7))

# Canek 0.2.0

- Initial release to CRAN
