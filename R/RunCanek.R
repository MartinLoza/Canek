#' RunCanek
#'
#' Runs Canek integration.
#'
#' @param x object with expression counts or list of matrices.
#' @param batches for S4 objects the column containing batch information.
#' @param slot slot used for Seurat objects (default: data).
#' @param assay assay used for Seurat objects (default: RNA).
#' @param features optional vector of features to use for correction.
#' @param selection.method method used for FindVariableFeatures on Seurat objects when features is NULL.
#' @param fvf.nfeatures function used to collapse variable features from different batches. Default is intersect.
#' @param ... additional arguments passed down to methods.
#'
#' @return An object of the appropriate type.
#' @export
#'
#' @rdname RunCanek
RunCanek <- function(x, ...) {
  UseMethod("RunCanek")
}

#' @rdname RunCanek
#' @export
RunCanek.Seurat <- function(x, batches = NULL, slot = "data", assay = "RNA", features = NULL, selection.method = "vst", fvf.nfeatures = 2000, ...) {

  # choose features.
  if (is.null(features)) {
    xl <- Seurat::SplitObject(Seurat::DietSeurat(x), split.by = batches)
    features <- Seurat::SelectIntegrationFeatures(xl, fvf.nfeatures = fvf.nfeatures, selection.method = selection.method, verbose = FALSE)
  }

  # make sure all features make sense.
  features <- features[features %in% rownames(x)]

  y <- x[features, ]

  counts <- Seurat::GetAssayData(y, slot = slot, assay = assay)

  batches <- split(colnames(y), y[[batches]])
  batches <- lapply(batches, function(batch) {
    counts[, batch]
  })

  counts <- RunCanek(batches, ...)
  integrated <- Seurat::CreateAssayObject(data = counts[["Batches Integrated"]])
  x[["Canek"]] <- integrated
  Seurat::DefaultAssay(x) <- "Canek"

  Seurat::VariableFeatures(x, assay = assay) <- features
  x
}

#' @rdname RunCanek
#' @export
RunCanek.SingleCellExperiment <- function(x, batches = NULL, assay = "counts", ...) {
  counts <- SummarizedExperiment::assay(x, assay)

  batches <- split(colnames(x), x[[batches]])
  batches <- lapply(batches, function(batch) {
    counts[, batch]
  })

  counts <- RunCanek(batches, ...)
  SummarizedExperiment::assay(x, "Canek") <- counts[["Batches Integrated"]]

  x
}

#' @rdname RunCanek
#' @export
RunCanek.list <- function(x, ...) {
  Correct_Batches(x, ...)
}
