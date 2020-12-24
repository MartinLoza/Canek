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
  x <- Seurat::SplitObject(x, split.by = batches)
  RunCanek(x, slot = slot, assay = assay, features = features, selection.method = selection.method, fvf.nfeatures = fvf.nfeatures, ...)
}

#' @rdname RunCanek
#' @export
RunCanek.SingleCellExperiment <- function(x, batches = NULL, assay = "counts", ...) {
  batches <- split(colnames(x), x[[batches]])
  x <- lapply(batches, function(batch) {
    x[, batch]
  })

  RunCanek(x, assay = assay, ...)
}

#' @rdname RunCanek
#' @export
RunCanek.list <- function(x, ...) {
  objtype <- unique(sapply(lapply(x, class), "[", x = 1))
  if (length(objtype) != 1) stop("Required list of identical object types.")
  switch(objtype,
    "Seurat" = RunCanekSeurat(x, ...),
    "SingleCellExperiment" = RunCanekSingleCellExperiment(x, ...),
    "matrix" = Correct_Batches(x, ...)
  )

}

RunCanekSeurat <- function(x, slot = "data", assay = "RNA", features = NULL, selection.method = "vst", nfeatures = 2000, fvf.nfeatures = 2000, ...) {

  if (is.null(features)) {
    features <- Seurat::SelectIntegrationFeatures(x, nfeatures = nfeatures, fvf.nfeatures = fvf.nfeatures, selection.method = selection.method, verbose = FALSE)
  }

  counts <- lapply(x, function(xx) {
    Seurat::GetAssayData(xx, slot = slot, assay = assay)[features, ]
  })

  counts <- Canek::Correct_Batches(counts, ...)
  integrated <- Seurat::CreateAssayObject(counts = counts)
  x <- Reduce(merge, x)

  x[["Canek"]] <- integrated
  Seurat::DefaultAssay(x) <- "Canek"

  Seurat::VariableFeatures(x, assay = "Canek") <- features
  Seurat::LogSeuratCommand(x)
}

RunCanekSingleCellExperiment <- function(x, assay = NULL, ...) {
  counts <- lapply(x, SummarizedExperiment::assay, i = assay)
  counts <- Canek::Correct_Batches(counts, ...)


  x <- Reduce(SummarizedExperiment::cbind, x)
  SummarizedExperiment::assays(x, withDimnames = FALSE)[["Canek"]] <- counts

  x
}
