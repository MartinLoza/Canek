#' RunCanek
#'
#' Runs Canek integration.
#'
#' @param x object with expression counts or list of matrices.
#' @param batches for S4 objects the column containing batch information.
#' @param slot slot used for Seurat objects (default: data).
#' @param assay assay used for Seurat objects.
#' @param features optional vector of features to use for correction.
#' @param selection.method method used for FindVariableFeatures on Seurat objects when features is NULL.
#' @param nfeatures  number of features returned by SelectIntegrationFeatures.
#' @param fvf.nfeatures number of features returned by FindVariableFeatures.
#' @param integration.name name for the integrated assay.
#' @param debug whether to store information about correction vector.
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
RunCanek.Seurat <- function(x, batches = NULL, slot = "data", assay = NULL, features = NULL, selection.method = "vst", nfeatures = 2000, fvf.nfeatures = 2000, integration.name = "Canek", debug = FALSE, ...) {

  #if not assay is selected, we used the default one
  if(is.null(assay)){
    assay <- Seurat::DefaultAssay(x)
  }

  Seurat::DefaultAssay(x) <- assay
  obj <- Seurat::DietSeurat(x, counts = TRUE, data = TRUE, scale.data = FALSE, assays = assay, misc = FALSE)
  Seurat::VariableFeatures(obj) <- NULL
  obj <- Seurat::SplitObject(obj, split.by = batches)

  if (is.null(features)) {
    features <- Seurat::SelectIntegrationFeatures(obj, nfeatures = nfeatures, fvf.nfeatures = fvf.nfeatures, selection.method = selection.method, verbose = FALSE)
  }

  counts <- lapply(obj, function(xx) {
    Seurat::GetAssayData(xx, slot = slot, assay = assay)[features, ]
  })

  counts <- Canek::CorrectBatches(counts, debug = debug, ...)

  if (debug) {
    info <- counts
    info[["Batches Integrated"]] <- NULL
    counts <- counts[["Batches Integrated"]]
  }

  integrated <- Seurat::CreateAssayObject(counts = counts)

  x[[integration.name]] <- integrated
  Seurat::DefaultAssay(x) <- integration.name

  Seurat::VariableFeatures(x, assay = integration.name) <- features

  if (debug) {
    Seurat::Tool(x) <- info
  }

  Seurat::LogSeuratCommand(x)
}

#' @rdname RunCanek
#' @export
RunCanek.SingleCellExperiment <- function(x, batches = NULL, assay = "logcounts", integration.name = "Canek", debug = FALSE, ...) {
  batches <- split(colnames(x), x[[batches]])
  obj <- lapply(batches, function(batch) {
     x[, batch]
  })

  counts <- lapply(obj, SummarizedExperiment::assay, i = assay)
  counts <- Canek::CorrectBatches(counts, debug = debug, ...)

  if (debug) {
    info <- counts
    info[["Batches Integrated"]] <- NULL
    counts <- counts[["Batches Integrated"]]
  }

  SummarizedExperiment::assays(x, withDimnames = FALSE)[[integration.name]] <- counts
  x
}

#' @rdname RunCanek
#' @export
RunCanek.list <- function(x, ...) {
  objtype <- unique(sapply(lapply(x, class), "[", x = 1))
  if (length(objtype) != 1) stop("Required list of identical object types.")
  switch(objtype,
    "matrix" = CorrectBatches(x, ...),
    stop("When input is a list, it should be a list of matrix objects.")
  )
}
