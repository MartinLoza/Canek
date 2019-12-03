#' RunCanek
#'
#' Runs Canek integration.
#'
#' @param x object with expression counts or list of matrices.
#' @param batches for S4 objects the column containing batch information.
#' @param slot slot used for Seurat objects (default: data).
#' @param assay assay used for Seurat objects (default: RNA).
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
RunCanek.Seurat <- function(x, batches = NULL, slot = "data", assay = "RNA", ...) {
  counts <- GetAssayData(x, slot = slot, assay = assay)
  features <- VariableFeatures(x)

  batches <- split(colnames(x), x[[batches]])
  batches <- lapply(batches, function(batch) {
    counts[, batch]
  })

  counts <- RunCanek(batches)
  integrated <- CreateAssayObject(data = counts[["Batches Integrated"]])
  x[["Canek"]] <- integrated
  DefaultAssay(x) <- "Canek"

  VariableFeatures(x) <- features
  x
}

#' @rdname RunCanek
#' @export
RunCanek.list <- function(x, ...) {
  Correct_Batches(x)
}
