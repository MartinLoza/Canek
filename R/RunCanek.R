#' RunCanek
#'
#' Runs Canek integration.
#'
#' @param x object with expression counts or list of matrices.
#' @param batches for S4 objects the column containing batch information.
#' @param slot slot used for Seurat objects (default: data).
#' @param assay assay used for Seurat objects. If NULL (default), an existing "SCT" assay is
#' preferred over the object's default assay.
#' @param features optional vector of features to use for correction.
#' @param selection.method method used for FindVariableFeatures on Seurat objects when features is NULL.
#' @param nfeatures  number of features returned by SelectIntegrationFeatures.
#' @param fvf.nfeatures number of features returned by FindVariableFeatures.
#' @param integration.name name for the integrated assay.
#' @param debug whether to store information about correction vector.
#' @param correctEmbeddings whether to perform the correction on PCA embeddings instead of gene expression (Seurat objects only).
#' @param pcaDim number of PCA dimensions to use when correctEmbeddings is TRUE. If NULL (default),
#' it is inferred from the object's existing "pca" reduction; if none is found, used 30 as default with a warning.
#' When correctEmbeddings is TRUE and a "pca" reduction already exists on the object, it is reused
#' directly instead of being recomputed. If the assay in use is SCTransform-normalized and no "pca"
#' reduction exists, an error is raised instead of computing one internally, since SCTransform fit
#' separately per batch produces residuals that are not on a comparable scale across batches.
#' @param maxLoop number of times to iterate the correction (correctEmbeddings = TRUE
#' only), using each iteration's corrected result as the input to the next. Defaults to 5.
#' @param loopTol average change in median correction magnitude used to stop iterating early. Defaults to 1e-3. Ignored if maxLoop = 1.
#' @param ... additional arguments passed down to methods, e.g. \code{ncores} to parallelize
#' MNN pair finding (see \code{\link{CorrectBatches}}; defaults to 1, sequential).
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
RunCanek.Seurat <- function(x, batches = NULL, slot = "data", assay = NULL, features = NULL, selection.method = "vst", nfeatures = 2000, fvf.nfeatures = 2000, integration.name = "Canek", debug = FALSE, correctEmbeddings = TRUE, pcaDim = NULL, maxLoop = 5, loopTol = 1e-3, ...) {

  #if not assay is selected, prefer an existing SCT assay, otherwise use the current default
  if(is.null(assay)){
    assay <- if("SCT" %in% Seurat::Assays(x)) "SCT" else Seurat::DefaultAssay(x)
  }

  isSCT <- methods::is(x[[assay]], "SCTAssay")

  #SCTransform fit per batch produces residuals on different scales, so a PCA computed internally
  #from those raw values wouldn't be comparable across batches. Require a joint PCA the user has
  #already computed on properly reconciled data instead of computing one ourselves.
  if(isSCT && correctEmbeddings && !("pca" %in% Seurat::Reductions(x))){
    stop(
      "No existing 'pca' reduction found, but assay '", assay, "' is SCTransform-normalized. ",
      "Canek can't safely compute its own PCA from SCTransform residuals fit per batch,
      since they are on different scales (see SCTransform documentation).\n",
      "To implement Canek, first reconcile the normalized data and compute a joint PCA first, e.g.:\n",
      "  features <- SelectIntegrationFeatures(object.list, nfeatures = 3000)\n",
      "  object.list <- PrepSCTIntegration(object.list, anchor.features = features)\n",
      "  x <- merge(object.list[[1]], object.list[[2]], ...)\n",
      "  x <- RunPCA(x)\n",
      "then call RunCanek() again.",
      call. = FALSE
    )
  }

  #if correcting on embeddings, infer pcaDim from an existing PCA reduction unless the user passed one directly
  if(correctEmbeddings && is.null(pcaDim)){
    if("pca" %in% Seurat::Reductions(x)){
      pcaDim <- ncol(Seurat::Embeddings(x, reduction = "pca"))
    } else {
      pcaDim <- 30
      warning("No existing 'pca' reduction found on the object; defaulting pcaDim to 30 for correctEmbeddings. Pass pcaDim explicitly to override.", call. = FALSE)
    }
  }

  Seurat::DefaultAssay(x) <- assay

  if(correctEmbeddings && "pca" %in% Seurat::Reductions(x)){
    #reuse the object's existing joint PCA instead of recomputing one internally 
    emb <- Seurat::Embeddings(x, "pca")[, seq_len(pcaDim), drop = FALSE]
    batchVec <- x[[batches, drop = TRUE]]
    splitCells <- split(rownames(emb), batchVec)
    counts <- lapply(splitCells, function(cn) t(emb[cn, , drop = FALSE]))

    counts <- Canek::CorrectBatches(counts, debug = debug, correctEmbeddings = TRUE,
                                     precomputedEmbeddings = TRUE, pcaDim = pcaDim,
                                     maxLoop = maxLoop, loopTol = loopTol, ...)
  } else {
    obj <- Seurat::DietSeurat(x, counts = TRUE, data = TRUE, scale.data = FALSE, assays = assay, misc = FALSE)
    Seurat::VariableFeatures(obj) <- NULL
    obj <- Seurat::SplitObject(obj, split.by = batches)

    if (is.null(features)) {
      features <- Seurat::SelectIntegrationFeatures(obj, nfeatures = nfeatures, fvf.nfeatures = fvf.nfeatures, selection.method = selection.method, verbose = FALSE)
    }

    counts <- lapply(obj, function(xx) {
      if (packageVersion("Seurat") >= "5.0.0")
        Seurat::GetAssayData(xx, layer = slot, assay = assay)[features, ]
      else
        Seurat::GetAssayData(xx, slot = slot, assay = assay)[features, ]
    })

    if(correctEmbeddings){
      counts <- Canek::CorrectBatches(counts, debug = debug, correctEmbeddings = TRUE, pcaDim = pcaDim, maxLoop = maxLoop, loopTol = loopTol, ...)
    } else {
      counts <- Canek::CorrectBatches(counts, debug = debug, correctEmbeddings = FALSE, ...)
    }
  }

  if (debug) {
    info <- counts
    info[["Batches Integrated"]] <- NULL
    counts <- counts[["Batches Integrated"]]
  }

  ### TEST TEST TEST
  #if correct embeddings, we create an embedding objects
  if(correctEmbeddings == TRUE){
    integrated <-  Seurat::CreateDimReducObject(embeddings = t(counts), assay = assay, key = "Canek_" )
    x[[tolower(integration.name)]] <- integrated
    Seurat::DefaultAssay(x) <- assay
    if(!is.null(features))
      Seurat::VariableFeatures(x, assay = assay) <- features

  }else{
    integrated <- Seurat::CreateAssayObject(counts = counts)
    x[[integration.name]] <- integrated
    Seurat::DefaultAssay(x) <- integration.name
    Seurat::VariableFeatures(x, assay = integration.name) <- features
  }

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
