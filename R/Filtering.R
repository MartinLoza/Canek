#' PairsFiltering
#'
#' Function to filter MNNs pairs
#'
#' @param refBatch Reference batch single-cell data.
#' @param queBatch Query's batch single-cell data.
#' @param pairs A matrix containing MNNs pairs. First column corresponds to query-batch cell indexes.
#' @param verbose Print output.
#'
#' @return A matrix containing the filtered pairs. First column corresponds to query-batch cell indexes.
#'
#' @details Filter MNN pairs by quantiles.
#'
PairsFiltering <- function(refBatch,queBatch, pairs, verbose = FALSE){

    if(verbose)
      cat("\n\nFitering pairs by quantiles")

    v1 <- t(refBatch[,pairs[,2]])
    v2 <- t(queBatch[,pairs[,1]])

    d <- as.matrix(dist(rbind(v1,v2)))
    d <- d[matrix(seq(from = 1, to = (nrow(v1)*2)), ncol = 2)]

    q <- quantile(d)

    IQR <- q["75%"] - q["25%"]

    outliers <- c(which(d < (q["25%"] - (1.5*IQR))), which(d > (q["75%"] + (1.5*IQR))))

    return(pairs[-outliers,])
}

#' FilterRowZeros
#'
#' @param matrix Numeric matrix to filter.
#' @param rowZeros If known, zeros row indexes.
#'
#' @return Matrix with no zeros rows.
FilterRowZeros <- function(matrix = NULL, rowZeros = NULL){

  ## Init
  rNames <- rownames(matrix)

  if(length(rNames) == 0)
    rNames <- seq(1:nrow(matrix))

  if(is.null(rowZeros))
    rowZeros <- WhichRowZeros(matrix = matrix, idxReturn = TRUE)$idx

  matrix <- matrix[-zeroIdx,, drop = FALSE]

  return(matrix)
}

#' WhichRowZeros
#'
#' @param matrix Numeric matrix to filter.
#' @param idxReturn Return the zeros row indexes.
#'
#' @return A list with the zeros rows' names and their indexes.
WhichRowZeros <- function(matrix = NULL, idxReturn = FALSE){

  ## Init
  rNames <- rownames(matrix)
  results <- list()

  if(length(rNames) == 0)
    rNames <- seq(1:nrow(matrix))

  ## Get row zeros names
  idx <- which(matrixStats::rowSums2(matrix) == 0)
  rowZeros <- rNames[idx]

  results[["names"]] <- rowZeros

  if(idxReturn)
    results[["idx"]] <- idx

  return(results)
}














