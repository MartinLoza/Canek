

#' Title PairsFiltering
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














