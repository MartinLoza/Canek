

#' Title Pairs_Selection
#'
#' Function to filter MNNs pairs
#'
#' @param B1 Reference batch single-cell data.
#' @param B2 Query's batch single-cell data.
#' @param Pairs A matrix containing MNNs pairs. First column corresponds to query-batch cell indexes.
#' @param Verbose Print output.
#'
#' @return A matrix containing the filtered pairs. First column corresponds to query-batch cell indexes.
#'
#' @details Filter pairs to used on batch effect correction by quantiles.
#'
Pairs_Selection <- function(B1,B2, Pairs, Verbose = FALSE){

    if(Verbose)
      cat("\n\nFitering pairs by quantiles")

    v1 <- t(B1[,Pairs[,2]])
    v2 <- t(B2[,Pairs[,1]])

    d <- as.matrix(dist(rbind(v1,v2)))
    d <- d[matrix(seq(from = 1, to = (nrow(v1)*2)), ncol = 2)]

    q <- quantile(d)

    IQR <- q["75%"] - q["25%"]

    outliers <- c(which(d < (q["25%"] - (1.5*IQR))), which(d > (q["75%"] + (1.5*IQR))))

    selPairs <- Pairs[-outliers,]

    return(list("Selected Pairs" = selPairs,
                "Clusters" = NULL,
                "Selected Cluster" = NULL))

}














