

#' Title EkfBE
#'
#' Estimation of batch-effect by Extended Kalman Filter
#'
#' @param refBatch Reference batch single-cell data.
#' @param queBatch Query's batch single-cell data.
#' @param pairs A matrix containing MNNs pairs. First column corresponds to query-batch cell indexes.
#' @param sampling Whether or not sampling of MNNs pairs is used on the estimation process.
#' @param numSamples Number of MNNs pairs samples used on the estimation process.
#' @param verbose Print output.
#'
#' @return A list containing the estimated correction vector and MNNs pair samples used on the estimation process.
#' The length of the correction vector is equal to the number of genes.
#'
#' @details refBatch is used as the reference batch and queBatch is used as the query batch.
#' Input batches need to have the same number of genes.
#'

EkfBE <- function(refBatch, queBatch, pairs, sampling=NULL, numSamples= NULL, verbose = FALSE){

  #INIT
  Epochs <- 1
  Pairs_Dim <- dim(pairs)
  Num_Pairs <- Pairs_Dim[1]
  B2_Dim <- dim(queBatch)
  B2_NumCells <- B2_Dim[2]
  Num_genes <- B2_Dim[1]
  Progress <- Num_genes/10
  Progress <- floor(Progress)
  Min_Samples <- 300
  Gain=0.5

  if(Num_Pairs < 20){
    sampling <- FALSE
    warning('\nWarning: Low number of pairs', call. = TRUE)
  }

  #Set number of samples
  if( is.null(numSamples) ){
    if(sampling == TRUE){
      if(Num_Pairs > Min_Samples){
        numSamples <- round(Num_Pairs*0.2)
      }else{
        numSamples <- Min_Samples
      }
    }else{
      numSamples <- Num_Pairs
    }
  }

  if(Num_Pairs < numSamples){
    Epochs <- ceiling(numSamples/Num_Pairs)
    sampling <- FALSE
    warning('\nWarning: Number of pairs is lower than number of samples', call. = TRUE)

    if(verbose)
      cat(paste( "\n\tNumber of epochs: ", Epochs ))
  }

  if (sampling) {
    Samples <- sample(c(1:Num_Pairs),size = numSamples, replace = FALSE)
    Samples <- matrix(c(pairs[Samples,1], pairs[Samples,2] ), nrow = numSamples)
    colnames(Samples) <- colnames(pairs)
  } else {
    numSamples <- Num_Pairs
    Samples <- pairs
  }

  if(verbose)
    cat(paste( '\n\tNumber of pairs used for the estimation:', numSamples ))

  ## INIT VARIABLES
  numSamples_Epoch <- numSamples*Epochs
  Q <- 25
  R <- 100
  # Q <- 2
  # R <- 5
  Correction_Vector <- rep(0,Num_genes)


  S <- rep(0, Num_genes)
  K <- rep(0, Num_genes)
  w <- matrix(0, nrow = Num_genes, ncol = numSamples_Epoch)
  P <- matrix(0, nrow = Num_genes, ncol = numSamples_Epoch)
  x <- matrix(0, nrow = Num_genes, ncol = numSamples_Epoch)
  xs <- matrix(0, nrow = Num_genes, ncol = numSamples_Epoch)


  # P[, 1] <- 10
  P[, 1] <- 100
  Index_Epoch <- 2

  #Estimation using Extendend Kalman Filter over the samples
  for (ep in 1:Epochs){

    if(ep == 1){
      Index_Vector <- c(2:numSamples)
    }else{
      Index_Vector <- c(1:numSamples)
    }

     for (k in Index_Vector) {

      #Reference
      x[, Index_Epoch] = refBatch[, Samples[k,2]]
      # Prediction
      xs[, Index_Epoch] = queBatch[, Samples[k,1]] + w[, Index_Epoch-1]
      #Error
      Error <- x[, Index_Epoch] - xs[, Index_Epoch]
      #Kalman Gain
      S <- P[, Index_Epoch-1] + R
      K <- Gain*P[, Index_Epoch-1]*(1/S)
      #Model Update
      w[, Index_Epoch]<- w[, Index_Epoch-1] + K*Error
      P[, Index_Epoch] <- P[, Index_Epoch-1]-K*P[, Index_Epoch-1]+Q

      Index_Epoch <- Index_Epoch + 1

      }

    }

    Correction_Vector = rowMeans(w)

  Estimation_Data <- list("Correction Vector" = Correction_Vector, "Sampled Pairs" = Samples)

  return(Estimation_Data)
}


#' Correction vector estimation
#'
#' @description Batch effect estimation using the MNNs pairs.
#'
#' @param refBatch Reference batch.
#' @param queBatch Query batch.
#' @param pairs A numerical matrix containing MNNs pairs cell indexes. First column corresponds to query batch cells.
#'
#' @return A list containing the estimated correction vector and the estimation data.
#' The length of the correction vector is equal to the number of genes.
#'
#' @details The input batches must have the same number of genes. The model used on the estimation has the form of g_ref = g_que + be, where
#' the batch effect is represented as a value added to the reference gene expression. The batch effect is estimated as
#' the median of the gene expression difference among the reference and the query batch, e.g. Median(g_ref - g_que).
#'
MedianBE <- function(refBatch, queBatch, pairs ){

  # Model -> g_ref = g_que + be
  # be -> g_ref - g_que

  #Get pairs index
  pQue <- pairs[,1]
  pRef <- pairs[,2]

  #Subsetting
  pRef <- refBatch[,pRef]
  pQue <- queBatch[,pQue]

  be <- pRef-pQue

  return(list("Correction Vector" = rowMedians(be), "Sampled Pairs" = NULL))
}
