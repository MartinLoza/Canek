

#' Title EKF_BE
#'
#' Estimation of batch-effect by Extended Kalman Filter
#'
#' @param B1 Reference batch single-cell data.
#' @param B2 Query's batch single-cell data.
#' @param Pairs A matrix containing MNNs pairs. First column corresponds to query-batch cell indexes.
#' @param Sampling Whether or not sampling of MNNs pairs is used on the estimation process.
#' @param Number_Samples Number of MNNs pairs samples used on the estimation process.
#' @param Gain Gain used when updating the estimated values.
#' @param Verbose Print output.
#'
#' @return A list containing the estimated correction vector and MNNs pair samples used on the estimation process.
#' The length of the correction vector is equal to the number of genes.
#'
#' @details B1 is used as the reference batch and B2 is used as the query batch.
#' Input batches need to have the same number of genes.
#'

EKF_BE <- function(B1,
                   B2,
                   Pairs,
                   Sampling=NULL,
                   Number_Samples= NULL,
                   Gain=0.5,
                   Verbose = FALSE
                   ){

  #INIT
  Epochs <- 1
  Pairs_Dim <- dim(Pairs)
  Num_Pairs <- Pairs_Dim[1]
  B2_Dim <- dim(B2)
  B2_NumCells <- B2_Dim[2]
  Num_genes <- B2_Dim[1]
  Progress <- Num_genes/10
  Progress <- floor(Progress)
  Min_Samples <- 300

  if(Num_Pairs < 20){
    Sampling <- FALSE
    warning('\nWarning: Low number of pairs', call. = TRUE)
  }

  #Set number of samples
  if( is.null(Number_Samples) ){
    if(Sampling == TRUE){
      if(Num_Pairs > Min_Samples){
        Number_Samples <- round(Num_Pairs*0.2)
      }else{
        Number_Samples <- Min_Samples
      }
    }else{
      Number_Samples <- Num_Pairs
    }
  }

  if(Num_Pairs < Number_Samples){
    Epochs <- ceiling(Number_Samples/Num_Pairs)
    Sampling <- FALSE
    warning('\nWarning: Number of pairs is lower than number of samples', call. = TRUE)

    if(Verbose)
      cat(paste( "\n\tNumber of epochs: ", Epochs ))
  }

  if (Sampling) {
    Samples <- sample(c(1:Num_Pairs),size = Number_Samples, replace = FALSE)
    Samples <- matrix(c(Pairs[Samples,1], Pairs[Samples,2] ), nrow = Number_Samples)
    colnames(Samples) <- colnames(Pairs)
  } else {
    Number_Samples <- Num_Pairs
    Samples <- Pairs
  }

  if(Verbose)
    cat(paste( '\n\tNumber of pairs used for the estimation:', Number_Samples ))

  ## INIT VARIABLES
  Number_Samples_Epoch <- Number_Samples*Epochs
  Q <- 25
  R <- 100
  # Q <- 2
  # R <- 5
  Correction_Vector <- rep(0,Num_genes)


  S <- rep(0, Num_genes)
  K <- rep(0, Num_genes)
  w <- matrix(0, nrow = Num_genes, ncol = Number_Samples_Epoch)
  P <- matrix(0, nrow = Num_genes, ncol = Number_Samples_Epoch)
  x <- matrix(0, nrow = Num_genes, ncol = Number_Samples_Epoch)
  xs <- matrix(0, nrow = Num_genes, ncol = Number_Samples_Epoch)


  # P[, 1] <- 10
  P[, 1] <- 100
  Index_Epoch <- 2

  #Estimation using Extendend Kalman Filter over the samples
  for (ep in 1:Epochs){

    if(ep == 1){
      Index_Vector <- c(2:Number_Samples)
    }else{
      Index_Vector <- c(1:Number_Samples)
    }

     for (k in Index_Vector) {

      #Reference
      x[, Index_Epoch] = B1[, Samples[k,2]]
      # Prediction
      xs[, Index_Epoch] = B2[, Samples[k,1]] + w[, Index_Epoch-1]
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

  Estimation_Data <- list("Sampled Pairs" = Samples ,"Correction Vector" = Correction_Vector)

  return(Estimation_Data)
}


Sub_BE <- function(B1,
                   B2,
                   Pairs,
                   Verbose = FALSE
                   ){

  # Model -> g_ref = g_que + be
  # be -> g_ref - g_que

  #Get pairs index
  pQue <- Pairs[,1]
  pRef <- Pairs[,2]

  #Subsetting
  pRef <- B1[,pRef]
  pQue <- B2[,pQue]

  be <- pRef-pQue

  return(list("Correction Vector" = rowMeans(be), "Sampled Pairs" = NULL))
}
