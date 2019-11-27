
##EKF_BE##
#Function to estimate the batch effect from two input batches. B1 is used as the reference batch and B2 is used as the query batch. These two
#input batches need to have the same genes.
# INPUT : B1 -> Batch 1 (Reference)
#         B2 -> Batch 2 (Query)
#         Pairs -> Cell pairs for batch effect estimations
#         Sampling -> A boolean value indicating if the estimator will use sampling from the cell pairs
#         Number_Samples -> Number of samples to use from the cell pairs
#         Gain -> Gain to use on the estimation process
#
# OUTPUT : Correction Vector containing the batch effect estimation. The vector is equal size as the number of genes.
EKF_BE <- function(B1,B2, Pairs, Sampling=NULL, Number_Samples= NULL, Gain=0.1){
  #Gain=0.8

  Epochs <- 1

  #Pairs Dimensions
  Pairs_Dim <- dim(Pairs)
  Num_Pairs <- Pairs_Dim[1]
  #Pairs <- matrix(c(as.vector(Pairs[1:Num_Pairs,1]), as.vector(Pairs[1:Num_Pairs,2])), nrow = Num_Pairs)

  if(Num_Pairs < 10){
    Sampling <- FALSE
    warning('\nWarning: Low number of pairs', call. = TRUE)

  }

  #Query Batch Dimensions
  B2_Dim <- dim(B2)
  B2_NumCells <- B2_Dim[2]
  Num_genes <- B2_Dim[1]

  Progress <- Num_genes/10
  Progress <- floor(Progress)

  #Check number of samples
  if( is.null(Number_Samples) ){
    if(Sampling == TRUE){
      Number_Samples <- round(Num_Pairs*0.2)
    }else{
      Number_Samples <- 1000
    }
  }

  if(Num_Pairs < Number_Samples){
    Epochs <- ceiling(Number_Samples/Num_Pairs)
    Number_Samples = Num_Pairs
    Sampling <- FALSE
    #Epochs <- ceiling(1000/Number_Samples)
    warning('\nWarning: Number of pairs is lower than number of samples', call. = TRUE)
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
  cat(paste( '\n\tNumber of pairs used for the estimation:', Number_Samples ))

  ## INIT VARIABLES
  Number_Samples_Epoch <- Number_Samples*Epochs
  #Init weight
  w <- c(rep(0,Number_Samples_Epoch))
  #Init P
  P <- c(rep(0,Number_Samples_Epoch))
  #Init Q
  Q <- 2
  R <- 5
  #Kalman Gain
  K <- 0
  #Extra
  S <- 0
  #Estimation variables
  x <- c(rep(0,Number_Samples_Epoch))
  xs <- c(rep(0,Number_Samples_Epoch))
  Correction_Vector <- rep(0,Num_genes)


  #Estimations for all genes
  for (i in 1:Num_genes)
  {
    #Init Variables
    P[1] <- 10
    x[1] <- 0
    xs[1] <- 0
    Index_Epoch <- 2

    #Gene to analize
    Gene <- i

    if (Gene == 1){
      cat( paste( "\n\tEstimating Correction Vector. Progress:", 0, "%" ) )
    }else{
      if(mod(Gene,Progress) == 0){
        cat( paste( "\r\tEstimating Correction Vector. Progress:", (Gene/Progress*10) , "%" ) )
      }
    }

    #Estimation using Extendend Kalman Filter over the samples
    for (ep in 1:Epochs){

      if(ep == 1){
        Index_Vector <- c(2:Number_Samples)
      }else{
        Index_Vector <- c(1:Number_Samples)
      }

      for (k in Index_Vector) {

        #Reference
        x[Index_Epoch] = B1[Gene,Samples[k,2]]
        # Prediction
        xs[Index_Epoch] = B2[Gene,Samples[k,1]] + w[Index_Epoch-1]
        #Error
        Error <- x[Index_Epoch] - xs[Index_Epoch]
        #Kalman Gain
        S <- P[Index_Epoch-1] + R
        K <- Gain*P[Index_Epoch-1]*(1/S)
        #Model Update
        w[Index_Epoch]<- w[Index_Epoch-1] + K*Error
        P[Index_Epoch] <- P[Index_Epoch-1]-K*P[Index_Epoch-1]+Q

        Index_Epoch <- Index_Epoch + 1

      }

    }

    Correction_Vector[i] = mean(w)
  }

  Estimation_Data <- list("Sampled Pairs" = Samples ,"Correction Vector" = Correction_Vector)

  return(Estimation_Data)
}
