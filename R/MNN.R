##FindMnnPairs##
#Find MNN pairs given two matrices containing nearest neighbors
# INPUT :
#
# OUTPUT :
FindMnnPairs <- function(B1_B2_NN = NULL, B2_B1_NN = NULL, B2_NCells = NULL ){
  # INIT
  m_Pairs <- NULL

  if ( is.null(B1_B2_NN) ){
    stop('B1_B2_NN, Nearest neighbors need to be defined')
  }

  if ( is.null(B2_B1_NN) ){
    stop('B2_B1_NN, Nearest neighbors need to be defined')
  }

  if ( is.null(B2_NCells) ){
    stop('B2_NCells, Number of queBatch cells needs to be defined.')
  }

  #CHECK PAIRS
  for (i in 1:B2_NCells) {

    p_B1_B2 <- which(B1_B2_NN[,2] == i)
    if(length(p_B1_B2) == 0){
      next()
    }
    B1_B2_sub <- matrix(B1_B2_NN[p_B1_B2,], ncol = 2)

    p_B2_B1 <- which(B2_B1_NN[,1] == i)
    if(length(p_B2_B1) == 0){
      next()
    }
    B2_B1_sub <- matrix( B2_B1_NN[p_B2_B1,], ncol = 2)

    d <- dim(B1_B2_sub)[1]

    for (j in 1:d) {

      m_pair <- which( B2_B1_sub[,2] == B1_B2_sub[j,1] )

      if( length( m_pair ) != 0){
        m_Pairs <- rbind( m_Pairs, B2_B1_sub[m_pair,]  )
      }
    }
  }

  colnames(m_Pairs) <- c('queBatch-Cells-Index', 'refBatch-Cells-Index' )

  return(list("Pairs" = m_Pairs))
}

##Get_MNN_Pairs##
#Get MNN given two batches
# INPUT :
#
# OUTPUT :
GetMnnPairs <- function(refBatch = NULL, queBatch = NULL, kNN = 25){

  if ( is.null(refBatch) ){
    stop('refBatch, Batch needs to be defined')
  }
  if ( is.null(queBatch) ){
    stop('queBatch, Batch needs to be defined')
  }

  B1_NCells <- ncol(refBatch)
  B2_NCells <- ncol(queBatch)
  Dim_B1 <- kNN*B1_NCells
  Dim_B2 <- kNN*B2_NCells
  B1_B2_NN <- matrix(0L, nrow = Dim_B1, ncol = 2)
  B2_B1_NN <- matrix(0L, nrow = Dim_B2, ncol = 2)
  colnames(B1_B2_NN) <- c("Batch-1", "Batch-2")
  colnames(B2_B1_NN) <- c("Batch-2", "Batch-1")

  NN <- get.knnx(data = t(queBatch), query = t(refBatch), k = kNN)
  NN_Index <- NN$nn.index
  B1_B2_NN[,1] <- rep(c(1:B1_NCells), each = kNN)
  B1_B2_NN[,2] <- t(NN_Index)[1:Dim_B1]

  NN <- get.knnx(data = t(refBatch), query = t(queBatch), k = kNN)
  NN_Index <- NN$nn.index
  B2_B1_NN[,1] <- rep(c(1:B2_NCells), each = kNN)
  B2_B1_NN[,2] <- t(NN_Index)[1:Dim_B2]

  Pairs <- FindMnnPairs(B1_B2_NN = B1_B2_NN, B2_B1_NN = B2_B1_NN, B2_NCells = B2_NCells )

  return(Pairs)

}


#' FindkNN
#'
#' @param m1 Main matrix.
#' @param m2 Matrix to compare.
#' @param k Number of neighbors.
#'
#' @return Data frame with kNN and distances.
FindkNN <- function(m1 = NULL, m2 = NULL,k = NULL){
  nObservations <- nrow(m1)
  df <- data.frame(matrix(data = integer(), ncol = (k+1), nrow = nObservations*k))
  colnames(df) <- c("m1", "m2", "distance")

  NN <- FNN::get.knnx(data = m2, query = m1, k = k)
  df[["m1"]] <- rep(x = seq_len(nObservations), each = k)
  df[["m2"]] <- Reduce(rbind, t(NN$nn.index))
  df[["distance"]] <- Reduce(rbind, t(NN$nn.dist))

  return(df)
}

#' FindMNN
#'
#' @param m1 Main matrix.
#' @param m2 Matrix to compare.
#' @param k Number of neighbors.
#'
#' @return Data frame with MNN and distances.
FindMNN <- function(m1 = NULL, m2 = NULL, k = NULL){

  ## Obtain the crossed kNN pairs
  m1m2NN <- FindkNN(m1 = m1, m2 = m2, k = k)
  m2m1NN <- FindkNN(m1 = m2, m2 = m1, k = k)
  ## Accomodate
  colnames(m2m1NN) <- c("m2","m1", "distance")
  m2m1NN <- m2m1NN[,colnames(m1m2NN)]
  df <- rbind(m1m2NN, m2m1NN)
  ## Find MNN
  df[["pair"]] <- paste(df$m1, df$m2, sep = "_")
  df <- df[duplicated(df$pair),]
  rownames(df) <- df$pair
  df <- df[,c("m1", "m2", "distance")]

  return(df)
}
