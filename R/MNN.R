##Get_Neighbors##
#Function to find nearest neighbors
# INPUT :
#
# OUTPUT :
Get_Neighbors <- function( Crossed_Distances = NULL, B1_NCells = NULL, B2_NCells = NULL, k_Neighbors = 20 ){

  if ( (B2_NCells < k_Neighbors) || (B1_NCells < k_Neighbors) ){
    stop('Number of cells is lower than k-Neighbors')
  }

  B1_B2_NN <- NULL

  for (i in 1:B1_NCells) {
    kpairs <- order(Crossed_Distances[i,1:B2_NCells])

    for (j in 1:k_Neighbors) {
      B1_B2_NN <- rbind(B1_B2_NN, c(i,kpairs[j]))
    }
  }

  colnames(B1_B2_NN) <- c("V1", "V2")
  return(list("NN vector" = B1_B2_NN))
}

##Find_MNN_Pairs##
#Find MNN pairs given two matrices containing nearest neighbors
# INPUT :
#
# OUTPUT :
Find_MNN_Pairs <- function(B1_B2_NN = NULL, B2_B1_NN = NULL, B2_NCells = NULL ){
  # INIT
  m_Pairs <- NULL

  if ( is.null(B1_B2_NN) ){
    stop('B1_B2_NN, Nearest neighbors need to be defined')
  }
  if ( is.null(B2_B1_NN) ){
    stop('B2_B1_NN, Nearest neighbors need to be defined')
  }

  if( is.null(B2_NCells) ){
    B2_NCells <- dim(B2)[2]
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

  colnames(m_Pairs) <- c('B2-Cells-Index', 'B1-Cells-Index' )

  return(list("Pairs" = m_Pairs))
}

##Get_MNN_Pairs##
#Get MNN given two batches
# INPUT :
#
# OUTPUT :
Get_MNN_Pairs <- function(B1 = NULL, B2 = NULL, k_Neighbors = 20){

  if ( is.null(B1) ){
    stop('B1, Batch needs to be defined')
  }
  if ( is.null(B2) ){
    stop('B2, Batch needs to be defined')
  }

  B1_NCells <- dim(B1)[2]
  B2_NCells <- dim(B2)[2]

  NN <- get.knnx(data = t(B2), query = t(B1), k = k_Neighbors)
  NN_Index <- NN$nn.index

  B1_B2_NN <- NULL
  for (i in 1:B1_NCells) {

    for (j in 1:k_Neighbors) {
      B1_B2_NN <- rbind( B1_B2_NN, c(i,NN_Index[i,j]) )
    }

  }
  colnames(B1_B2_NN) <- c("V1", "V2")

  NN <- get.knnx(data = t(B1), query = t(B2), k = k_Neighbors)
  NN_Index <- NN$nn.index

  B2_B1_NN <- NULL
  for (i in 1:B2_NCells) {

    for (j in 1:k_Neighbors) {
      B2_B1_NN <- rbind( B2_B1_NN, c(i,NN_Index[i,j]) )
    }

  }
  colnames(B2_B1_NN) <- c("V1", "V2")

  Pairs <- Find_MNN_Pairs(B1_B2_NN = B1_B2_NN, B2_B1_NN = B2_B1_NN, B2_NCells = B2_NCells )

  return(Pairs)

}
