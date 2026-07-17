##FindMnnPairs##
#Find MNN pairs given two matrices containing nearest neighbors
# INPUT :
#
# OUTPUT :
FindMnnPairs <- function(B1_B2_NN = NULL, B2_B1_NN = NULL, B2_NCells = NULL ){
  if ( is.null(B1_B2_NN) ){
    stop('B1_B2_NN, Nearest neighbors need to be defined')
  }

  if ( is.null(B2_B1_NN) ){
    stop('B2_B1_NN, Nearest neighbors need to be defined')
  }

  if ( is.null(B2_NCells) ){
    stop('B2_NCells, Number of queBatch cells needs to be defined.')
  }

  # Precompute, which rows of each NN table belong to each query cell
  # index.
  idx_B1_B2 <- split(seq_len(nrow(B1_B2_NN)), factor(B1_B2_NN[,2], levels = seq_len(B2_NCells)))
  idx_B2_B1 <- split(seq_len(nrow(B2_B1_NN)), factor(B2_B1_NN[,1], levels = seq_len(B2_NCells)))

  # Collect each query cell's matched row-blocks in a list and rbind once at
  # the end to optimize running time.
  cellPairs <- vector("list", B2_NCells)

  #CHECK PAIRS
  for (i in 1:B2_NCells) {

    p_B1_B2 <- idx_B1_B2[[i]]
    if(length(p_B1_B2) == 0){
      next()
    }
    B1_B2_sub <- matrix(B1_B2_NN[p_B1_B2,], ncol = 2)

    p_B2_B1 <- idx_B2_B1[[i]]
    if(length(p_B2_B1) == 0){
      next()
    }
    B2_B1_sub <- matrix( B2_B1_NN[p_B2_B1,], ncol = 2)

    d <- dim(B1_B2_sub)[1]
    matchedRows <- vector("list", d)

    for (j in 1:d) {

      m_pair <- which( B2_B1_sub[,2] == B1_B2_sub[j,1] )

      if( length( m_pair ) != 0){
        matchedRows[[j]] <- B2_B1_sub[m_pair, , drop = FALSE]
      }
    }
    cellPairs[[i]] <- do.call(rbind, matchedRows)
  }

  m_Pairs <- do.call(rbind, cellPairs)
  colnames(m_Pairs) <- c('queBatch-Cells-Index', 'refBatch-Cells-Index' )

  return(list("Pairs" = m_Pairs))
}

##Get_MNN_Pairs##
#Get MNN given two batches
# INPUT :
#
# OUTPUT :
GetMnnPairs <- function(refBatch = NULL, queBatch = NULL, kNN = 25, ncores = 1){

  if ( is.null(refBatch) ){
    stop('refBatch, Batch needs to be defined')
  }
  if ( is.null(queBatch) ){
    stop('queBatch, Batch needs to be defined')
  }

  if (ncores > 1) {
    if (.Platform$OS.type == "windows") {
      stop("ncores > 1 is not supported on Windows: parallel::mclapply() requires ",
           "fork(), which Windows doesn't provide. Use ncores = 1.", call. = FALSE)
    }
    if (ncores > parallel::detectCores()) {
      stop("ncores (", ncores, ") exceeds the number of available cores (",
           parallel::detectCores(), ").", call. = FALSE)
    }
  }

  B1_NCells <- ncol(refBatch)
  B2_NCells <- ncol(queBatch)
  Dim_B1 <- kNN*B1_NCells
  Dim_B2 <- kNN*B2_NCells
  B1_B2_NN <- matrix(0L, nrow = Dim_B1, ncol = 2)
  B2_B1_NN <- matrix(0L, nrow = Dim_B2, ncol = 2)
  colnames(B1_B2_NN) <- c("Batch-1", "Batch-2")
  colnames(B2_B1_NN) <- c("Batch-2", "Batch-1")

  # The two directional searches are independent of each other, so with
  # ncores > 1 they run concurrently (only 2 tasks exist here, so requesting
  # more than 2 cores doesn't speed up this step any further).
  if (ncores > 1) {
    NNs <- parallel::mclapply(1:2, function(i) {
      if (i == 1) get.knnx(data = t(queBatch), query = t(refBatch), k = kNN)
      else get.knnx(data = t(refBatch), query = t(queBatch), k = kNN)
    }, mc.cores = min(ncores, 2))
    NN_B1B2 <- NNs[[1]]
    NN_B2B1 <- NNs[[2]]
  } else {
    NN_B1B2 <- get.knnx(data = t(queBatch), query = t(refBatch), k = kNN)
    NN_B2B1 <- get.knnx(data = t(refBatch), query = t(queBatch), k = kNN)
  }

  NN_Index <- NN_B1B2$nn.index
  B1_B2_NN[,1] <- rep(c(1:B1_NCells), each = kNN)
  B1_B2_NN[,2] <- t(NN_Index)[1:Dim_B1]

  NN_Index <- NN_B2B1$nn.index
  B2_B1_NN[,1] <- rep(c(1:B2_NCells), each = kNN)
  B2_B1_NN[,2] <- t(NN_Index)[1:Dim_B2]

  Pairs <- FindMnnPairs(B1_B2_NN = B1_B2_NN, B2_B1_NN = B2_B1_NN, B2_NCells = B2_NCells )

  return(Pairs)

}
