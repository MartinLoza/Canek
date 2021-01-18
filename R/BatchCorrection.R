
#' CorrectBatches
#'
#' Batch-Effect correction over a list of single cell batches
#'
#' @param lsBatches List of batches to integrate. Batches should contain the same number of genes as rows.
#' @param queNumCelltypes A number indicating the expected number of cells types on the batches to integrate. The default value is set as a string "Surprise-me" on which an estimation of the cell types is defined.
#' @param sampling Whether or not sampling of MNNs pairs is used on the estimation process.
#' @param numSamples Number of MNNs pairs samples used on the estimation process.
#' @param kNN Number of k-nearest-neighbors used to find MNNs pairs.
#' @param pcaDim PCA dimensions used to find MNNs pairs.
#' @param maxMem Maximum number of memberships used when memberships are automatically defined.
#' @param fuzzy Whether or not a fuzzy logic join is used on the local correction vectors.
#' @param hierarchical Whether or not a hiearchical integration scheme is used when correcting more than two batches.
#' @param verbose Print output
#' @param estMethod TODO
#' @param pairsFilter whether to perform pair filtering (default: FALSE)
#' @param perCellMNN TODO
#' @param debug TODO
#' @param ... pass down methods from RunCanek().
#'
#' @details CorrectBatches is non-linear/linear hybrid method for single-cell batch-effect correction that couples identification of similar cells
#'  between datasets using Mutual Nearest Neighbors (MNNs) with an Extended Kalman Filter (EKF).
#'
#'  A non-linear correction is performed using fuzzy logic to join a set of linear correction vectors which are cell-type locally estimated.
#'
#' @examples
#' Batches <- SimBatches$batches
#' z <- CorrectBatches(Batches)
#'
#' Uncorrected_PCA <- prcomp(t(cbind(Batches[[1]], Batches[[2]])))
#' plot(Uncorrected_PCA$x[,1:2])
#' Corrected_PCA <- prcomp(t(z))
#' plot(Corrected_PCA$x[,1:2])
#'
#' @return A list containing the integrated datasets as matrix and the correction data .
#' @export
#'
CorrectBatches <- function(lsBatches, hierarchical = TRUE,
                           queNumCelltypes = NULL, maxMem = 5,
                           sampling = NULL, numSamples = NULL,
                           kNN = 30, pcaDim = 50,
                           pairsFilter = FALSE, perCellMNN = 0.08,
                           fuzzy = TRUE, estMethod = "Median",
                           debug = FALSE, verbose = FALSE, ... ){

  if(verbose)
    tic("\nTotal correction time ")

  #Init
  namesInBatches <- names(lsBatches)
  numBatches <- length(lsBatches)
  lsCorrection <- list()

  #Check input batches as matrices
  lsBatches <- lapply(lsBatches, as.matrix)

  # First batch is the one with highest number of cells
  nCells <- as.data.frame(lapply(lsBatches, ncol))
  lsBatches <- lsBatches[sort(t(nCells), decreasing = TRUE, index.return = TRUE)$ix]

  # Cosine normalize the input batches
  cnBatches <- lapply(lsBatches, batchelor::cosineNorm)

  order <- rep(namesInBatches[1], ncol(lsBatches[[1]]))

  for(i in 2:numBatches){

    namesBatches <- names(lsBatches)

    # hierarchical selection
    if(hierarchical == TRUE & length(lsBatches) > 2){

      # pca with all the other batches
      nCellsRef <- ncol(lsBatches[[1]])

      pcaBatches <- lapply(cnBatches[-1], function(x){prcomp_irlba(t(cbind(cnBatches[[1]],x)))})

      nPairs <- lapply(pcaBatches, function(x){
        pairs <- GetMnnPairs(refBatch = t(x$x[1:nCellsRef,]),
                             queBatch = t(x$x[(nCellsRef+1):nrow(x$x),]),
                             kNN = 30)$Pairs
        return(nrow(pairs))
      })

      rm(pcaBatches)

      # Select query batch (max number of pairs)
      nPairs <- as.data.frame(nPairs)
      Query <- (1+which(nPairs == max(nPairs, na.rm = TRUE)))

    }else{  #If the integration is not hierarchical
      Query <- 2
    }

    order <- c(order, rep(namesBatches[Query], ncol(lsBatches[[Query]])))

    # Correction
    if(verbose)
      cat(paste('\nINTEGRATING', namesBatches[Query],"INTO", namesBatches[1],"\n", sep = " "))

    Correction <- CorrectBatch(refBatch = lsBatches[[1]], queBatch = lsBatches[[Query]],
                               queNumCelltypes = queNumCelltypes, pcaDim = pcaDim,
                               maxMem = maxMem, kNN = kNN,
                               fuzzy = fuzzy, estMethod = estMethod,
                               pairsFilter = pairsFilter, perCellMNN = perCellMNN,
                               sampling = sampling, numSamples = numSamples,
                               cnRef = cnBatches[[1]], cnQue = cnBatches[[Query]],
                               verbose = verbose)

    # new ref at the beggining
    lsBatches <- lsBatches[-Query]
    cnBatches <- cnBatches[-Query]
    lsBatches[[1]] <- cbind(lsBatches[[1]], Correction[["Corrected Query Batch"]])
    # new cb at the beggining
    cnBatches[[1]] <- batchelor::cosineNorm(lsBatches[[1]])
    names(lsBatches)[1] <- paste(namesBatches[1],namesBatches[Query],sep = "/")

    if(debug == TRUE){
      lsCorrection[[names(lsBatches)[1]]] <- Correction
    }
  }

  if(debug == TRUE){
    lsCorrection[["Batches Integrated"]] <- lsBatches[[1]]
  }

  #order output dataset
  outOrder <- integer()
  for(i in namesInBatches){
    outOrder <- c(outOrder, which(order == i))
  }

  lsBatches[[1]] <- lsBatches[[1]][,outOrder]

  if(verbose)
    toc()

  return(if(debug == FALSE) lsBatches[[1]] else lsCorrection)
}

#' CorrectBatch
#'
#' Function to correct batch effect over two batches
#'
#' @param refBatch Batch to use as reference for the integration.
#' @param queBatch Batch to correct.
#' @param queNumCelltypes A number indicating the expected number of cells types on the batches to integrate. The default value is set as a string "Surprise-me" on which an estimation of the cell types is defined.
#' @param sampling Whether or not samples MNNs pairs' samples are used on the estimation process.
#' @param numSamples Number of MNNs pairs' samples used on the estimation process.
#' @param pairs A matrix containing MNNs pairs. First column corresponds to query-batch cell indexes.
#' @param idxQuery Index of cells from the query-batch used as observations of the batch-effect.
#' @param idxRef Index of cells from the reference-batch used as observations of the batch-effect.
#' @param kNN Number of k-nearest-neighbors used to find MNNs pairs.
#' @param pcaDim PCA dimensions used to find MNNs pairs.
#' @param maxMem Maximum number of memberships used when memberships are automatically defined.
#' @param fuzzy Whether or not a fuzzy logic join is used on the local correction vectors.
#' @param verbose Print output
#' @param estMethod TODO
#' @param pairsFilter whether to perform pair filtering (default: FALSE)
#' @param perCellMNN TODO
#' @param cnRef TODO
#' @param cnQue TODO
#'
#' @details Canek, a new non-linear/linear hybrid method for batch-effect correction that couples identification of similar cells
#'  between datasets using Mutual Nearest Neighbors (MNNs) with an Extended Kalman Filter (EKF).
#'
#'  A non-linear correction is performed by using fuzzy logic to join a set of linear correction vectors which are cell-type locally estimated.
#'
#' @examples
#' x <- SimBatches$batches[[1]]
#' y <- SimBatches$batches[[2]]
#' z <- Canek:::CorrectBatch(x,y)
#' Corrected <- z$`Corrected Query Batch`
#'
#' Uncorrected_PCA <- prcomp(t(cbind(x,y)))
#' plot(Uncorrected_PCA$x[,1:2])
#' Corrected_PCA <- prcomp(t(cbind(x,z$`Corrected Query Batch`)))
#' plot(Corrected_PCA$x[,1:2])
#'
#' @return A list containing the corrected batch as a matrix and correction data
#' @export
#'
#'
CorrectBatch <- function(refBatch, queBatch,
                         cnRef = NULL, cnQue = NULL,
                         queNumCelltypes = NULL, maxMem = 5,
                         pairs = NULL, kNN = 30,
                         sampling = NULL, numSamples = NULL,
                         idxQuery = NULL, idxRef = NULL,
                         pcaDim = 50, perCellMNN = 0.08,
                         fuzzy = TRUE, estMethod = "Median",
                         pairsFilter = FALSE, verbose = FALSE){

  if(verbose)
    tic("\n Correction time")

  memPairs <- NULL
  memCorrData <- list()
  corGene <- NULL
  fuzzyData <- NULL
  nMem <- NULL

  if(!is.null(queNumCelltypes)){
    nMem <- queNumCelltypes
  }

  nCellsRef <- ncol(refBatch)
  nCellsQue <- ncol(queBatch)
  nCells <-nCellsRef + nCellsQue

  if(is.null(pairs)){

    pcaBatches <- prcomp_irlba(t(cbind(if(!is.null(cnRef)) cnRef else batchelor::cosineNorm(refBatch),
                                       if(!is.null(cnQue)) cnQue else batchelor::cosineNorm(queBatch))),
                               n = pcaDim)

    pcaRef <- pcaBatches$x[1:nCellsRef,]
    pcaQue <- pcaBatches$x[(nCellsRef+1):nCells,]

    rm(pcaBatches)

    if(verbose)
      cat(paste("\n\nFinding mutual nearest neighbors from", kNN,"nearest neighbors"))

    pairs <- GetMnnPairs(refBatch = if(is.null(idxRef)) t(pcaRef) else t(pcaRef[,idxRef]),
                         queBatch = if(is.null(idxQuery)) t(pcaQue) else t(pcaRef[,idxQuery]),
                         kNN = kNN)

    pairs <- pairs$Pairs
  }else{

    pcaQue <- prcomp_irlba(t(queBatch),n = 10)
    pcaQue <- pcaQue$x
  }

 if(verbose)
  cat(paste('\n\tNumber of MNN pairs:', nrow(pairs)))

 if(is.null(nMem)){

   if(verbose)
    cat("\n\nFinding number of memberships")

   nMem <- pamk(pcaQue[,1:10],
                krange = 1:maxMem,
                usepam = (if(nCellsQue < 2000) TRUE else FALSE))$nc

   if(verbose)
    cat(paste('\n\tNumber of memberships found:', nMem) )
 }

 #Cluster in memberships
 cluMem <- kmeans(pcaQue[,1:10],nMem)

 #INIT Correction Matrix
 corGene <- matrix(0, nrow = nrow(refBatch), ncol = nMem)

 zeroCorrection <- rep(FALSE, nMem)

 for(mem in 1:nMem){

   if(verbose)
    cat(paste('\n\nAnalyzing Membership ', mem) )

   #Membership cell index
   idxCells <- which(cluMem$cluster == mem)
   #Membership cells number
   numCellMem <- ncol(queBatch[,idxCells])

   #########################
   ###Pairs by membership###
   #########################
   memPairs <- integer()
   for (j in 1:length(idxCells)) {
     memPairs <- c(memPairs,which(pairs[,1] == idxCells[j]))
   }
   #Subset of pairs corresponding to the membership
   memPairs <- pairs[memPairs,]

   #########################
   ###Pairs by clustering###
   #########################
   if(pairsFilter){

     if (nrow(memPairs) > 10){
       memPairs <- PairsFiltering(refBatch = t(pcaRef), queBatch = t(pcaQue),
                                  pairs = memPairs, verbose = verbose)

        if(verbose)
          cat( paste('\n\tNumber of selected pairs:', nrow(memPairs) ) )
      }else {
        warning('\nWarning: Cannot perform pairs filtering due to a low number of pairs from this Membership', call. = TRUE)
      }
   }

   norNumPairs <- (ceiling(nrow(memPairs)/kNN))/(numCellMem)

   if (norNumPairs > perCellMNN){

      ####################################
      ###Estimation of the batch effect###
      ####################################

     if(estMethod == "EKF"){
       if(verbose)
         cat("\n\n\tEXTENDED KALMAN FILTER")

       corVector <- EkfBE(refBatch = refBatch, queBatch = queBatch,
                          pairs = memPairs, sampling = sampling,
                          numSamples = numSamples, verbose = verbose)[["Correction Vector"]]
     }

     if(estMethod == "Median"){
       if(verbose)
         cat("\n\n\tMedian Method")

       corVector <- MedianBE(refBatch = refBatch, queBatch = queBatch,
                             pairs = memPairs)[["Correction Vector"]]
     }

     corGene[,mem] <- corVector

    }else{
      warning('\nWarning: Not enough pairs found for this Membership. No correction is performed', call. = TRUE)
      corVector <- corGene[,mem]
      zeroCorrection[mem] <- TRUE
    }

   memCorrData[[paste("Membership", mem)]] <- list("Cells Index" = idxCells, "Correction Vector" = corVector)
 }

 #################
 ###   FUZZY   ###
 #################

 ####INIT####
 corCell <- matrix(0, nrow = nCellsQue, ncol = nMem )

 # Set column names according to number of memberships
 colnames(corCell) <- paste0("Mem-", seq_len(nMem))

 # Init membership's cells (1 to the cell's membership and 0 to the other memberships)
 for (Mem in seq_len(nMem)){
   corCell[which(cluMem$cluster == Mem), Mem] <- 1
 }

 #fuzzy process and Correction
 if(fuzzy == TRUE & nMem > 1){

   if(verbose)
    cat('\n\nFUZZY ')

   fuzzyData <- Fuzzy(cluMem = cluMem, pcaQue = pcaQue,
                      corCell = corCell, verbose = verbose)

   corCell <- fuzzyData$`Fuzzy Memberships`
   MST <- fuzzyData$MST

 }else{
   MST <- mst(dist(cluMem$centers[,1:2]))
 }

 #No Zero Correction Vectors
 isZero <- which(zeroCorrection == TRUE)
 if(length(isZero) == nMem){
   warning('\nWarning: No correction vectors where found.\nConsider using a higher number of kNN or a lower number of clusters to filter pairs', call. = TRUE)
 }else if(length(isZero) != 0){

   noZeroCV <- CheckZeroCV(MST = MST, cluMem = cluMem,
                           memCorrData = memCorrData, corGene = corGene,
                           zeroCorrection = zeroCorrection)

   memCorrData <- noZeroCV[["memCorrData"]]
   corGene <- noZeroCV[["Correction_Matrix"]]
 }

 queCorrected <-  queBatch + (corGene  %*% t(corCell/rowSums(corCell)) )

 ### Set data lists to return
 memData <- list("Cluster Membership" = nMem, "Membership Correction Data" = memCorrData)

 correctionData <- list("Correction Matrix" = corGene, "MNN Pairs" = pairs,
                        "Membership Data" = memData, "Fuzzy Data" = fuzzyData)

 if(verbose)
   toc()

 return(list("Reference Batch (B1)" = refBatch, "Query Batch (B2)" = queBatch,
             "Corrected Query Batch"= queCorrected, "Correction Data" = correctionData))
}










