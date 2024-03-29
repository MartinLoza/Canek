
#' CorrectBatches
#'
#' Batch-effect correction over a list of single cell batches
#'
#' @param lsBatches List of batches to integrate. Batches should contain the same number of genes as rows.
#' @param queNumCelltypes Number of cell types in the query batch. By default Canek searches the number of
#' cell types using an heuristic algorithm. Change this parameter if you know the number of cell types in advanced.
#' @param sampling Use MNNs pairs sampling when using a Kalman filter to estimate the correction vector.
#' @param numSamples If sampling. Number of MNNs pairs samples to use on the estimation process.
#' @param kNN Number of k-nearest-neighbors used to define the MNNs pairs.
#' @param pcaDim Number of PCA dimensions to use.
#' @param maxMem Maximum number of memberships from the query batch. This parameter is used on the heuristic algorithm to find the number of cell types.
#' @param fuzzy Use fuzzy logic to join the local correction vectors.
#' @param fuzzyPCA Number of PCs to use in the fuzzy process.
#' @param hierarchical Use hierarchical integration scheme when correcting more than two batches.
#' If set to FALSE, the input batches are sorted by number of cells and integrated on descending order.
#' @param verbose Print output.
#' @param estMethod Method to use when estimating the correction vectors:
#' \itemize{
#'   \item{Median. Use the cells median distance}
#'   \item{EKF. Use an extended Kalman filter}
#' }
#' @param pairsFilter Filter MNNs pairs before estimating the correction vectors. If TRUE,
#' the pairs are filtered from outliers using an interquartile range method.
#' @param perCellMNN Threshold value to decide if a membership's correction value is calculated.
#' As a rough interpretation, this values can be thought as the proportion of cells from a membership
#' with an associated MNN pair. If the proportion is low, an specific correction vectors is
#' not calculated for this membership.
#' @param doCosNorm Whether to do cosine normalization.
#' @param fracSampling Fraction of cells to sample in the hierarchical selection (default is NULL, no sampling).
#' @param clusterMethod Method used to identify memberships.
#' @param debug Return correction's information
#' @param ... Pass down methods from RunCanek().
#'
#' @details CorrectBatches is a method to correct batch-effect from two or more single-cell batches.
#' Batch-effects observations are defined using mutual nearest neighbors (MNNs) pairs and cell
#' groups from the query batch are distinguished using clustering. We estimate a correction vector
#' for each cluster using its MNNs pairs and use these vectors to remove the batch effect from the query batch in two ways:
#' \itemize{
#'    \item{A linear correction is performed by equally correcting the cells from the same cluster.}
#'    \item{A non-linear correction is performed by differently correcting each cell using fuzzy logic.}
#' }
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
                           sampling = FALSE, numSamples = NULL,
                           kNN = 30, pcaDim = 50,
                           pairsFilter = FALSE, perCellMNN = 0.08,
                           fuzzy = TRUE, fuzzyPCA = 10,
                           estMethod = "Median", clusterMethod = "louvain",
                           doCosNorm = FALSE, fracSampling = NULL,
                           debug = FALSE, verbose = FALSE, ... ){

  if(debug || verbose){
    tTotal <- Sys.time()
  }

  #Init
  namesInBatches <- names(lsBatches)
  numBatches <- length(lsBatches)
  lsCorrection <- list()
  inCellNames <- lapply(lsBatches, colnames)
  inCellNames <- Reduce(c,inCellNames) #input cell names

  ## Check batches names
  if(is.null(namesInBatches)){
    warning('\nWarning: Batches names not found. Assigning names.', call. = TRUE)
    namesInBatches <- 1:numBatches
    names(lsBatches) <- namesInBatches
  }

  #Check unique cell names
  if(any(duplicated(inCellNames)))
    stop(' Cell names are not unique. Ensure all cell names are unique across batches.', call. = TRUE)

  #Check input batches as matrices
  lsBatches <- lapply(lsBatches, as.matrix)

  # First batch is the one with highest number of cells
  nCells <- sapply(lsBatches, ncol)
  lsBatches <- lsBatches[order(nCells, decreasing = TRUE)]

  # Cosine normalize the input batches
  if (doCosNorm)
    cnBatches <- lapply(lsBatches, CosNorm)
  else
    cnBatches <- lsBatches

  #order <- rep(namesInBatches[1], ncol(lsBatches[[1]]))

  for(i in 2:numBatches){

    namesBatches <- names(lsBatches)

    # hierarchical selection
    if (hierarchical && length(lsBatches) > 2) {

      if(!is.null(fracSampling) && is.numeric(fracSampling)) {

        if(verbose)
          cat(paste("\nHierarchical mode with sampling"))

        sampIdx <- lapply(namesBatches, function(n){
          sample(1:ncol(lsBatches[[n]]),
                 size = floor(ncol(lsBatches[[n]])*fracSampling),
                 replace = FALSE)
        })
        names(sampIdx) <- namesBatches

        nCellsRef <- length(sampIdx[[1]])
      }else{
        nCellsRef <- ncol(lsBatches[[1]])
      }

      # pca with all the other batches
      #nCellsRef <- ncol(lsBatches[[1]])

      # pcaBatches <- lapply(cnBatches[-1], function(x) {
      #   prcomp_irlba(t(cbind(cnBatches[[1]], x)))
      # })

      # Test. Sampling in hierarchical mode
      pcaBatches <- lapply(names(cnBatches)[-1], function(n){
        if (exists("sampIdx")) {
          m <- t(cbind(cnBatches[[1]][, sampIdx[[1]]],
                  cnBatches[[n]][, sampIdx[[n]]]))
        } else {
          m <- t(cbind(cnBatches[[1]], cnBatches[[n]]))
        }
        prcomp_irlba(m)
      })

      nPairs <- lapply(pcaBatches, function(x){
        pairs <- GetMnnPairs(refBatch = t(x$x[1:nCellsRef, ]),
                             queBatch = t(x$x[(nCellsRef+1):nrow(x$x), ]),
                             kNN = 30)$Pairs
        return(nrow(pairs))
      })

      rm(pcaBatches)

      # Select query batch (max number of pairs)
      # We select the first element, in case there are two datasets with the same number of pairs
      Query <- which.max(nPairs) + 1

    } else {  #If the integration is not hierarchical
      Query <- 2
    }

    # Correction
    if(verbose)
      cat(paste('\nINTEGRATING', namesBatches[Query],"INTO", namesBatches[1],"\n", sep = " "))

    Correction <- CorrectBatch(refBatch = lsBatches[[1]], queBatch = lsBatches[[Query]],
                               queNumCelltypes = queNumCelltypes, pcaDim = pcaDim,
                               maxMem = maxMem, kNN = kNN,
                               fuzzy = fuzzy, fuzzyPCA = fuzzyPCA, estMethod = estMethod,
                               pairsFilter = pairsFilter, perCellMNN = perCellMNN,
                               sampling = sampling, numSamples = numSamples,
                               cnRef = cnBatches[[1]], cnQue = cnBatches[[Query]],
                               doCosNorm = doCosNorm, clusterMethod = clusterMethod,
                               verbose = verbose)

    # new ref at the beginning
    lsBatches <- lsBatches[-Query]
    cnBatches <- cnBatches[-Query]
    lsBatches[[1]] <- cbind(lsBatches[[1]], Correction[["Corrected Query Batch"]])
    # new cb at the beginning
    if (doCosNorm)
      cnBatches[[1]] <- CosNorm(lsBatches[[1]])
    else
      cnBatches[[1]] <- lsBatches[[1]]
    names(lsBatches)[1] <- paste(namesBatches[1],namesBatches[Query],sep = "/")

    if(debug == TRUE){
      lsCorrection[[names(lsBatches)[1]]] <- Correction
    }
  }

  #order output dataset
  lsBatches[[1]] <- lsBatches[[1]][,inCellNames]

  if(debug == TRUE){
    lsCorrection[["Batches Integrated"]] <- lsBatches[[1]]
  }

  if(debug || verbose){
    tTotal <- difftime(Sys.time(), tTotal, units = "min")

    if(verbose)
      cat(paste0('\nTotal correction time: ', tTotal, " seconds"))

    if(debug)
      lsCorrection[["Total_Correction_Time"]] <- tTotal
  }

  return(if(debug == FALSE) lsBatches[[1]] else lsCorrection)
}

#' CorrectBatch
#'
#' Batch effect correction on two single-cell batches
#'
#' @param refBatch Reference batch.
#' @param queBatch Query batch (batch to correct).
#' @param queNumCelltypes Number of cell types in the query batch. By default Canek searches the number of cell
#' types using an heuristic algorithm. Change this parameter if you know the number of cell types in advanced.
#' @param sampling Use MNNs pairs sampling when using a Kalman filter to estimate the correction vector.
#' @param numSamples If sampling. Number of MNNs pairs samples to use on the estimation process.
#' @param pairs A numerical matrix containing MNNs pairs cell indexes. First column corresponds to query batch cell indexes.
#' @param idxQuery Numerical vector indicating the index of the cells from the query batch to use
#' on the correction vector estimation.
#' @param idxRef Numerical vector indicating the index of the cells from the reference batch to use
#' on the correction vector estimation.
#' @param kNN Number of k-nearest-neighbors used to define the MNNs pairs.
#' @param pcaDim Number of PCA dimensions to use.
#' @param fuzzyPCA Number of PCs to use in the fuzzy process.
#' @param maxMem Maximum number of memberships from the query batch. This parameter is used on the
#' heuristic algorithm to find the number of cell types.
#' @param fuzzy Use fuzzy logic to join the local correction vectors.
#' @param verbose Print output.
#' @param estMethod Method to use when estimating the correction vectors:
#' \itemize{
#'   \item{Median. Use the cells median distance.}
#'   \item{EKF. Use an extended Kalman filter.}
#' }
#' @param pairsFilter Filter MNNs pairs before estimating the correction vectors. If TRUE,
#' the pairs are filtered from outliers using an interquartile range method.
#' @param perCellMNN Threshold value to decide if a membership's correction value is calculated.
#' As a rough interpretation, this values can be thought as the proportion of cells from a membership
#' with an associated MNN pair. If the proportion is low, an specific correction vectors is
#' not calculated for this membership.
#' @param doCosNorm Whether to do cosine normalization.
#' @param cnRef Cosine normalization of the reference batch.
#' @param cnQue Cosine normalization of the query batch.
#' @param clusterMethod Method used to identify memberships.
#'
#'
#' @details CorrectBatch is a method to correct batch-effect from two single-cell batches.
#' Batch-effects observations are defined using mutual nearest neighbors (MNNs) pairs and cell
#' groups from the query batch are distinguished using clustering. We estimate a correction vector
#' for each cluster using its MNNs pairs and use these vectors to remove the batch effect from the query batch in two ways:
#' \itemize{
#'    \item{A linear correction is performed by equally correcting the cells from the same cluster.}
#'    \item{A non-linear correction is performed by differently correcting each cell using fuzzy logic.}
#' }
#'
#' @examples
#' x <- SimBatches$batches[[1]]
#' y <- SimBatches$batches[[2]]
#' z <- CorrectBatch(x, y)
#' Corrected <- z$`Corrected Query Batch`
#'
#' Uncorrected_PCA <- prcomp(t(cbind(x,y)))
#' plot(Uncorrected_PCA$x[,1:2])
#' Corrected_PCA <- prcomp(t(cbind(x,z$`Corrected Query Batch`)))
#' plot(Corrected_PCA$x[,1:2])
#'
#' @return A list containing the input batches, the corrected query batch, and the correction data
#' @export
#'
CorrectBatch <- function(refBatch, queBatch,
                         cnRef = NULL, cnQue = NULL,
                         queNumCelltypes = NULL, maxMem = 5,
                         pairs = NULL, kNN = 30,
                         sampling = FALSE, numSamples = NULL,
                         idxQuery = NULL, idxRef = NULL,
                         pcaDim = 50, perCellMNN = 0.08,
                         fuzzy = TRUE, fuzzyPCA = 10,
                         estMethod = "Median", clusterMethod = "louvain",
                         pairsFilter = FALSE, doCosNorm = FALSE,
                         verbose = FALSE) {

  tBatch <- Sys.time()

  debugData <- list(info = list(), membership = list())

  memPairs <- NULL
  memCorrData <- list()
  corGene <- NULL
  nMem <- NULL

  if(!is.null(queNumCelltypes)){
    nMem <- queNumCelltypes
  }

  nCellsRef <- ncol(refBatch)
  nCellsQue <- ncol(queBatch)
  nCells <- nCellsRef + nCellsQue

  # FIND MNN pairs ----
  if(is.null(pairs)){
    if (doCosNorm) {
      if (is.null(cnRef)) cnRef <- CosNorm(refBatch)
      if (is.null(cnQue)) cnQue <- CosNorm(queBatch)
    } else {
      cnRef <- refBatch
      cnQue <- queBatch
    }

    pcaBatches <- prcomp_irlba(t(cbind(cnRef, cnQue)), n = pcaDim)

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

  debugData$pairs <- data.frame(ref = colnames(refBatch)[pairs[, 2]], query = colnames(queBatch)[pairs[, 1]])

 if(verbose)
  cat(paste('\n\tNumber of MNN pairs:', nrow(pairs)))

  # FIND memberships ----
  switch(clusterMethod,
    "kmeans" = {
      cluster <- ClusterKMeans(pcaQue[, 1:10], maxMem = maxMem, nMem = nMem, usepam = nCellsQue < 2000, verbose = verbose)
    },
    "louvain" = {
      cluster <- ClusterLouvain(pcaQue[, 1:10], k = kNN, verbose = verbose)
    },
    stop("cluster method unknown.")
  )

  cluMem <- cluster$result
  nMem <- cluster$nMem

 # INIT correction matrix ----
 corGene <- matrix(0, nrow = nrow(refBatch), ncol = nMem)
 colnames(corGene) <- seq_len(nMem)
 zeroCorrection <- rep(FALSE, nMem)

 for(mem in 1:nMem){

   if(verbose)
    cat(paste('\n\nAnalyzing Membership ', mem))

   ## FIND pairs by membership ----
   # Membership cell index
   idxCells <- which(cluMem$cluster == mem)

   debugData[["membership"]][[as.character(mem)]] <- data.frame(cells = colnames(queBatch)[idxCells], membership = mem, corrected = TRUE)

   # Membership cells number
   numCellMem <- ncol(queBatch[,idxCells])

   memPairs <- pairs[pairs[, 1] %in% idxCells, , drop = FALSE]

   # FILTER pairs ----
   if(pairsFilter){

     if (nrow(memPairs) > 10){
       memPairs <- PairsFiltering(refBatch = t(pcaRef), queBatch = t(pcaQue),
                                  pairs = memPairs, verbose = verbose)

        if(verbose)
          cat( paste('\n\tNumber of selected pairs:', nrow(memPairs) ) )
      }else {
        if (verbose)
          cat("\n\n\tNot enough pairs found for membership ", mem, ": no filtering done.")
      }
   }

   # ESTIMATE correction ----
   norNumPairs <- (ceiling(nrow(memPairs)/kNN))/(numCellMem)

   if (norNumPairs > perCellMNN){
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

     if(estMethod == "Mean") {
       if(verbose)
         cat("\n\n\tMean Method")

       corVector <- MeanBE(refBatch = refBatch, queBatch = queBatch,
                             pairs = memPairs)[["Correction Vector"]]
     }

     corGene[,mem] <- corVector

    } else {
      if (verbose)
        cat("\n\n\tNot enough pairs found for membership ", mem, ": no correction vector is computed.")
      corVector <- corGene[,mem]
      zeroCorrection[mem] <- TRUE
      debugData[["membership"]][[as.character(mem)]]$corrected[debugData[["membership"]][[as.character(mem)]]$cells %in% colnames(queBatch)[idxCells]] <- FALSE
    }

   memCorrData[[paste("Membership", mem)]] <- list("Cells Index" = idxCells, "Correction Vector" = corVector)
 }

 # CALCULATE minimum spanning tree ----
 MST <- CalculateMST(cluMem$centers[,1:fuzzyPCA])

 # CHECK No Zero Correction Vectors ----
 isZero <- which(zeroCorrection == TRUE)
 if(length(isZero) == nMem){
   warning('\nWarning: No correction vectors where found.\nConsider using a higher number of kNN or a lower number of clusters to filter pairs', call. = TRUE)
 }else if(length(isZero) != 0){

   noZeroCV <- CheckZeroCV(cluMem = cluMem, corGene = corGene, fuzzyPCA = fuzzyPCA,
                           memCorrData = memCorrData, zeroCorrection = zeroCorrection, MST = MST)

   memCorrData <- noZeroCV$memCorrData
   corGene <- noZeroCV$corGene
   MST <- noZeroCV$MST
   cluMem <- noZeroCV$cluMem
   nMem <- ncol(corGene)
 }


 # FUZZY correction ----
 # Init
 corCell <- matrix(0, nrow = nCellsQue, ncol = ncol(corGene))

 # Set column names according to number of memberships
 colnames(corCell) <- colnames(corGene)

 # Init membership's cells (1 to the cell's membership and 0 to the other memberships)
 for (Mem in colnames(corGene)){
   corCell[which(cluMem$cluster == as.integer(Mem)), Mem] <- 1
 }

 # Fuzzy process and Correction
 if(fuzzy && nMem > 1){

   if(verbose)
     cat('\n\n Fuzzy process ')
   fuzzyData <- Fuzzy(cluMem = cluMem, pcaQue = pcaQue, MST = MST,
                      fuzzyPCA = fuzzyPCA, corCell = corCell,
                      verbose = verbose)

   corCell <- fuzzyData$`Fuzzy Memberships`

 }else{
   fuzzyData <- list("Fuzzy Memberships" = corCell, "MST" = MST,
                     "Fuzzied" =  NULL, "Edges Data" = NULL)
 }

 corMatrix <- (corGene  %*% t(corCell/rowSums(corCell)) )
 queCorrected <-  queBatch + corMatrix

 debugData$matrix$corMatrix <- corMatrix
 debugData$matrix$corCell <- corCell
 debugData$matrix$corGene <- corGene
 debugData$matrix$cells <- colnames(queCorrected)
 debugData$matrix$features <- rownames(queCorrected)

  # SET data lists to return ----
 memData <- list("Cluster Membership" = nMem, "Membership Correction Data" = memCorrData)

 correctionData <- list("Correction Matrix" = corMatrix, "MNN Pairs" = pairs,
                        "Membership Data" = memData, "Fuzzy Data" = fuzzyData,
                        "Clusters" = cluMem)

 tBatch <- difftime(Sys.time(), tBatch, units = "min")


 debugData$info$cputime <- tBatch
 debugData$info$clusterMethod <- clusterMethod
 debugData$info$doCosNorm <- doCosNorm
 debugData$info$fuzzy <- fuzzy
 debugData$info$pairsFilter <- pairsFilter
 debugData$info$pcaQue <- pcaQue
 debugData$info$pcaRef <- pcaRef


 if(verbose)
   cat(paste0('\nBatch correction time: ', tBatch, " seconds"))


 return(list("Reference Batch (B1)" = refBatch, "Query Batch (B2)" = queBatch,
             "Corrected Query Batch"= queCorrected, "Correction Data" = correctionData,
             "Correction_Time" = tBatch, debug = debugData))
}











