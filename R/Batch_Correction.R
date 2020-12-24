
#' Correct_Batches
#'
#' Batch-Effect correction over a list of single cell batches
#'
#' @param Batches List of batches to integrate. Batches should contain the same number of genes as rows.
#' @param queNumCelltypes A number indicating the expected number of cells types on the batches to integrate. The default value is set as a string "Surprise-me" on which an estimation of the cell types is defined.
#' @param Sampling Whether or not sampling of MNNs pairs is used on the estimation process.
#' @param Number_Samples Number of MNNs pairs samples used on the estimation process.
#' @param k_Neighbors Number of k-nearest-neighbors used to find MNNs pairs.
#' @param Dimensions PCA dimensions used to find MNNs pairs.
#' @param Max_Membership Maximum number of memberships used when memberships are automatically defined.
#' @param Fuzzy Whether or not a fuzzy logic join is used on the local correction vectors.
#' @param Hierarchical Whether or not a hiearchical integration scheme is used when correcting more than two batches.
#' @param Verbose Print output
#' @param Estimation TODO
#' @param FilterPairs whether to perform pair filtering (default: FALSE)
#' @param perCellMNN TODO
#' @param ... pass down methods from RunCanek().
#'
#' @details Correct_Batches is non-linear/linear hybrid method for single-cell batch-effect correction that couples identification of similar cells
#'  between datasets using Mutual Nearest Neighbors (MNNs) with an Extended Kalman Filter (EKF).
#'
#'  A non-linear correction is performed using fuzzy logic to join a set of linear correction vectors which are cell-type locally estimated.
#'
#' @examples
#' Batches <- SimBatches$batches
#' Corrected <- Correct_Batches(Batches)
#' #' z <- Correct_Batches(Batches)
#' #Corrected <- z$`Batches Integrated`
#'
#' Uncorrected_PCA <- prcomp(t(cbind(Batches[[1]], Batches[[2]])))
#' plot(Uncorrected_PCA$x[,1:2])
#' Corrected_PCA <- prcomp(t(Corrected))
#' # Corrected_PCA <- prcomp(t(z$`Batches Integrated`))
#' plot(Corrected_PCA$x[,1:2])
#'
#' @return A list containing the integrated datasets as matrix and the correction data .
#' @export
#'
Correct_Batches <- function(Batches, queNumCelltypes = NULL,
                            Sampling = NULL,
                            Number_Samples = NULL,
                            k_Neighbors = 30,
                            #PCA = TRUE,
                            Dimensions = 50,
                            Max_Membership = 5,
                            Fuzzy = TRUE,
                            Hierarchical = TRUE,
                            Verbose = FALSE,
                            #Cosine_Norm = TRUE,
                            Estimation = "Median",
                            FilterPairs = FALSE,
                            perCellMNN = 0.08,
                            ...
                            ){

  if(Verbose)
    tic("\nTotal correction time ")

  #Init
  namesInBatches <- names(Batches)
  Num_Batches <- length(Batches)

  #Check input batches as matrices
  Batches <- lapply(Batches, as.matrix)

  #order first batch is the one with highest number of cells
  nCells <- as.data.frame(lapply(Batches, ncol))
  Batches <- Batches[sort(t(nCells), decreasing = TRUE, index.return = TRUE)$ix]

  # Cosine normalize the input batches
  cnBatches <- lapply(Batches, batchelor::cosineNorm)

  #In hierarchical mode, the pairs of the batches are checked in order to decide which batches are integrated first. The logic
  #is that more similar batches would share a higher number of pairs
  if(Hierarchical == TRUE & Num_Batches >2 ){

    Order <- rep(namesInBatches[1], ncol(Batches[[1]]))

    for(i in 2:Num_Batches){
      # ref at the beggining


      # Get pairs
      if(length(Batches) > 2){

        # pca with all the other batches
        nCellsRef <- ncol(Batches[[1]])

        pcaBatches <- lapply(cnBatches[-1], function(x){prcomp_irlba(t(cbind(cnBatches[[1]],x)))})

        nPairs <- lapply(pcaBatches, function(x){
          pairs <- Get_MNN_Pairs(B1 = t(x$x[1:nCellsRef,]),
                                 B2 = t(x$x[(nCellsRef+1):nrow(x$x),]),
                                 k_Neighbors = 30)$Pairs
          return(nrow(pairs))
        })

        rm(pcaBatches)

        # Select query batch (max number of pairs)
        nPairs <- as.data.frame(nPairs)
        Query <- (1+which(nPairs == max(nPairs, na.rm = TRUE)))
      }else{
        Query <- 2
      }

      Order <- c(Order, rep(names(Batches)[Query], ncol(Batches[[Query]])))

      # Eliminate unnecesary data
      cnBatches <- cnBatches[-Query]

      # correct
      if(Verbose)
        cat(paste('\nINTEGRATING', names(Batches)[Query],"INTO", names(Batches)[1],"\n", sep = " "))

      Correction <- Correct_Batch(refBatch = Batches[[1]], queBatch = Batches[[Query]],
                                  queNumCelltypes = queNumCelltypes, Dimensions = Dimensions,
                                  Max_Membership = Max_Membership, k_Neighbors = k_Neighbors,
                                  Fuzzy = Fuzzy, Estimation = Estimation,
                                  FilterPairs = FilterPairs, perCellMNN = perCellMNN,
                                  Sampling = Sampling, Number_Samples = Number_Samples,
                                  #cnRef = cnBatches[[Ref]], cnQue = cnBatches[[Query]],
                                  Verbose = Verbose)[["Corrected Query Batch"]]

      # new ref at the beggining
      Batches <- Batches[-Query]
      Batches[[1]] <- cbind(Batches[[1]], Correction)
      # new cb at the beggining
      cnBatches[[1]] <- batchelor::cosineNorm(Batches[[1]])
      names(Batches)[1] <- paste(names(Batches)[1],namesInBatches[Query],sep = "/")
      # repeat
    }

    #order output dataset
    mOrder <- integer()
    for(i in namesInBatches){
      mOrder <- c(mOrder, which(Order == i))
    }

    Corrected_Batches[["Batches Integrated"]] <- Batches[[1]][,mOrder]



  }else{  #If the integration is not hierarchical

    for(i in 2:Num_Batches){


      if(Verbose)
        cat(paste('\nINTEGRATING', namesInBatches[Query],"INTO", namesInBatches[Ref],"\n", sep = " ") )

      Correction <- Correct_Batch(refBatch = Batches[[Ref]], queBatch = Batches[[Query]],
                                  queNumCelltypes = queNumCelltypes, Dimensions = Dimensions,
                                  Max_Membership = Max_Membership, k_Neighbors = k_Neighbors,
                                  Fuzzy = Fuzzy, Estimation = Estimation,
                                  FilterPairs = FilterPairs, perCellMNN = perCellMNN,
                                  Sampling = Sampling, Number_Samples = Number_Samples,


    }

    Corrected_Batches[["Batches Integrated"]] <- Batches[[1]]
  }

  if(Verbose)
    toc()

  return(Corrected_Batches)
}

#' Correct_Batch
#'
#' Function to correct batch effect over two batches
#'
#' @param refBatch Batch to use as reference for the integration.
#' @param queBatch Batch to correct.
#' @param queNumCelltypes A number indicating the expected number of cells types on the batches to integrate. The default value is set as a string "Surprise-me" on which an estimation of the cell types is defined.
#' @param Sampling Whether or not samples MNNs pairs' samples are used on the estimation process.
#' @param Number_Samples Number of MNNs pairs' samples used on the estimation process.
#' @param Pairs A matrix containing MNNs pairs. First column corresponds to query-batch cell indexes.
#' @param Cells_Index_Query Index of cells from the query-batch used as observations of the batch-effect.
#' @param Cells_Index_Reference Index of cells from the reference-batch used as observations of the batch-effect.
#' @param k_Neighbors Number of k-nearest-neighbors used to find MNNs pairs.
#' @param Dimensions PCA dimensions used to find MNNs pairs.
#' @param Max_Membership Maximum number of memberships used when memberships are automatically defined.
#' @param Fuzzy Whether or not a fuzzy logic join is used on the local correction vectors.
#' @param Verbose Print output
#' @param Estimation TODO
#' @param FilterPairs whether to perform pair filtering (default: FALSE)
#' @param perCellMNN TODO
#'
#' @details Canek, a new non-linear/linear hybrid method for batch-effect correction that couples identification of similar cells
#'  between datasets using Mutual Nearest Neighbors (MNNs) with an Extended Kalman Filter (EKF).
#'
#'  A non-linear correction is performed by using fuzzy logic to join a set of linear correction vectors which are cell-type locally estimated.
#'
#' @examples
#' x <- SimBatches$batches[[1]]
#' y <- SimBatches$batches[[2]]
#' z <- Canek:::Correct_Batch(x,y)
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
Correct_Batch <- function(refBatch, queBatch,
                          queNumCelltypes = NULL,
                          Sampling = NULL,
                          Number_Samples = NULL,
                          Pairs = NULL,
                          Cells_Index_Query = NULL,
                          Cells_Index_Reference = NULL,
                          k_Neighbors = 30,
                          Dimensions = 50,
                          Max_Membership = 5,
                          Fuzzy = TRUE,
                          Verbose = FALSE,
                          #PCA = TRUE, #TODO:implement for different input data
                          #Cosine_Norm = TRUE, # TODO: implement for different input data
                          Estimation = "Median",
                          FilterPairs = FALSE,
                          perCellMNN = 0.08,
                          cnRef = NULL, cnQue = NULL
                          ){

  if(Verbose)
    tic("\n Correction time")

  memPairs <- NULL
  Membership_Correction_Data <- list()
  corGene <- NULL
  Fuzzy_Data <- NULL
  nMem <- NULL

  if(!is.null(queNumCelltypes)){
    nMem <- queNumCelltypes
  }

  # TODO: correguir esto para no crear nuevos datasets, solo crear un indice. USAR INDICE SOLO EN ENCONTRAR MNN

  # if( !is.null(Cells_Index_Reference) ){
  #   if ( nMem > 1 ){
  #     warning('\nWarning: CELLS INDEX NOT USED. Cannot use cells index for more than one membership function', call. = TRUE)
  #     B1_Selected <- refBatch
  #   } else{
  #     B1_Selected <- refBatch[,Cells_Index_Reference]
  #   }
  # } else{
  #   B1_Selected <- refBatch
  # }
  #
  # if( !is.null(Cells_Index_Query) ){
  #   if ( nMem > 1 ){
  #     warning('\nWarning: CELLS INDEX NOT USED. Cannot use cells index for more than one membership function', call. = TRUE)
  #     B2_Selected <- queBatch
  #   }else{
  #     B2_Selected <- queBatch[,Cells_Index_Query]
  #   }
  # }else{
  #   B2_Selected <- queBatch
  # }

  nCellsRef <- ncol(refBatch)
  nCellsQue <- ncol(queBatch)
  Num_Cells <-nCellsRef + nCellsQue

  if(is.null(Pairs)){

    PCA_Batches <- prcomp_irlba(t(cbind(if(!is.null(cnRef)) cnRef else batchelor::cosineNorm(refBatch),
                                        if(!is.null(cnQue)) cnQue else batchelor::cosineNorm(queBatch))),
                                n = Dimensions)

    pcaRef <- PCA_Batches$x[1:nCellsRef,]
    pcaQue <- PCA_Batches$x[(nCellsRef+1):Num_Cells,]

    rm(PCA_Batches)

    if(Verbose)
      cat(paste("\n\nFinding mutual nearest neighbors from", k_Neighbors,"nearest neighbors"))

    Pairs <- Get_MNN_Pairs(B1 = t(pcaRef),
                           B2 = t(pcaQue),
                           k_Neighbors = k_Neighbors)

    Pairs <- Pairs$Pairs
  }else{

    pcaQue <- prcomp_irlba(t(queBatch),n = 10)
    pcaQue <- pcaQue$x
  }

 if(Verbose)
  cat(paste( '\n\tNumber of MNN pairs:', nrow(Pairs) ))

 if(is.null(nMem)){

   if(Verbose)
    cat("\n\nFinding number of memberships")

   nMem <- pamk(pcaQue[,1:10],
                krange = 1:Max_Membership,
                usepam = (if(nCellsQue < 2000) TRUE else FALSE))$nc

   if(Verbose)
    cat(paste('\n\tNumber of memberships found:', nMem) )
 }

 #Cluster in memberships
 cluMem <- kmeans(pcaQue[,1:10],nMem)

 #INIT Correction Matrix
 corGene <- matrix(0, nrow = nrow(refBatch), ncol = nMem)

 Zero_Correction <- rep(FALSE, nMem)

 for(Membership in 1:nMem){

   if(Verbose)
    cat(paste('\n\nAnalyzing Membership ', Membership) )

   #Membership cell index
   Membership_Cells_Index <- which(cluMem$cluster == Membership)
   #Membership cells number
   numCellMembership <- ncol(queBatch[,Membership_Cells_Index])

   #########################
   ###Pairs by membership###
   #########################
   Membership_Pairs_Index <- integer()
   for (j in 1:length(Membership_Cells_Index)) {
     Membership_Pairs_Index <- c(Membership_Pairs_Index,which(Pairs[,1] == Membership_Cells_Index[j]))
   }
   #Subset of pairs corresponding to the membership
   memPairs <- Pairs[Membership_Pairs_Index, ]

   #########################
   ###Pairs by clustering###
   #########################
   if(FilterPairs){

     if (nrow(memPairs) > 10){
       memPairs <- Pairs_Selection(B1 = t(pcaRef), B2 = t(pcaQue),
                                   Pairs = memPairs, Verbose = Verbose)

        if(Verbose)
          cat( paste( '\n\tNumber of selected pairs:', nrow(memPairs) ) )
      }else {
        warning('\nWarning: Not enough pairs found for this Membership.\nNo pairs selection is performed', call. = TRUE)
      }
   }

   norNumPairs <- (ceiling(nrow(memPairs)/k_Neighbors))/(numCellMembership)

   if (norNumPairs > perCellMNN){

      ####################################
      ###Estimation of the batch effect###
      ####################################

     if(Estimation == "EKF"){
       if(Verbose)
         cat("\n\n\tEXTENDED KALMAN FILTER")

       corVector <- EKF_BE(B1 = refBatch, B2 = queBatch,
                           Pairs = memPairs, Sampling = Sampling,
                           Number_Samples = Number_Samples, Verbose = Verbose)[["Correction Vector"]]
     }

     if(Estimation == "Median"){
       if(Verbose)
         cat("\n\n\tMedian Method")

       corVector <- Average_BE(B1 = refBatch, B2 = queBatch,
                               Pairs = memPairs)[["Correction Vector"]]
     }

     corGene[,Membership] <- corVector

    }else{
      warning('\nWarning: Not enough pairs found for this Membership. No correction is performed', call. = TRUE)
      corVector <- corGene[,Membership]
      Zero_Correction[Membership] <- TRUE
    }

   Membership_Correction_Data[[paste("Membership", Membership)]] <- list("Cells Index" = Membership_Cells_Index,
                                                                         "Correction Vector" = corVector)
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

 #Fuzzy process and Correction
 if(Fuzzy == TRUE & nMem > 1){

   if(Verbose)
    cat('\n\nFUZZY ')

   Fuzzy_Data <- Fuzzy(cluMem = cluMem, Cells_PCA = pcaQue,
                       corCell = corCell, Verbose = Verbose)

   corCell <- Fuzzy_Data$`Fuzzy Memberships`
   MST <- Fuzzy_Data$MST

 }else{
   MST <- mst(dist(cluMem$centers[,1:2]))
 }

 #No Zero Correction Vectors
 Is_Zero <- which(Zero_Correction == TRUE)
 if(length(Is_Zero) == nMem){
   warning('\nWarning: No correction vectors where found.\nConsider using a higher number of kNN or a lower number of clusters to filter pairs', call. = TRUE)
 }else if(length(Is_Zero) != 0){

   No_Zero_CV <- CheckZeroCV(MST = MST,
                             cluMem = cluMem,
                             Membership_Correction_Data = Membership_Correction_Data,
                             corGene = corGene,
                             Zero_Correction = Zero_Correction
                             )

   Membership_Correction_Data <- No_Zero_CV[["Membership_Correction_Data"]]
   corGene <- No_Zero_CV[["Correction_Matrix"]]

 }

 B2_Corrected <-  queBatch + (corGene  %*% t(corCell/rowSums(corCell)) )

 ### Set data lists to return
  Membership_Data <- list("Cluster Membership" = nMem, "Membership Correction Data" = Membership_Correction_Data )

  Correction_Data <- list("Correction Matrix" = corGene, "MNN Pairs" = Pairs,
                          "Membership Data" = Membership_Data, "Fuzzy Data" = Fuzzy_Data )

  Corrected_Batches <- list("Reference Batch (B1)" = refBatch, "Query Batch (B2)" = queBatch,
                            "Corrected Query Batch"= B2_Corrected, "Correction Data" = Correction_Data )

  if(Verbose)
    toc()

  return(Corrected_Batches)
}











