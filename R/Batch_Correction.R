
#' Correct_Batches
#'
#' Batch-Effect correction over a list of single cell batches
#'
#' @param Batches List of batches to integrate. Batches should contain the same number of genes as rows.
#' @param Query_Batch_Cell_Types A number indicating the expected number of cells types on the batches to integrate. The default value is set as a string "Surprise-me" on which an estimation of the cell types is defined.
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
#' z <- Correct_Batches(Batches)
#' Corrected <- z$`Batches Integrated`
#'
#' Uncorrected_PCA <- prcomp(t(cbind(Batches[[1]], Batches[[2]])))
#' plot(Uncorrected_PCA$x[,1:2])
#' Corrected_PCA <- prcomp(t(z$`Batches Integrated`))
#' plot(Corrected_PCA$x[,1:2])
#'
#' @return A list containing the integrated datasets as matrix and the correction data .
#' @export
#'
Correct_Batches <- function(Batches, Query_Batch_Cell_Types = "Surprise-me",
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
                            Estimation = "Average",
                            FilterPairs = FALSE,
                            perCellMNN = 0.08,
                            ...
                            ){

  if(Verbose)
    tic("\nTotal correction time ")

  #Init
  Names_Batches <- names(Batches)
  Num_Batches <- length(Batches)
  Corrected_Batches <- list()
  Batches_Integrated <- NULL
  Order <- NULL
  change <- FALSE
  Was_Integrated <- matrix(FALSE, nrow = Num_Batches, ncol = 1 )
  rownames(Was_Integrated) <- c( as.character(1:nrow(Was_Integrated)) )

  #Check input batches as matrices
  for (i in 1:Num_Batches) {
    Batches[[i]] <- as.matrix(Batches[[i]])
  }

  #In hierarchical mode, the pairs of the batches are checked in order to decide which batches are integrated first. The logic
  #is that more similar batches would share a higher number of pairs
  if(Hierarchical == TRUE & Num_Batches >2 ){

    i <- 1

    while (i < nrow(Was_Integrated)) {

      if(Was_Integrated[i] == FALSE){

        nCells_Bi <- ncol(Batches[[i]])
        N_Pairs_Bj <- as.vector(NULL)
        Last_Batch <- nrow(Was_Integrated)
        Num_Batches_2_Integrates <- length(which(Was_Integrated == FALSE))

        for(j in (i+1):Last_Batch){

          if( (Was_Integrated[j] == FALSE) &  (Num_Batches_2_Integrates > 2) ){

            nCells_Bj <- ncol(Batches[[j]])

            #Cosine normalization before getting pairs
            cnRefBatch <- batchelor::cosineNorm(Batches[[i]])
            cnQueBatch <- batchelor::cosineNorm(Batches[[j]])

            PCA_Batches <- prcomp_irlba(  t( cbind(cnRefBatch,cnQueBatch) ) )
            PCA_Bi <- PCA_Batches$x[1:nCells_Bi,]
            PCA_Bj <- PCA_Batches$x[(nCells_Bi+1):(nCells_Bi + nCells_Bj),]

            Pairs <- Get_MNN_Pairs(B1 = t(PCA_Bi),B2 = t(PCA_Bj),  k_Neighbors = 30)

            #N_Pairs_Bj <- rbind(N_Pairs_Bj, (nrow(Pairs$Pairs)/nCells_Bj))
            N_Pairs_Bj <- rbind(N_Pairs_Bj, nrow(Pairs$Pairs))
            rownames(N_Pairs_Bj)[nrow(N_Pairs_Bj)] <- j
          }

        }

        if(!is.null(N_Pairs_Bj)){
          Query <- as.integer( rownames(N_Pairs_Bj)[ which(N_Pairs_Bj == max(N_Pairs_Bj) )] )
        }else{
          Query <- Last_Batch
        }

        Ref <- i
        Names_Batches <- names(Batches)

        # TODO: ver si sirve, sino borrarlo
        if(ncol(Batches[[Ref]]) < ncol(Batches[[Query]]) ){
          temp <- Ref
          Ref <- Query
          Query <- temp
          rm(temp)
        }

        if(Verbose)
          cat(paste('\nINTEGRATING', Names_Batches[Query],"INTO", Names_Batches[Ref],"\n", sep = " ") )

        Correction <- Correct_Batch(refBatch = Batches[[Ref]],
                                    queBatch = Batches[[Query]],
                                    Query_Batch_Cell_Types = Query_Batch_Cell_Types,
                                    Sampling = Sampling,
                                    Number_Samples = Number_Samples,
                                    k_Neighbors = k_Neighbors,
                                    #PCA = PCA,
                                    Dimensions = Dimensions,
                                    Max_Membership = Max_Membership,
                                    Fuzzy = Fuzzy,
                                    Verbose = Verbose,
                                    #Cosine_Norm = Cosine_Norm,
                                    Estimation = Estimation,
                                    FilterPairs = FilterPairs,
                                    perCellMNN = perCellMNN
                                    )

        New_Name <- paste(Names_Batches[Ref],Names_Batches[Query],sep = "/")
        Corrected_Batches[[New_Name]] <- Correction
        Batches[[New_Name]] <- cbind( Batches[[Ref]], Correction[["Corrected Query Batch"]] )
        Names_Batches <- c(Names_Batches, New_Name)
        Was_Integrated[Ref] <- TRUE
        Was_Integrated[Query] <- TRUE
        Was_Integrated <- rbind(Was_Integrated, FALSE)

      }
      i <- i+1
    }

    Batches_Integrated <- Batches[[length(Batches)]]
    for (Batch in 1:Num_Batches) {
      Order <- c(Order, colnames(Batches[[Batch]]))
    }
    Batches_Integrated <- Batches_Integrated[,Order]

    Corrected_Batches[["Batches Integrated"]] <- Batches_Integrated

  }else{  #If the integration is not hierarchical

    for(i in 2:Num_Batches){

      Ref <- Batches[[1]]
      Query <- Batches[[i]]

      # TODO: ver si sirve, sino borrarlo
      if(ncol(Ref) < ncol(Query) ){
        temp <- Ref
        Ref <- Query
        Query <- temp
        rm(temp)
        change = TRUE
      }

      if(Verbose){
       if(change == FALSE){
         cat(paste('\nINTEGRATING', Names_Batches[i],"INTO", Names_Batches[1],"\n", sep = " ") )
       }else{
         cat(paste('\nINTEGRATING', Names_Batches[1] ,"INTO", Names_Batches[i] ,"\n", sep = " ") )
       }
      }


      Correction <- Correct_Batch(refBatch = Ref,
                                  queBatch = Query,
                                  Query_Batch_Cell_Types = Query_Batch_Cell_Types,
                                  Sampling = Sampling,
                                  Number_Samples = Number_Samples,
                                  k_Neighbors = k_Neighbors,
                                  #PCA = PCA,
                                  Dimensions = Dimensions,
                                  Max_Membership = Max_Membership,
                                  Fuzzy = Fuzzy,
                                  Verbose = Verbose,
                                  #Cosine_Norm = Cosine_Norm,
                                  Estimation = Estimation,
                                  FilterPairs = FilterPairs,
                                  perCellMNN = perCellMNN
                                  )

      New_Name <- paste(Names_Batches[1],Names_Batches[i],sep = "/")
      Corrected_Batches[[New_Name]] <- Correction
      names(Corrected_Batches[[New_Name]]) <- c(paste("Reference Batch (",Names_Batches[1] ,")", sep = ""),
                                                paste("Query Batch (",Names_Batches[i] , ")", sep = ""),
                                                "Corrected Query Batch", "Correction Data")
      Batches[[1]] <- cbind( Ref, Correction[["Corrected Query Batch"]] )
      Names_Batches[1] <- New_Name

      change = FALSE

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
#' @param Query_Batch_Cell_Types A number indicating the expected number of cells types on the batches to integrate. The default value is set as a string "Surprise-me" on which an estimation of the cell types is defined.
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
                          Query_Batch_Cell_Types = "Surprise-me",
                          Sampling = NULL,
                          Number_Samples = NULL,
                          Pairs = NULL,
                          Cells_Index_Query = NULL,
                          Cells_Index_Reference = NULL,
                          k_Neighbors = 30,
                          #PCA = TRUE, #TODO:implement for different input data
                          Dimensions = 50,
                          Max_Membership = 5,
                          Fuzzy = TRUE,
                          Verbose = FALSE,
                          #Cosine_Norm = TRUE, # TODO: implement for different input data
                          Estimation = "Average",
                          FilterPairs = FALSE,
                          perCellMNN = 0.08
                          ){

  if(Verbose)
    tic("\n Correction time")

  #Check input batches as matrices
  if(!is.matrix(refBatch)){
    refBatch <- as.matrix(refBatch)
  }
  if(!is.matrix(queBatch)){
    queBatch <- as.matrix(queBatch)
  }

  Membership_Pairs <- NULL
  Membership_Correction_Data <- list()
  Correction_Matrix <- NULL
  Fuzzy_Data <- NULL
  Num_Memberships <- NULL

  if( is.numeric(Query_Batch_Cell_Types) ){
    Num_Memberships <- Query_Batch_Cell_Types
  }else{
    if(Query_Batch_Cell_Types != "Surprise-me"){
      warning('\nWarning: Query_Batch_Cell_Types set value not recognized. Using "Surprise-me" instead', call. = TRUE)
    }
  }

  # TODO: correguir esto para no crear nuevos datasets, solo crear un indice. USAR INDICE SOLO EN ENCONTRAR MNN

  # if( !is.null(Cells_Index_Reference) ){
  #   if ( Num_Memberships > 1 ){
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
  #   if ( Num_Memberships > 1 ){
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

    PCA_Batches <- prcomp_irlba(t(cbind(batchelor::cosineNorm(refBatch),
                                        batchelor::cosineNorm(queBatch))),
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

 if(is.null(Num_Memberships)){

   if(Verbose)
    cat("\n\nFinding number of memberships")

   Num_Memberships <- pamk(pcaQue[,1:10],
                          krange = 1:Max_Membership,
                          usepam = (if(nCellsQue < 2000) TRUE else FALSE))$nc

   if(Verbose)
    cat(paste('\n\tNumber of memberships found:', Num_Memberships) )
 }

 #Cluster in memberships
 Cluster_Membership <- kmeans(pcaQue[,1:10],Num_Memberships)

 #INIT Correction Matrix
 Correction_Matrix <- matrix(0, nrow = nrow(refBatch), ncol = Num_Memberships)

 Zero_Correction <- rep(FALSE, Num_Memberships)

 for(Membership in 1:Num_Memberships){

   if(Verbose)
    cat(paste('\n\nAnalyzing Membership ', Membership) )

   #Membership cell index
   Membership_Cells_Index <- which(Cluster_Membership$cluster == Membership)
   #Membership cells number
   numCellMembership <- ncol(queBatch[,Membership_Cells_Index])

   #########################
   ###Pairs by membership###
   #########################
   Membership_Pairs_Index <- integer()
   for (j in 1:length(Membership_Cells_Index)) {
     Membership_Pairs_Index <- c(Membership_Pairs_Index,which(Pairs[,1]==Membership_Cells_Index[j]))
   }
   #Subset of pairs corresponding to the membership
   Membership_Pairs <- Pairs[Membership_Pairs_Index, ]

   #########################
   ###Pairs by clustering###
   #########################
   if(FilterPairs){

     if (nrow(Membership_Pairs) > 10){
        Pairs_Select <- Pairs_Selection(B1 = t(pcaRef),
                                        B2 = t(pcaQue),
                                        Pairs = Membership_Pairs,
                                        Verbose = Verbose)

        Selected_Pairs <- Pairs_Select[['Selected Pairs']]

        if(Verbose)
          cat( paste( '\n\tNumber of selected pairs:', nrow(Selected_Pairs) ) )
      }else {
        warning('\nWarning: Not enough pairs found for this Membership. No pairs selection is performed', call. = TRUE)
        Selected_Pairs <- Membership_Pairs
        Pairs_Select <- NULL
      }

   }else{
      Selected_Pairs <- Membership_Pairs
      Pairs_Select <- NULL
   }

   norNumPairs <- (ceiling(nrow(Selected_Pairs)/k_Neighbors))/(numCellMembership)

   if (norNumPairs > perCellMNN){

      ####################################
      ###Estimation of the batch effect###
      ####################################

     if(Estimation == "EKF"){
       if(Verbose)
         cat("\n\n\tEXTENDED KALMAN FILTER")

       Estimation_Data <- EKF_BE(B1 = refBatch,
                                 B2 = queBatch,
                                 Pairs = Selected_Pairs,
                                 Sampling = Sampling,
                                 Number_Samples = Number_Samples,
                                 Verbose = Verbose
                                 )
     }

     if(Estimation == "Average"){
       if(Verbose)
         cat("\n\n\tMedian Method")

       Estimation_Data <- Average_BE(B1 = refBatch,
                                     B2 = queBatch,
                                     Pairs = Selected_Pairs
                                     )
     }

     Correction_Vector <- Estimation_Data[["Correction Vector"]]
     Correction_Matrix[,Membership] <- Correction_Vector

    }else{
      warning('\nWarning: Not enough pairs found for this Membership. No correction is performed', call. = TRUE)
      Estimation_Data <- NULL
      Correction_Vector <- Correction_Matrix[,Membership]
      Zero_Correction[Membership] <- TRUE
    }

   Membership_Correction_Data[[paste("Membership", Membership)]] <- list( "Cells Index" = Membership_Cells_Index, "Pairs Selection Data" = Pairs_Select,
                                                                          "Sampled MNN Pairs" = Estimation_Data[["Sampled Pairs"]],
                                                                          "Correction Vector" = Correction_Vector   )

 }

 #################
 ###   FUZZY   ###
 #################

 ####INIT####
 Correction_Memberships <- matrix(0, nrow = nCellsQue, ncol = Num_Memberships )

 #Set column names according to number of memberships
 Name_Col <- NULL
 for (Mem in 1:Num_Memberships) {
   Name_Col <- rbind(Name_Col, paste("Mem",Mem,sep = "-") )
 }
 colnames(Correction_Memberships) <- Name_Col

 #Each cell is initialized according to its membership. Initilization is 1 to its membership and 0 to the other memberships
 for (Cell in 1:nCellsQue){
   Cell_Mem <- Cluster_Membership$cluster[Cell]
   Correction_Memberships[Cell,Cell_Mem] <- 1
 }

 #Fuzzy process and Correction
 if(Fuzzy == TRUE & Num_Memberships > 1 ){

   if(Verbose)
    cat('\n\nFUZZY ')

   Fuzzy_Data <- Fuzzy(Cluster_Membership = Cluster_Membership,
                       Cells_PCA = pcaQue,
                       Correction_Memberships = Correction_Memberships,
                       Verbose = Verbose
                       )

   Correction_Memberships <- Fuzzy_Data$`Fuzzy Memberships`
   MST <- Fuzzy_Data$MST

 }else{
   MST <- mst(dist(Cluster_Membership$centers[,1:2] ) )
 }

 #No Zero Correction Vectors
 Is_Zero <- which(Zero_Correction == TRUE)
 if(length(Is_Zero) == Num_Memberships){
   warning('\nWarning: No correction vectors where found. Consider using a higher number of kNN or a lower number of clusters to filter pairs', call. = TRUE)
 }else if( (length(Is_Zero) != 0) ){

   No_Zero_CV <- CheckZeroCV(MST = Fuzzy_Data$MST,
                             Cluster_Membership = Cluster_Membership,
                             Membership_Correction_Data = Membership_Correction_Data,
                             Correction_Matrix = Correction_Matrix,
                             Zero_Correction = Zero_Correction
                             )

   Membership_Correction_Data <- No_Zero_CV[["Membership_Correction_Data"]]
   Correction_Matrix <- No_Zero_CV[["Correction_Matrix"]]

 }

 B2_Corrected <-  queBatch + (Correction_Matrix  %*% t(Correction_Memberships/rowSums(Correction_Memberships)) )

 ### Set data lists to return
  Membership_Data <- list("Cluster Membership" = Cluster_Membership, "Membership Correction Data" = Membership_Correction_Data )

  Correction_Data <- list("Correction Matrix" = Correction_Matrix, "MNN Pairs" = Pairs,
                          "Membership Data" = Membership_Data, "Fuzzy Data" = Fuzzy_Data )

  Corrected_Batches <- list("Reference Batch (B1)" = refBatch, "Query Batch (B2)" = queBatch,
                            "Corrected Query Batch"= B2_Corrected, "Correction Data" = Correction_Data )

  if(Verbose)
    toc()

  return(Corrected_Batches)
}











