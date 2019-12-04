
#' Correct_Batches
#'
#' Batch-Effect correction over a list of single cell batches
#'
#' @param Batches List of batches to integrate. Batches should contain the same number of genes as rows.
#' @param Query_Batch_Cell_Types A number indicating the expected number of cells types on the batches to integrate. The default value is set as a string "Surprise-me" on which an estimation of the cell types is defined.
#' @param Similar_Cells A string value indicating in a semi-supervised the way MNNs pairs should be filtered. Accepted input values are "Low", "Medium" and "High".
#' @param Num_Clusters Number of clusters used to filter MNNs pairs.
#' @param Sampling Whether or not sampling of MNNs pairs is used on the estimation process.
#' @param Number_Samples Number of MNNs pairs samples used on the estimation process.
#' @param k_Neighbors Number of k-nearest-neighbors used to find MNNs pairs.
#' @param PCA Whether or not MNNs pairs are found under a principal components representation.
#' @param Dimensions PCA dimensions used to find MNNs pairs.
#' @param Max_Membership Maximum number of memberships used when memberships are automatically defined.
#' @param Fuzzy Whether or not a fuzzy logic join is used on the local correction vectors.
#' @param Hierarchical Whether or not a hiearchical integration scheme is used when correcting more than two batches.
#' @param Verbose Print output
#'
#' @details Correct_Batches is non-linear/linear hybrid method for single-cell batch-effect correction that couples identification of similar cells
#'  between datasets using Mutual Nearest Neighbors (MNNs) with an Extended Kalman Filter (EKF).
#'
#'  A non-linear correction is performed using fuzzy logic to join a set of linear correction vectors which are cell-type locally estimated.
#'
#' @examples
#' Batches <- SimBatches$Batches
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
                            Similar_Cells = "High", Num_Clusters = NULL, Sampling = NULL,
                            Number_Samples = NULL, k_Neighbors = 20, PCA = TRUE,
                            Dimensions = 30, Max_Membership = 5, Fuzzy = TRUE,
                            Hierarchical = TRUE, Verbose = FALSE ){

  if(Verbose)
    tic("\nTotal correction time ")

  #Init
  Names_Batches <- names(Batches)
  Num_Batches <- length(Batches)
  Corrected_Batches <- list()
  Batches_Integrated <- NULL
  Order <- NULL
  Was_Integrated <- matrix( rep(FALSE, Num_Batches ), ncol = 1 )
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

            PCA_Batches <- prcomp_irlba(  t( cbind(Batches[[i]], Batches[[j]]) ) )
            PCA_Bi <- PCA_Batches$x[1:nCells_Bi,]
            PCA_Bj <- PCA_Batches$x[(nCells_Bi+1):(nCells_Bi + nCells_Bj),]

            Pairs <- Get_MNN_Pairs(B1 = t(PCA_Bi),B2 = t(PCA_Bj),  k_Neighbors = 5)

            N_Pairs_Bj <- rbind(N_Pairs_Bj, (nrow(Pairs$Pairs)/nCells_Bj))
            #N_Pairs_Bj <- rbind(N_Pairs_Bj, nrow(Pairs$Pairs))
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

        if(Verbose)
          cat(paste('\nINTEGRATING', Names_Batches[Query],"INTO", Names_Batches[Ref],"\n", sep = " ") )

        Correction <- Correct_Batch(Reference_Batch = Batches[[Ref]], Query_Batch = Batches[[Query]],
                                    Query_Batch_Cell_Types = Query_Batch_Cell_Types,
                                    Similar_Cells = Similar_Cells, Num_Clusters = Num_Clusters,  Sampling = Sampling,
                                    Number_Samples = Number_Samples, k_Neighbors = k_Neighbors, PCA = PCA,
                                    Dimensions = Dimensions, Max_Membership = Max_Membership, Fuzzy = Fuzzy)


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

      if(Verbose)
        cat(paste('\nINTEGRATING', Names_Batches[i],"INTO", Names_Batches[1],"\n", sep = " ") )

      Correction <- Correct_Batch(Reference_Batch = Ref, Query_Batch = Query, Query_Batch_Cell_Types = Query_Batch_Cell_Types,
                                  Similar_Cells = Similar_Cells, Num_Clusters = Num_Clusters,  Sampling = Sampling,
                                  Number_Samples = Number_Samples, k_Neighbors = k_Neighbors, PCA = PCA,
                                  Dimensions = Dimensions, Max_Membership = Max_Membership, Fuzzy = Fuzzy, Verbose = Verbose)

      New_Name <- paste(Names_Batches[1],Names_Batches[i],sep = "/")
      Corrected_Batches[[New_Name]] <- Correction
      names(Corrected_Batches[[New_Name]]) <- c(paste("Reference Batch (",Names_Batches[1] ,")", sep = ""),
                                                paste("Query Batch (",Names_Batches[i] , ")", sep = ""),
                                                "Corrected Query Batch", "Correction Data")
      Batches[[1]] <- cbind( Ref, Correction[["Corrected Query Batch"]] )
      Names_Batches[1] <- New_Name

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
#' @param Reference_Batch Batch to use as reference for the integration.
#' @param Query_Batch Batch to correct.
#' @param Query_Batch_Cell_Types A number indicating the expected number of cells types on the batches to integrate. The default value is set as a string "Surprise-me" on which an estimation of the cell types is defined.
#' @param Similar_Cells A string value indicating in a semi-supervised the way MNNs pairs should be filtered. Accepted input values are "Low", "Medium" and "High".
#' @param Num_Clusters Number of clusters used to filter MNNs pairs.
#' @param Sampling Whether or not samples MNNs pairs' samples are used on the estimation process.
#' @param Number_Samples Number of MNNs pairs' samples used on the estimation process.
#' @param Pairs A matrix containing MNNs pairs. First column corresponds to query-batch cell indexes.
#' @param Cells_Index_Query Index of cells from the query-batch used as observations of the batch-effect.
#' @param Cells_Index_Reference Index of cells from the reference-batch used as observations of the batch-effect.
#' @param k_Neighbors Number of k-nearest-neighbors used to find MNNs pairs.
#' @param PCA Whether or not MNNs pairs are found under a principal components representation.
#' @param Dimensions PCA dimensions used to find MNNs pairs.
#' @param Max_Membership Maximum number of memberships used when memberships are automatically defined.
#' @param Fuzzy Whether or not a fuzzy logic join is used on the local correction vectors.
#' @param Verbose Print output
#'
#' @details Canek, a new non-linear/linear hybrid method for batch-effect correction that couples identification of similar cells
#'  between datasets using Mutual Nearest Neighbors (MNNs) with an Extended Kalman Filter (EKF).
#'
#'  A non-linear correction is performed by using fuzzy logic to join a set of linear correction vectors which are cell-type locally estimated.
#'
#' @examples
#' x <- SimBatches$Batches[[1]]
#' y <- SimBatches$Batches[[2]]
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
Correct_Batch <- function(Reference_Batch, Query_Batch, Query_Batch_Cell_Types = "Surprise-me",
                          Similar_Cells = "High", Num_Clusters = NULL, Sampling = NULL,
                          Number_Samples = NULL, Pairs = NULL, Cells_Index_Query = NULL,
                          Cells_Index_Reference = NULL, k_Neighbors = 20, PCA = TRUE,
                          Dimensions = 30,   Max_Membership = 5, Fuzzy = TRUE, Verbose = FALSE){

  if(Verbose)
    tic("\n Correction time")

  #Check input batches as matrices
  if(!is.matrix(Reference_Batch)){
    Reference_Batch <- as.matrix(Reference_Batch)
  }
  if(!is.matrix(Query_Batch)){
    Query_Batch <- as.matrix(Query_Batch)
  }

  #Initialization of variables
  B1 <- Reference_Batch
  B2 <- Query_Batch

  B2_Corrected <- B2
  Membership_Pairs <- NULL
  Membership_Correction_Data <- list()
  Correction_Matrix <- NULL
  Fuzzy_Data <- NULL
  Num_Memberships = NULL

  if( is.numeric(Query_Batch_Cell_Types) ){
    Num_Memberships <- Query_Batch_Cell_Types
  }else{
    if(Query_Batch_Cell_Types != "Surprise-me"){
      warning('\nWarning: Query_Batch_Cell_Types set value not recognized. Using "Surprise-me" instead', call. = TRUE)
    }
  }

  if( !is.null(Cells_Index_Reference) ){
    if ( Num_Memberships > 1 ){
      warning('\nWarning: CELLS INDEX NOT USED. Cannot use cells index for more than one membership function', call. = TRUE)
      B1_Selected <- B1
    } else{
      B1_Selected <- B1[,Cells_Index_Reference]
    }
  } else{
    B1_Selected <- B1
  }

  if( !is.null(Cells_Index_Query) ){
    if ( Num_Memberships > 1 ){
      warning('\nWarning: CELLS INDEX NOT USED. Cannot use cells index for more than one membership function', call. = TRUE)
      B2_Selected <- B2
    }else{
      B2_Selected <- B2[,Cells_Index_Query]
    }
  }else{
    B2_Selected <- B2
  }

  B1_Selected_Num_Cells <- ncol(B1_Selected)
  B2_Selected_Num_Cells <- ncol(B2_Selected)
  Num_Cells <-B1_Selected_Num_Cells + B2_Selected_Num_Cells
  Num_genes <- nrow(B1_Selected)

  #Setting the Number of Clusters
  if(is.null(Num_Clusters)){
    if(Similar_Cells == "Low"){
      if(B2_Selected_Num_Cells > 100){
        Num_Clusters <- 100
      }else{
        Num_Clusters <- 10
      }
      if(is.null(Sampling))
        Sampling <- FALSE
    }else if(Similar_Cells == "Medium"){
      Num_Clusters <- 3
      if(is.null(Sampling))
        Sampling <- FALSE
    }else{
      if(Similar_Cells != "High")
          warning('\nWarning: Similar_Cells set value not recognized. Using "High" instead', call. = TRUE)
      Num_Clusters <- 1
      if(is.null(Sampling))
        Sampling <- TRUE
    }
  }

  if( !is.null(Pairs) ){
    Pairs <- Pairs
  }else{

    if(PCA == TRUE){

      PCA_Batches <- prcomp_irlba( rbind(t(B1_Selected), t(B2_Selected) ), n = Dimensions)
      PCA_B1 <- PCA_Batches$x[1:B1_Selected_Num_Cells,]
      PCA_B2 <- PCA_Batches$x[(B1_Selected_Num_Cells+1):Num_Cells,]

      if(Verbose)
        cat( paste("\n\nFinding mutual nearest neighbors from ", k_Neighbors,"nearest neighbors") )

      Pairs <- Get_MNN_Pairs(B1 = t(PCA_B1),B2 = t(PCA_B2),  k_Neighbors = k_Neighbors)

    }else{
      Pairs <- Get_MNN_Pairs(B1 = B1_Selected, B2 = B2_Selected , k_Neighbors = k_Neighbors)
    }

    Pairs = Pairs$Pairs
  }

 if(Verbose)
  cat(paste( '\n\tNumber of MNN pairs found:', nrow(Pairs) ))

 if( is.null(Num_Memberships) ){

   if(Verbose)
    cat("\n\nFinding number of memberships")

   if( !exists("PCA_B2") ){
     PCA_B2 <- prcomp_irlba( t(B2_Selected) )
     PCA_B2 <- PCA_B2$x
   }

   if(ncol(B2_Selected) < 2000){
     usepam <- TRUE
   }else{
     usepam <- FALSE
   }

   Num_Memberships <- pamk( PCA_B2[,1:3], krange = 2:Max_Membership, usepam = usepam )
   Num_Memberships <- Num_Memberships$nc

   if(Verbose)
    cat(paste('\n\tNumber of memberships found:', Num_Memberships) )
 }

 #Cluster in memberships
 Cluster_Membership <- kmeans(PCA_B2[,1:2],Num_Memberships)

 #INIT Correction Matrix
 Correction_Matrix <- matrix(rep(0,Num_genes*Num_Memberships), ncol = Num_Memberships)

 for(Membership in 1:Num_Memberships){

   if(Verbose)
    cat(paste('\n\nAnalyzing Membership ', Membership) )

   #Membership cell index
   Membership_Cells_Index <- which(Cluster_Membership$cluster == Membership)
   #Membership cells subset
   B2_Membership <- B2_Selected[,Membership_Cells_Index]

   #########################
   ###Pairs by membership###
   #########################
   Membership_Pairs_Index <- which(Pairs[,1]==Membership_Cells_Index[1])
   for (j in 2:length(Membership_Cells_Index) ) {
     Membership_Pairs_Index <- c(Membership_Pairs_Index,which(Pairs[,1]==Membership_Cells_Index[j]))
   }
   #Subset of pairs corresponding to the membership
   Membership_Pairs <- Pairs[Membership_Pairs_Index,]

   #########################
   ###Pairs by clustering###
   #########################
   if( Num_Clusters != 1 ){

     if (length(Membership_Pairs)>20){
        Pairs_Select <- Pairs_Selection(B1 = t(PCA_B1), B2 = t(PCA_B2), Pairs = Membership_Pairs, Num_Clusters = Num_Clusters, Verbose = Verbose)
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

   if ( length(Selected_Pairs) > 20 ){

      #########################################################
      ###Extended Kalman Filter to estimate the batch effect###
      #########################################################
      if(Verbose)
        cat("\n\n\tEXTENDED KALMAN FILTER")

      Estimation_Data <- EKF_BE(B1 = B1_Selected,B2 = B2_Selected, Pairs = Selected_Pairs, Sampling = Sampling,
                                  Number_Samples = Number_Samples, Verbose = Verbose)

      Correction_Vector <- Estimation_Data[["Correction Vector"]]

      Correction_Matrix[,Membership] <- Correction_Vector

    }else{
      warning('\nWarning: Not enough pairs found for this Membership. No correction is performed', call. = TRUE)
      Estimation_Data <- NULL
      Correction_Vector <- Correction_Matrix[,Membership]
    }

   Membership_Correction_Data[[paste("Membership", Membership)]] <- list( "Cells Index" = Membership_Cells_Index, "Pairs Selection Data" = Pairs_Select,
                                                                          "Sampled MNN Pairs" = Estimation_Data[["Sampled Pairs"]],
                                                                          "Correction Vector" = Correction_Vector   )

 }

 #################
 ###   FUZZY   ###
 #################

 ####INIT####
 Correction_Memberships <- matrix( rep(0, B2_Selected_Num_Cells*Num_Memberships), ncol = Num_Memberships )

 #Set column names according to number of memberships
 Name_Col <- NULL
 for (Mem in 1:Num_Memberships) {
   Name_Col <- rbind(Name_Col, paste("Mem",Mem,sep = "-") )
 }
 colnames(Correction_Memberships) <- Name_Col

 #Each cell is initialized according to its membership. Initilization is 1 to its membership and 0 to the other memberships
 for (Cell in 1:B2_Selected_Num_Cells){
   Cell_Mem <- Cluster_Membership$cluster[Cell]
   Correction_Memberships[Cell,Cell_Mem] <- 1
 }

 #Fuzzy process and Correction
 if(Fuzzy == TRUE & Num_Memberships > 1 ){

   if(Verbose)
    cat('\n\nFUZZY ')

   Fuzzy_Data <- Fuzzy(Cluster_Membership = Cluster_Membership, Cells_PCA = PCA_B2, Correction_Memberships = Correction_Memberships,
                       Verbose = Verbose)
   B2_Corrected <-  B2 + (Correction_Matrix  %*% t(Fuzzy_Data$`Fuzzy Memberships`) )

 }else{

   B2_Corrected <-  B2 + (Correction_Matrix  %*% t(Correction_Memberships) )

 }

 ### Set data lists to return
  Membership_Data <- list("Cluster Membership" = Cluster_Membership, "Membership Correction Data" = Membership_Correction_Data )

  Correction_Data <- list("Correction Matrix" = Correction_Matrix, "MNN Pairs" = Pairs,
                          "Membership Data" = Membership_Data, "Fuzzy Data" = Fuzzy_Data )

  Corrected_Batches <- list("Reference Batch (B1)" = B1, "Query Batch (B2)" = B2,
                            "Corrected Query Batch"= B2_Corrected, "Correction Data" = Correction_Data )

  if(Verbose)
    toc()

  return(Corrected_Batches)
}











