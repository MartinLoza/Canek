
#' Correct_Batches
#'
#' Function to correct batch effect over two batches
#'
#' @param Batches List of batches to integrate. Batches should contain the same number of genes as rows.
#' @param Query_Batch_Cell_Types A number indicating the expected number of cells types on the batches to integrate. The default value is set as a string "Surprise-me" on which an estimation of the cell types is defined.
#' @param Similar_Cells A string value indicating in a semi-supervised the way MNNs pairs should be filtered. Accepted input values are "Low", "Medium" and "High".
#' @param Num_Clusters Number of clusters used to filter MNNs pairs.
#' @param Sampling Whether or not sampling of MNNs pairs is used on the estimation process.
#' @param Number_Samples A number defining the number of MNNs pairs samples to use on the estimation process.
#' @param k_Neighbors Number of k-nearest-neighbors used to find MNNs pairs.
#' @param PCA Whether or not MNNs pairs are found under a principal components representation.
#' @param Dimensions PCA dimensions used to find MNNs pairs.
#' @param Max_Membership Maximum number of memberships used when memberships are automatically defined.
#' @param Fuzzy Whether or not a fuzzy logic join is used on the local correction vectors.
#' @param Hierarchical Whether or not a hiearchical integration scheme is used when correcting more than two batches.
#'
#' @return A list containing a matrix with the integrated datasets.
#' @export
#'
#' @examples
Correct_Batches <- function(Batches, Query_Batch_Cell_Types = "Surprise-me",
                            Similar_Cells = "High", Num_Clusters = NULL, Sampling = NULL,
                            Number_Samples = NULL, k_Neighbors = 20, PCA = TRUE,
                            Dimensions = 50, Max_Membership = 5, Fuzzy = TRUE,
                            Hierarchical = TRUE ){

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

      cat(paste('\nINTEGRATING', Names_Batches[i],"INTO", Names_Batches[1],"\n", sep = " ") )

      Correction <- Correct_Batch(Reference_Batch = Ref, Query_Batch = Query, Query_Batch_Cell_Types = Query_Batch_Cell_Types,
                                  Similar_Cells = Similar_Cells, Num_Clusters = Num_Clusters,  Sampling = Sampling,
                                  Number_Samples = Number_Samples, k_Neighbors = k_Neighbors, PCA = PCA,
                                  Dimensions = Dimensions, Max_Membership = Max_Membership, Fuzzy = Fuzzy)

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

  toc()
  return(Corrected_Batches)
}

##Correct_Batch##
#Function to correct batch effect over two batches
# INPUT : Reference_Batch -> Reference batch
#         Query_Batch -> Query batch (This batch will be modified)
#         Pairs -> Cell pairs for batch effect estimations
#         Num_Clusters -> Number of clusters to use when selecting pairs
#         Sampling -> A boolean value indicating if the estimator will use sampling from the cell pairs
#         Number_Samples -> Number of samples to use from the cell pairs
# OUTPUT : Correction Vector containing the batch effect estimation. The vector is equal size as the number of genes.
Correct_Batch <- function(Reference_Batch, Query_Batch, Query_Batch_Cell_Types = "Surprise-me",
                          Similar_Cells = "High", Num_Clusters = NULL, Sampling = NULL,
                          Number_Samples = NULL, Pairs = NULL, Cells_Index_Query = NULL,
                          Cells_Index_Reference = NULL, k_Neighbors = 20, PCA = TRUE,
                          Dimensions = 50,   Max_Membership = 5, Fuzzy = TRUE){

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
    }else{
      B1_Selected <- B1[,Cells_Index_Reference]
    }
  }else{
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
    if(Similar_Cells == "High"){
      Num_Clusters <- 1
      if(is.null(Sampling))
        Sampling <- TRUE
    }else if(Similar_Cells == "Very Low"){
      #Num_Clusters <- floor(B2_Selected_Num_Cells/200)
      Num_Clusters <- 100
      if(is.null(Sampling))
        Sampling <- FALSE
    }else{
      if(Similar_Cells != "Low")
        warning('\nWarning: Similar_Cells set value not recognized. Using "Low instead', call. = TRUE)
      Num_Clusters <- 3
      if(is.null(Sampling))
        Sampling <- FALSE
    }
  }

  if( !is.null(Pairs) ){
    Pairs <- Pairs
  }else{

    if(PCA == TRUE){

      PCA_Batches <- prcomp_irlba( rbind(t(B1_Selected), t(B2_Selected) ), n = Dimensions)
      PCA_B1 <- PCA_Batches$x[1:B1_Selected_Num_Cells,]
      PCA_B2 <- PCA_Batches$x[(B1_Selected_Num_Cells+1):Num_Cells,]

      cat( paste("\n\nFinding mutual nearest neighbors from ", k_Neighbors,"nearest neighbors") )
      Pairs <- Get_MNN_Pairs(B1 = t(PCA_B1),B2 = t(PCA_B2),  k_Neighbors = k_Neighbors)

    }else{
      Pairs <- Get_MNN_Pairs(B1 = B1_Selected, B2 = B2_Selected , k_Neighbors = k_Neighbors)
    }

    Pairs = Pairs$Pairs
  }

 cat(paste( '\n\tNumber of MNN pairs found:', nrow(Pairs) ))

 if( is.null(Num_Memberships) ){

   cat("\n\nFinding number of memberships")

   if( !exists("PCA_B2") ){
     PCA_B2 <- prcomp_irlba( t(B2_Selected) )
     PCA_B2 <- PCA_B2$x
   }

   Num_Memberships <- pamk( PCA_B2[,1:3], krange = 2:Max_Membership )
   Num_Memberships <- Num_Memberships$nc
   cat(paste('\n\tNumber of memberships found:', Num_Memberships) )
 }

 #Cluster in memberships
 Cluster_Membership <- kmeans(PCA_B2[,1:2],Num_Memberships)

 #INIT Correction Matrix
 Correction_Matrix <- matrix(rep(0,Num_genes*Num_Memberships), ncol = Num_Memberships)

 for(Membership in 1:Num_Memberships){

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
        Pairs_Select <- Pairs_Selection(B1 = B1_Selected, B2 = B2_Membership, Pairs = Membership_Pairs, Num_Clusters = Num_Clusters)
        Selected_Pairs <- Pairs_Select[['Selected Pairs']]
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
      cat("\n\n\tEXTENDED KALMAN FILTER")
      Estimation_Data <- EKF_BE(B1 = B1_Selected,B2 = B2_Selected, Pairs = Selected_Pairs, Sampling = Sampling,
                                  Number_Samples = Number_Samples)

      Correction_Vector <- Estimation_Data[["Correction Vector"]]

      Correction_Matrix[,Membership] <- Correction_Vector

    }else{
      warning('\nWarning: Not enough pairs found for this Membership. No correction is performed', call. = TRUE)
      Estimation_Data <- NULL
      Correction_Vector <- NULL
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

   cat('\n\nFUZZY ')
   Fuzzy_Data <- Fuzzy(Cluster_Membership = Cluster_Membership, Cells_PCA = PCA_B2, Correction_Memberships = Correction_Memberships )
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

  toc()

  return(Corrected_Batches)
}

##Fuzzy##
# This function stablishes the fuzzification for the memberships. A minimum spanning tree (MST) is created among memberships,
# and the fuzzification is performed for each of the edges of the MST.
# INPUT : Fuzzy data list
# OUTPUT : Fuzzy data list
Fuzzy <- function(Cluster_Membership = NULL, Cells_PCA = NULL, Correction_Memberships = NULL){

  #INIT
  Num_Cells <- nrow(Cells_PCA)
  Num_Memberships <- nrow(Cluster_Membership$centers)
  PCA_Max <- 2
  Fuzzy <- rep(FALSE, Num_Cells)
  Fuzzy_Memberships <- Correction_Memberships
  Edges_Data <- list()


  # Create Minimum spanning tree (MST) by using centers of Memberships as nodes
  cat( '\n\tObtaining Minimum Spanning Tree' )
  Mst <- mst(dist( Cluster_Membership$centers[,1:PCA_Max] ) )  #TODO: If possible change to higher dimensions

  #Get edges from MST
  Edges <- Get_Edges(Mst)

  #Fuzzy process
  cat('\n\tFuzzificating cells from each edge')
  for(Edge in 1:nrow(Edges)){


    #Init nodes
    IN_Membership <- list()
    OUT_Membership <- list()
    IN_Node <- Edges[ Edge, 1 ]
    OUT_Node <- Edges[ Edge, 2 ]
    IN_Node_PCA <- Cluster_Membership$centers[ IN_Node, 1:PCA_Max ]
    OUT_Node_PCA <- Cluster_Membership$centers[ OUT_Node, 1:PCA_Max ]

    #Translate OUT node according to IN node PCA coordinates
    OUT_Node_PCA_Transformed <- OUT_Node_PCA - IN_Node_PCA

    #Find angle between the current nodes (we set it negative because we would correct this angle in further steps)
    Alpha <- (-Get_Rotation_Angle ( OUT_Node_PCA_Transformed ) )
    names(Alpha)<-"Alpha"

    #Get cells from both IN and OUT memberships
    IN_Membership_Cells_Index <- which(Cluster_Membership$cluster == IN_Node)
    IN_Membership_Cells <- Cells_PCA[IN_Membership_Cells_Index, 1:PCA_Max ]
    rownames(IN_Membership_Cells) <- IN_Membership_Cells_Index

    OUT_Membership_Cells_Index <- which(Cluster_Membership$cluster == OUT_Node)
    OUT_Membership_Cells <- Cells_PCA[which(Cluster_Membership$cluster == OUT_Node), 1:PCA_Max ]
    rownames(OUT_Membership_Cells) <- OUT_Membership_Cells_Index

    #Translate according cells from both memberships according to IN node PCA coordinates
    IN_Membership_Cells_Transformed <- sweep(IN_Membership_Cells,2,IN_Node_PCA,"-")
    OUT_Membership_Cells_Transformed <- sweep(OUT_Membership_Cells,2,IN_Node_PCA,"-")

    #Paso 6 Rotate the cells and OUT Node
    IN_Membership_Cells_Transformed <- Rotation(IN_Membership_Cells_Transformed, angle = Alpha)
    OUT_Membership_Cells_Transformed <- Rotation(OUT_Membership_Cells_Transformed, angle = Alpha)
    OUT_Node_PCA_Transformed <- c( (OUT_Node_PCA_Transformed['PC1']*cos(Alpha) - OUT_Node_PCA_Transformed['PC2']*sin(Alpha) ),
                                   (OUT_Node_PCA_Transformed['PC1']*sin(Alpha) + OUT_Node_PCA_Transformed['PC2']*cos(Alpha) ) )

    #Filter cells. Only cells that are located between the two nodes PCA coordinates are transformed.
    IN_Membership_Cells_Transformed_Selected_Index <- which(IN_Membership_Cells_Transformed[,1] >= 0)
    OUT_Membership_Cells_Transformed_Selected_Index <- which(OUT_Membership_Cells_Transformed[,1] < OUT_Node_PCA_Transformed[1])

    IN_Membership_Cells_Filtered <- IN_Membership_Cells_Transformed[IN_Membership_Cells_Transformed_Selected_Index, 1:PCA_Max, drop = FALSE]
    OUT_Membership_Cells_Filtered <- OUT_Membership_Cells_Transformed[OUT_Membership_Cells_Transformed_Selected_Index, 1:PCA_Max, drop = FALSE]

    Cells_Filtered <- rbind( IN_Membership_Cells_Filtered,OUT_Membership_Cells_Filtered )

    #Get slope for fuzzification
    Slope <- (-1/OUT_Node_PCA_Transformed[1])
    Fuzzification <- NULL

    #Fuzzification
    Cells_Filtered_RowNames <- as.numeric(rownames(Cells_Filtered))

    for(Cell in 1:nrow(Cells_Filtered)){

      IN_Fuzzification <- 1 + (Slope*Cells_Filtered[Cell,1])
      OUT_Fuzzification <- ( 1 - IN_Fuzzification )
      Fuzzification <- rbind( Fuzzification, c(IN_Fuzzification, OUT_Fuzzification) )

      #If this cells has not been fuzzified before set fuzzification values
      if(Fuzzy[Cells_Filtered_RowNames[Cell]] == FALSE){

        Fuzzy_Memberships[Cells_Filtered_RowNames[Cell], IN_Node] <- IN_Fuzzification
        Fuzzy_Memberships[Cells_Filtered_RowNames[Cell], OUT_Node] <- OUT_Fuzzification

      }else{    ##If this cells was already fuzzified, set the fuzzy value as the average of fuzzy values

        Fuzzy_Memberships[Cells_Filtered_RowNames[Cell], IN_Node] <-
          mean( c(Fuzzy_Memberships[Cells_Filtered_RowNames[Cell], IN_Node], IN_Fuzzification) )
        Fuzzy_Memberships[Cells_Filtered_RowNames[Cell], OUT_Node] <-
          mean( c(Fuzzy_Memberships[Cells_Filtered_RowNames[Cell], OUT_Node],OUT_Fuzzification) )

      }

    }

    rownames(Fuzzification) <- Cells_Filtered_RowNames
    colnames(Fuzzification) <- c("IN-Mem-Fuzzy", "OUT-Mem-Fuzzy")

    IN_Membership <- list("Cells" = IN_Membership_Cells, "Transformed" = IN_Membership_Cells_Transformed,
                          "Filtered" = IN_Membership_Cells_Filtered)
    OUT_Membership <- list("Cells" = OUT_Membership_Cells, "Transformed" = OUT_Membership_Cells_Transformed,
                          "Filtered" = OUT_Membership_Cells_Filtered)

    Edges_Data[[paste("Edge-", Edge, sep = "")]] <- list( "IN Node" = IN_Node, "OUT Node" = OUT_Node,
                                                          "Angle" = Alpha, "IN-Mem-Cells Data" = IN_Membership,
                                                          "OUT-Mem-Cells Data" = OUT_Membership, "Slope" = Slope,
                                                          "Fuzzification" = Fuzzification)

  }

  Fuzzy_Data <- list( "Fuzzy Memberships" = Fuzzy_Memberships, "MST" = Mst, "Fuzzy" =   Fuzzy, "Edges Data" = Edges_Data)

  return(Fuzzy_Data)
}

##Get_Rotation_Angle##
#Get correction angle according to quadrant
# INPUT : Node PCA coordinates
#
# OUTPUT : Angle in radians
Get_Rotation_Angle <- function( PCA_Coordinates = NULL ){

  Rotation_Angle <- atan(PCA_Coordinates['PC2']/PCA_Coordinates['PC1'])

  #For second and third quadrant we add pi to rotation angle
  if( (PCA_Coordinates['PC1'] < 0 & PCA_Coordinates['PC2'] >= 0) |
      (PCA_Coordinates['PC1'] < 0 & PCA_Coordinates['PC2'] < 0) ){   #Second Quadrant
    Rotation_Angle <- Rotation_Angle + pi
  }

  if( (PCA_Coordinates['PC1'] >= 0 & PCA_Coordinates['PC2'] < 0) )
    Rotation_Angle <- Rotation_Angle + (2*pi)

  return(Rotation_Angle)

}


##Pairs_Selection##
#Function to select pairs to be used on batch effect correction. The pairs are selected by clustering the query batch
#and analyzing the mean distance of pair cells on each cluster. The pairs of the cluster with minimum mean distance are returned.
# INPUT : B1 -> Batch 1 (Reference)
#         B2 -> Batch 2 (Query)
#         Pairs -> MNN Pairs
#         Num_Clusters -> Number of clusters to use
# OUTPUT : Matrix containing the selected pairs.
Pairs_Selection <- function(B1,B2, Pairs, Num_Clusters = 1){

  cat("\n\nSelecting Pairs by Clusters")
  cat(paste('\n\tClusters for selecting pairs:',Num_Clusters))
  #Initialization of varibles for clustering and selection of cell pairs
  Mu <- matrix(0L, nrow = Num_Clusters, ncol = 1)
  SD <- matrix(0L, nrow = Num_Clusters, ncol = 1)

  # Clustering to select pairs
  Cluster <- kmeans(t(B2),Num_Clusters)

  #Initialization of variables for the correction
  Pairs_Dim <- dim(Pairs)
  Num_Pairs <- Pairs_Dim[1]
  Pairs <- matrix(c(as.vector(Pairs[1:Num_Pairs,1]), as.vector(Pairs[1:Num_Pairs,2])), nrow = Num_Pairs)
  #Samples <- Pairs

  for (i in 1:Num_Clusters) {

    #Cluster cells index
    Cl_Cells_Index<- which(Cluster$cluster ==i)

    #Cluster cells subset
    Cl_Cells<- B2[,Cl_Cells_Index]

    #Cluster cells pairs
    Cl_Pairs_Index <- which(Pairs[,1]==Cl_Cells_Index[1])
    for (j in 2:length(Cl_Cells_Index) ) {
      Cl_Pairs_Index <- c(Cl_Pairs_Index,which(Pairs[,1]==Cl_Cells_Index[j]))
    }

    #Subset of the pairs corresponding to the cluster
    Cl_Pairs <- Pairs[Cl_Pairs_Index,]

    #Number of pairs for this cluster
    Pairs_Dim <- dim(Cl_Pairs)
    Num_Pairs <- Pairs_Dim[1]
    Num_Pairs

    if (Num_Pairs == 0 || is.null(Num_Pairs) ){
      Mu[i] <- NA
      SD[i] <- NA

    } else {

      #Distance analysis of the pairs per cluster
      Distances <- matrix(0L, nrow = Num_Pairs, ncol = 1)

      for (j in 1:Num_Pairs) {
        Distances[j] <- dist(rbind(B1[,Cl_Pairs[j,2]],B2[,Cl_Pairs[j,1]]))
        Mu[i] <- mean(Distances)
        SD[i] <- sd(Distances)
      }
    }
  }

  #Select the cluster with minimum distance mean
  Min_Cl <- which.min(Mu)

  #Cluster cells index
  Cl_Cells_Index<- which(Cluster$cluster ==Min_Cl)

  #Cluster cells subset
  Cl_Cells<- B2[,Cl_Cells_Index]

  #Cluster cells pairs
  Cl_Pairs_Index <- which(Pairs[,1]==Cl_Cells_Index[1])
  for (j in 2:length(Cl_Cells_Index) ) {
    Cl_Pairs_Index <- c(Cl_Pairs_Index,which(Pairs[,1]==Cl_Cells_Index[j]))
  }

  #Subset of the pairs of the selected cluster
  Cl_Pairs <- Pairs[Cl_Pairs_Index,]

  return(list( "Selected Pairs" =Cl_Pairs, "Clusters" = Cluster, "Selected Cluster" = Min_Cl ) )
}

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

##Simulate_Batches##
#Function to simulate two batches containing an orthogonal batch effect.
# INPUT :
#
# OUTPUT :
Simulate_Batches <- function (B1_Ncells, B2_Ncells, Ngenes = 10, B1_mus, B2_mus, B1_CellsProb,
                              B2_CellsProb, B1_sds, B2_sds, B1ref_cols = c("red","yellow",'blue'),
                              B2ref_cols = c("orange", "black", "green"), Noise = FALSE )  {

  library(scran)
  require(Rtsne)
  require(ggplot2)


  ###############################################
  #   INIT                  #####################
  ###############################################
  B1ref_cols <- c("#000099", "#CC0000","#FFFF00")
  B2ref_cols <- c("#000099", "#CC0000","#00FF00")

  B1_shape <- rep(16, B1_Ncells)
  B2_shape <- rep(2, B2_Ncells)

  ###############################################
  #   SIMULATE BATCH 1      #####################
  ###############################################

  comp1 <- sample(1:3, prob=B1_CellsProb, size=B1_Ncells, replace=TRUE)

  # Sampling locations for cells in each component.
  set.seed(0)
  samples1 <- cbind(rnorm(n=B1_Ncells, mean=B1_mus[1,comp1], sd=B1_sds[1,comp1]),
                    rnorm(n=B1_Ncells, mean=B1_mus[2,comp1], sd=B1_sds[2,comp1]))

  # Set the colours for batch 1
  B1ref_cols <- B1ref_cols[comp1]

  # Random projection to D dimensional space, to mimic high-dimensional expression data.
  set.seed(0)
  proj <- matrix(rnorm(Ngenes*B1_Ncells), nrow=Ngenes, ncol=2)
  A1 <- samples1 %*% t(proj)

  # Add normally distributed noise.
  if (Noise == TRUE){
    A1 <- A1 + rnorm(Ngenes*B1_Ncells)
  }
  #Batch tags
  rownames(A1) <- paste0("Cell", seq_len(B1_Ncells), "-1")
  colnames(A1) <- paste0("Gene", seq_len(Ngenes))

  B1 <- t(A1)

  ###############################################
  #   SIMULATE BATCH 2      #####################
  ###############################################

  comp2 <- sample(1:3, prob=B2_CellsProb, size=B2_Ncells, replace=TRUE)

  # Sampling locations for cells in each component.
  set.seed(0)
  samples2 <- cbind(rnorm(n=B2_Ncells, mean=B2_mus[1,comp2],sd=B2_sds[1,comp2]),
                    rnorm(n=B2_Ncells, mean=B2_mus[2,comp2],sd=B2_sds[2,comp2]))

  # Set the colours for batch 2
  B2ref_cols <- B2ref_cols[comp2]

  # Random projection to D dimensional space, to mimic high-dimensional expression data.
  set.seed(0)
  proj <- matrix(rnorm(Ngenes*B2_Ncells), nrow=Ngenes, ncol=2)
  A2 <- samples2 %*% t(proj)

  # Add batch effect and normally distributed noise.
  Batch_Effect_Vector <- rnorm(Ngenes)
  Batch_Effect <- matrix(rep(Batch_Effect_Vector, each=B2_Ncells), ncol=Ngenes);
  A2 <- A2 + Batch_Effect # gene-specific batch effect (genes are columns)

  # Add normally distributed noise.
  if (Noise == TRUE){
    A2 <- A2 + rnorm(Ngenes*B2_Ncells) # noise
  }

  #Batch tags
  rownames(A2) <- paste0("Cell", seq_len(B2_Ncells), "-2")
  colnames(A2) <- paste0("Gene", seq_len(Ngenes))

  B2 <- t(A2)

  ###############################################
  #   OUTPUT FILE           #####################
  ###############################################

  All_Colours <- c( t(B1ref_cols), t(B2ref_cols) )
  All_Shapes <- c(B1_shape, B2_shape)

  Batch_Simulated <- list("Batch_1" = B1, "Batch_2" = B2, "Batch_Effect" = Batch_Effect, "Colours" = All_Colours, "Shapes" = All_Shapes)

  return(Batch_Simulated)

}

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

##Plot_Results##
#Function to plot results
# INPUT :
#
# OUTPUT :
Plot_Results <- function( Data = NULL, Data_Type = NULL, Title = "Data Plot", X_Axis = "X", Y_Axis = "Y",
                          Size = 1, Colour = 1, Shape = 16 ){

  if(Data_Type == "tSNE"){
    df <- data.frame(x = Data$Y[,1], y = Data$Y[,2] )
    X_Axis = "tSNE 1"
    Y_Axis = "tSNE 2"
  }
  if(Data_Type == "UMAP"){
    df <- data.frame(x = Data[,1], y = Data[,2] )
    X_Axis = "UMAP 1"
    Y_Axis = "UMAP 2"
  }
  if(Data_Type == "PCA"){
    df <- data.frame(x = Data$x[,1], y = Data$x[,2] )
    X_Axis = "PCA 1"
    Y_Axis = "PCA 2"
  }

  ggplot(df) +
    geom_point(aes(x=x, y=y, colour=Colour),shape=Shape, size = Size) +
    ggtitle(Title) + xlab(X_Axis) + ylab(Y_Axis) + ggthemes::theme_base() +
    guides(colour = guide_legend(override.aes = list(size=3)))


}


##Get_Edges##
#Get a vector with the edges from a minimum spanning tree
# INPUT : Minimum spanning tree (MST object)
#
# OUTPUT : Vector with edges
Get_Edges <- function( Mst = NULL ){

  Edges <- NULL
  Num_Nodes <- ncol(Mst)

  for(Node in 1:(Num_Nodes-1)){

    Node_Edges <- (which(Mst[,Node] != 0))

    for(Edge in 1:length(Node_Edges)){

      if( Node_Edges[Edge] > Node)
        Edges <- rbind(Edges, c(Node, Node_Edges[Edge]))

    }
  }

  colnames(Edges) <- c("IN", "OUT")

  return(Edges)
}

##Rad2Deg##
#Convert radians to degrees
# INPUT : Angle in radiasn
#
# OUTPUT : Angle in degree
Rad2Deg <- function(Angle_rad = NULL){
  return(Angle_rad*180/(pi))
}

