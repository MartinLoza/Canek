
#' Title Fuzzy
#'
#' Function to score cell's memberships by fuzzy logic
#'
#' @param cluMem Memberships' clustering data.
#' @param pcaQue PCA representation of the cells.
#' @param corCell Matrix containing the initial membership assignment.
#' Matrix dimensions are expected as #Cell x #Memberships, with each row sum equal to 1.
#' @param verbose Print output.
#'
#' @details This function stablishes the fuzzification for the cells' membership.
#'  A minimum spanning tree (MST) is created among memberships, and the fuzzification is performed
#'   for each of the edges of the MST.#'
#'
Fuzzy <- function(cluMem = NULL, pcaQue = NULL, corCell = NULL, verbose = FALSE){

  #INIT
  Num_Cells <- nrow(pcaQue)
  Num_Memberships <- nrow(cluMem$centers)
  PCA_Max <- 2
  Fuzzied <- rep(FALSE, Num_Cells)
  Fuzzy_Memberships <- corCell
  Edges_Data <- list()


  # Create Minimum spanning tree (MST) by using centers of Memberships as nodes
  if(verbose)
    cat( '\n\tObtaining Minimum Spanning Tree' )

  Mst <- CalculateMST(cluMem$centers[, 1:PCA_Max])

  #Get edges from MST
  Edges <- igraph::as_edgelist(Mst, names = FALSE)

  #Fuzzy process
  if(verbose)
    cat('\n\tFuzzificating cells from each edge')

  for(Edge in 1:nrow(Edges)){


    #Init nodes
    IN_Membership <- list()
    OUT_Membership <- list()
    IN_Node <- Edges[Edge,1]
    OUT_Node <- Edges[Edge,2]
    IN_Node_PCA <- cluMem$centers[IN_Node,1:PCA_Max]
    OUT_Node_PCA <- cluMem$centers[OUT_Node,1:PCA_Max]

    #Translate OUT node according to IN node PCA coordinates
    OUT_Node_PCA_Transformed <- OUT_Node_PCA - IN_Node_PCA

    #Find angle between the current nodes (we set it negative because we would correct this angle in further steps)
    Alpha <- (-GetRotationAngle ( OUT_Node_PCA_Transformed ) )
    names(Alpha)<-"Alpha"

    #Get cells from both IN and OUT memberships
    IN_Membership_Cells_Index <- which(cluMem$cluster == IN_Node)
    IN_Membership_Cells <- pcaQue[IN_Membership_Cells_Index, 1:PCA_Max ]
    rownames(IN_Membership_Cells) <- IN_Membership_Cells_Index

    OUT_Membership_Cells_Index <- which(cluMem$cluster == OUT_Node)
    OUT_Membership_Cells <- pcaQue[which(cluMem$cluster == OUT_Node), 1:PCA_Max ]
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

      if(Cells_Filtered[Cell,1] < 0){
        IN_Fuzzification <- 0
        OUT_Fuzzification <- 1
      }else if(Cells_Filtered[Cell,1] > OUT_Node_PCA_Transformed[1]){
        IN_Fuzzification <- 0
        OUT_Fuzzification <- 1
      }else{
        IN_Fuzzification <- 1 + (Slope*Cells_Filtered[Cell,1])
        OUT_Fuzzification <- ( 1 - IN_Fuzzification )
      }

      Fuzzification <- rbind( Fuzzification, c(IN_Fuzzification, OUT_Fuzzification) )

      #If this cells has not been fuzzified before set fuzzification values
      if(Fuzzied[Cells_Filtered_RowNames[Cell]] == FALSE){

        Fuzzy_Memberships[Cells_Filtered_RowNames[Cell], IN_Node] <- IN_Fuzzification
        Fuzzy_Memberships[Cells_Filtered_RowNames[Cell], OUT_Node] <- OUT_Fuzzification
        Fuzzied[Cells_Filtered_RowNames[Cell]] = TRUE

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

  Fuzzy_Data <- list( "Fuzzy Memberships" = Fuzzy_Memberships, "MST" = Mst, "Fuzzied" =   Fuzzied, "Edges Data" = Edges_Data)

  return(Fuzzy_Data)
}

##GetRotationAngle##
#Get correction angle according to quadrant
# INPUT : Node PCA coordinates
#
# OUTPUT : Angle in radians
GetRotationAngle <- function( PCA_Coordinates = NULL ){

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

#' CheckZeroCV
#'
#' @param MST Minimum Spanning Tree
#' @param cluMem Clusters used on MST
#' @param memCorrData Data to correct
#' @param corGene Data to correct
#' @param zeroCorrection Vector indicating which membership has a zero correction vector
#'
CheckZeroCV <-function(MST = NULL, cluMem = NULL,
                       memCorrData = NULL, corGene = NULL,
                       zeroCorrection = NULL){

  Is_Zero <- which(zeroCorrection == TRUE)
  Cluster_Dist <- as.matrix(dist(cluMem$centers,upper = TRUE))

  Node = 1
  while(length(Is_Zero) != 0){

    Related_Edges <- which(MST[Is_Zero[Node],] !=0)
    Related_Edges_No_Zero <- Related_Edges[which(zeroCorrection[Related_Edges] == FALSE)]
    if(length(Related_Edges_No_Zero) != 0){
      #if there are various, we select the one with the minimum distance
      if(length(Related_Edges_No_Zero) != 1){
        Related_Edges_No_Zero <- which(Cluster_Dist[Is_Zero[Node],] == min(Cluster_Dist[Is_Zero[Node],Related_Edges_No_Zero]))
      }
      #Assign correction vector
      memCorrData[[Is_Zero[Node]]]$`Correction Vector` <- memCorrData[[Related_Edges_No_Zero]]$`Correction Vector`
      corGene[,Is_Zero[Node]] <- memCorrData[[Related_Edges_No_Zero]]$`Correction Vector`
      zeroCorrection[Is_Zero[Node]] <- FALSE

      Node = 1
      Is_Zero <- which(zeroCorrection == TRUE)

    }else{ #If we don't find any related node with no zero correction vector, we analize the next node
      Node = Node + 1
    }
  }

  return(list("memCorrData" = memCorrData,
              "Correction_Matrix" = corGene))
}



