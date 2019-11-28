
#' Title Fuzzy
#'
#' Function to score cell's memberships by fuzzy logic
#'
#' @param Cluster_Membership Memberships' clustering data
#' @param Cells_PCA PCA representation of the cells
#' @param Correction_Memberships Matrix containing the initial membership assignment. Matrix dimensions are expected as #Cell x #Memberships, with each row sum equal to 1.
#'
#' @details This function stablishes the fuzzification for the cells' membership. A minimum spanning tree (MST) is created among memberships, and the fuzzification is performed for each of the edges of the MST.
#'
#'
#' @examples
Fuzzy <- function(Cluster_Membership = NULL, Cells_PCA = NULL, Correction_Memberships = NULL){

  #INIT
  Num_Cells <- nrow(Cells_PCA)
  Num_Memberships <- nrow(Cluster_Membership$centers)
  PCA_Max <- 2
  Fuzzied <- rep(FALSE, Num_Cells)
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
