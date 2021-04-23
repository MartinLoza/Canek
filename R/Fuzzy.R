
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

  Mst <- CalculateMST(cluMem$centers[, 1:fuzzyPCA])

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
    IN_Node_PCA <- cluMem$centers[IN_Node,1:fuzzyPCA]
    OUT_Node_PCA <- cluMem$centers[OUT_Node,1:fuzzyPCA]

    #Translate OUT node according to IN node PCA coordinates
    OUT_Node_PCA_Transformed <- OUT_Node_PCA - IN_Node_PCA

    #Find angle between the current nodes (we set it negative because we would correct this angle in further steps)
    Alpha <- (-GetRotationAngle ( OUT_Node_PCA_Transformed ) )
    names(Alpha)<-"Alpha"

    #Get cells from both IN and OUT memberships
    IN_Membership_Cells_Index <- which(cluMem$cluster == IN_Node)
    IN_Membership_Cells <- pcaQue[IN_Membership_Cells_Index, 1:fuzzyPCA ]
    rownames(IN_Membership_Cells) <- IN_Membership_Cells_Index

    OUT_Membership_Cells_Index <- which(cluMem$cluster == OUT_Node)
    OUT_Membership_Cells <- pcaQue[which(cluMem$cluster == OUT_Node), 1:fuzzyPCA ]
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

    IN_Membership_Cells_Filtered <- IN_Membership_Cells_Transformed[IN_Membership_Cells_Transformed_Selected_Index, 1:fuzzyPCA, drop = FALSE]
    OUT_Membership_Cells_Filtered <- OUT_Membership_Cells_Transformed[OUT_Membership_Cells_Transformed_Selected_Index, 1:fuzzyPCA, drop = FALSE]

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

FuzzyNew <- function(cluMem = NULL, pcaQue = NULL, corCell = NULL, fuzzyPCA = 2, verbose = FALSE){

  #INIT
  nCells <- nrow(pcaQue)
  nMem <- nrow(cluMem$centers)
  #PCA_Max <- 2 # for tests
  Fuzzied <- rep(FALSE, nCells)
  Edges_Data <- list()
  corCell <- as.data.frame(corCell)
  corCell[["Fuzzified"]] <- FALSE

  # Create Minimum spanning tree (MST) by using centers of Memberships as nodes
  if(verbose)
    cat( '\n\tObtaining Minimum Spanning Tree' )

  mst <- CalculateMST(cluMem$centers[,1:fuzzyPCA])

  #Get edges from MST
  edges <- igraph::as_edgelist(mst, names = FALSE)

  # now we analize one edge
  #edge = 1
  for(edge in seq_len(nrow(edges))){

    #Get the memberships of the selected edge
    inNode <- edges[edge,1]
    outNode <- edges[edge,2]

    #Get the cells related with the memberships
    idxCells <- which( cluMem$cluster == inNode | cluMem$cluster == outNode )
    edgeCells <- pcaQue[idxCells, 1:fuzzyPCA, drop = FALSE]

    #Get the centers of the memberships
    inMemCenter <- cluMem$centers[inNode, 1:fuzzyPCA, drop = FALSE]
    outMemCenter <- cluMem$centers[outNode, 1:fuzzyPCA, drop = FALSE]

    #Translate the cells and the centers
    edgeCells <- sweep(edgeCells, 2, inMemCenter, "-")
    outMemCenter <- outMemCenter - inMemCenter
    inMemCenter <- inMemCenter - inMemCenter

    #Obtain the norm2 and the unit vector of the edge
    #On this case, because we already translate the vectors, the edge vector corresponds to the outMemCenter
    normEV <- norm(outMemCenter, type = "2")
    unitEV <- outMemCenter/normEV

    # Obtain the components of the other vectors
    edgeCellsComp <- edgeCells %*% t(unitEV)

    # Compare with the edge's norm
    edgeCellsComp <- edgeCellsComp/normEV

    #Fuzzy comparison
    for(cell in seq_len(length(edgeCellsComp))){

      iCell <- idxCells[cell]

      if(edgeCellsComp[cell] >= 0 & edgeCellsComp[cell] <= 1){ #if the cell's projections is withing the edge, we assign a fuzzy score
        if(corCell$Fuzzified[iCell] == FALSE){
          corCell[iCell,outNode] <- edgeCellsComp[cell]
          corCell[iCell,inNode] <- 1 - edgeCellsComp[cell]
        }else{
          corCell[iCell,outNode] <-mean(corCell[iCell,outNode], edgeCellsComp[cell])
          corCell[iCell,inNode] <- mean(corCell[iCell,inNode], 1 - edgeCellsComp[cell])
        }
      }else{ #for now we just assign the cell's membership score, a possible improvement to check would be to change teh score to the other membership
        if(corCell$Fuzzified[iCell] == TRUE){
          if(outNode == cluMem$cluster[iCell]){
            corCell[iCell,outNode] <-mean(corCell[iCell,outNode], 1)
            corCell[iCell,inNode] <- mean(corCell[iCell,inNode], 0)
          }else{
            corCell[iCell,outNode] <-mean(corCell[iCell,outNode], 0)
            corCell[iCell,inNode] <- mean(corCell[iCell,inNode], 1)
          }
        }
      }

      corCell$Fuzzified[iCell] <- TRUE
    }

    # prepare debug info
    IN_Membership <- list("Cells" = NULL,
                          "Transformed" = NULL,
                          "Filtered" = NULL)
    OUT_Membership <- list("Cells" = NULL, "Transformed" = NULL,
                           "Filtered" = NULL)

    Edges_Data[[paste("Edge-", edge, sep = "")]] <- list( "IN Node" = inNode, "OUT Node" = outNode,
                                                          "Angle" = NULL, "IN-Mem-Cells Data" = IN_Membership,
                                                          "OUT-Mem-Cells Data" = OUT_Membership, "Slope" = NULL,
                                                          "Fuzzification" = NULL)
  }

  Fuzzied <- corCell[,ncol(corCell)]
  corCell <- as.matrix(corCell[,-ncol(corCell)])

  Fuzzy_Data <- list("Fuzzy Memberships" = corCell, "MST" = mst,
                     "Fuzzied" = Fuzzied, "Edges Data" = Edges_Data)

  return(Fuzzy_Data)
}

FuzzyNew2 <- function(cluMem = NULL, pcaQue = NULL, corCell = NULL, fuzzyPCA = 2, verbose = FALSE){

  #INIT
  nCells <- nrow(pcaQue)
  nMem <- nrow(cluMem$centers)
  #fuzzyPCA <- 2 # for tests
  Fuzzied <- rep(FALSE, nCells)
  Edges_Data <- list()
  corCell <- as.data.frame(corCell)
  corCell[["Fuzzified"]] <- FALSE

  # Create Minimum spanning tree (MST) by using centers of Memberships as nodes
  if(verbose)
    cat( '\n\tObtaining Minimum Spanning Tree' )

  mst <- CalculateMST(cluMem$centers[,1:fuzzyPCA])

  #Get edges from MST
  edges <- igraph::as_edgelist(mst, names = FALSE)

  # now we analize one edge
  #edge = 1
  for(edge in seq_len(nrow(edges))){

    #Get the memberships of the selected edge
    inNode <- edges[edge,1]
    outNode <- edges[edge,2]

    #Get the cells related with the memberships
    idxCells <- which( cluMem$cluster == inNode | cluMem$cluster == outNode )
    edgeCells <- pcaQue[idxCells, 1:fuzzyPCA, drop = FALSE]

    #Get the centers of the memberships
    inMemCenter <- cluMem$centers[inNode, 1:fuzzyPCA, drop = FALSE]
    outMemCenter <- cluMem$centers[outNode, 1:fuzzyPCA, drop = FALSE]

    #Translate the cells and the centers
    edgeCells <- sweep(edgeCells, 2, inMemCenter, "-")
    outMemCenter <- outMemCenter - inMemCenter
    inMemCenter <- inMemCenter - inMemCenter

    #Obtain the norm2 and the unit vector of the edge
    #On this case, because we already translate the vectors, the edge vector corresponds to the outMemCenter
    normEV <- norm(outMemCenter, type = "2")
    unitEV <- outMemCenter/normEV

    # Obtain the components of the other vectors
    edgeCellsComp <- edgeCells %*% t(unitEV)

    # new test
    # obtenemos el minimo
    minComp <- min(edgeCellsComp, na.rm = TRUE)

    # trasladar los componentes respecto al minimo (para comparar con el minimo y maximo en vez de con los centros)
    edgeCellsComp <- edgeCellsComp - minComp

    # comparar con el maximo
    edgeCellsComp <- edgeCellsComp/max(edgeCellsComp, na.rm = TRUE)

    #TODO: ver si esto se puede hacer, el problema es que necestiamos ver si fue fuzzificada para hacer la media
    # assign memberships
    # corCell[idxCells, outNode] <- edgeCellsComp
    # corCell[idxCells, inNode] <- 1 - edgeCellsComp

    #Fuzzy comparison
    for(cell in seq_len(length(edgeCellsComp))){

      iCell <- idxCells[cell]

      if(corCell$Fuzzified[iCell] == FALSE){
        corCell[iCell,outNode] <- edgeCellsComp[cell]
        corCell[iCell,inNode] <- 1 - edgeCellsComp[cell]
      }else{
        corCell[iCell,outNode] <-mean(corCell[iCell,outNode], edgeCellsComp[cell])
        corCell[iCell,inNode] <- mean(corCell[iCell,inNode], 1 - edgeCellsComp[cell])
      }
      corCell$Fuzzified[iCell] <- TRUE
    }

    # prepare debug info
    IN_Membership <- list("Cells" = NULL,
                          "Transformed" = NULL,
                          "Filtered" = NULL)
    OUT_Membership <- list("Cells" = NULL, "Transformed" = NULL,
                           "Filtered" = NULL)

    Edges_Data[[paste("Edge-", edge, sep = "")]] <- list( "IN Node" = inNode, "OUT Node" = outNode,
                                                          "Angle" = NULL, "IN-Mem-Cells Data" = IN_Membership,
                                                          "OUT-Mem-Cells Data" = OUT_Membership, "Slope" = NULL,
                                                          "Fuzzification" = NULL)
  }

  Fuzzied <- corCell[,ncol(corCell)]
  corCell <- as.matrix(corCell[,-ncol(corCell)])

  Fuzzy_Data <- list("Fuzzy Memberships" = corCell, "MST" = mst,
                     "Fuzzied" = Fuzzied, "Edges Data" = Edges_Data)

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

  names(zeroCorrection) <- 1:ncol(corGene)
  colnames(corGene) <- 1:ncol(corGene)

  isZero <- which(zeroCorrection == TRUE)
  Cluster_Dist <- as.matrix(dist(cluMem$centers,upper = TRUE))
  adjMST <- igraph::as_adjacency_matrix(MST)

  idx = 1
  while(length(isZero) != 0){

    Node <- as.character(isZero[idx])

    Related_Edges <- names(which(adjMST[Node,] != 0))
    Related_Edges_No_Zero <- Related_Edges[which(zeroCorrection[Related_Edges] == FALSE)]
    if(length(Related_Edges_No_Zero) != 0){
      #if there are various, we select the one with the minimum distance
      if(length(Related_Edges_No_Zero) > 1){
        Related_Edges_No_Zero <- which(Cluster_Dist[Node,] == min(Cluster_Dist[Node,Related_Edges_No_Zero]))
      }
      #Assign correction vector
      memCorrData[[as.integer(Node)]]$`Correction Vector` <- memCorrData[[as.integer(Related_Edges_No_Zero)]]$`Correction Vector`
      corGene[,Node] <- memCorrData[[as.integer(Related_Edges_No_Zero)]]$`Correction Vector`
      zeroCorrection[Node] <- FALSE

      idx = 1
      isZero <- which(zeroCorrection == TRUE)

    }else{ #If we don't find any related node with no zero correction vector, we analize the next node
      idx = idx + 1
    }
  }

  return(list("memCorrData" = memCorrData,
              "Correction_Matrix" = corGene))
}



