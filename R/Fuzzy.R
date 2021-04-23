
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
Fuzzy <- function(cluMem = NULL, pcaQue = NULL, corCell = NULL, fuzzyPCA = 10, verbose = FALSE){

  #INIT
  nCells <- nrow(pcaQue)
  nMem <- nrow(cluMem$centers)
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

    # Using the minimum and maximum to compare the components
    ## Translate the components according with the minimum
    edgeCellsComp <- edgeCellsComp - minComp

    ## Compare with the maximum
    edgeCellsComp <- edgeCellsComp/max(edgeCellsComp, na.rm = TRUE)

    #TODO: this could be a future improvement, but we have to find a different way to check the Fuzzified flag.
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
                                                          "Fuzzification" = NULL, "CellsComponents" = edgeCellsComp)
  }

  Fuzzied <- corCell[,ncol(corCell)]
  corCell <- as.matrix(corCell[,-ncol(corCell)])

  Fuzzy_Data <- list("Fuzzy Memberships" = corCell, "MST" = mst,
                     "Fuzzied" = Fuzzied, "Edges Data" = Edges_Data)

  return(Fuzzy_Data)
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



