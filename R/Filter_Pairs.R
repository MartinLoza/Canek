

#' Title Pairs_Selection
#'
#' Function to filter MNNs pairs
#'
#' @param B1 Reference batch single-cell data.
#' @param B2 Query's batch single-cell data.
#' @param Pairs A matrix containing MNNs pairs. First column corresponds to query-batch cell indexes.
#' @param Num_Clusters Number of clusters used to filter pairs.
#' @param Verbose Print output.
#'
#' @return A matrix containing the filtered pairs. First column corresponds to query-batch cell indexes.
#'
#' @details Function to select pairs to used on batch effect correction.
#'  The pairs are selected by clustering the query batch and analyzing the mean distance of pair cells on each cluster.
#'  The pairs of the cluster with minimum mean distance are returned.
#'
Pairs_Selection <- function(B1,B2, Pairs, Num_Clusters = 1, Verbose = FALSE){

  if(Verbose){
    cat("\n\nSelecting Pairs by Clusters")
    cat(paste('\n\tClusters for selecting pairs:',Num_Clusters))
  }

  #Initialization of varibles for clustering and selection of cell pairs
  Mu <- matrix(0L, nrow = Num_Clusters, ncol = 1)
  SD <- matrix(0L, nrow = Num_Clusters, ncol = 1)

  # Clustering to select pairs
  #Cluster <- kmeans(t(B2)[,1:10],Num_Clusters)
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
