
##GetEdges##
#Get a vector with the edges from a minimum spanning tree
# INPUT : Minimum spanning tree (MST object)
#
# OUTPUT : Vector with edges
GetEdges <- function( Mst = NULL ){

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
