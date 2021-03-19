
CalculateMST <- function(x) {
  x <- as.matrix(dist(x))
  g <- igraph::graph_from_adjacency_matrix(x, mode = "undirected", weighted = "weight")
  igraph::mst(g)
}
