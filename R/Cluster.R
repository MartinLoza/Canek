ClusterKMeans <- function(x, maxMem = 10, nMem = NULL, verbose = TRUE, usepam = FALSE) {
  if(is.null(nMem)){
    if(verbose)
      cat("\n\nFinding number of memberships")

    nMem <- pamk(x,
                 krange = seq_len(maxMem),
                 usepam = usepam)$nc

    if (verbose)
      cat(paste('\n\tNumber of memberships found:', nMem))
  }
  res <- kmeans(x, nMem)
  res <- list(cluster = res$cluster, centers = res$centers)

  list(result = res, nMem = nMem)
}

ClusterLouvain <- function(x, k = 10, verbose = TRUE) {
  if (verbose)
    cat("\n\nFinding number of memberships")

  g <- bluster::makeSNNGraph(x, k = k)
  res <- igraph::cluster_louvain(g)

  memberships <- igraph::membership(res)
  nMem <- length(unique(memberships))

  if (verbose)
    cat(paste('\n\tNumber of memberships found:', nMem))

  centers <- CalculateCenters(x, memberships)

  res <- list(cluster = memberships, centers = centers)
  list(result = res, nMem = nMem)
}

CalculateCenters <- function(x, memberships) {
  lMem <- unique(memberships)
  nMem <- length(lMem)


  centers <- sapply(seq_len(nMem), function(n) {
    membership <- lMem[n]
    colMeans(x[memberships == membership, ])
  })
  colnames(centers) <- seq_len(nMem)
  t(centers)
}
