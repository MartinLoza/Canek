set.seed(0)

Batches <- SimBatches$batches
x <- CorrectBatches(Batches, queNumCelltypes = 2, doCosNorm = TRUE, clusterMethod = "kmeans", debug = TRUE)
y <- CorrectBatches(Batches, doCosNorm = TRUE, clusterMethod = "louvain", debug = TRUE)
dataKmeans <- x$`B2/B1`$`Correction Data`$Clusters
dataLouvain <- y$`B2/B1`$`Correction Data`$Clusters

test_that("Clustering with kmeans works", {
  expect_false(is.null(dataKmeans))
  expect_equal(names(dataKmeans), c("cluster", "centers"))
  expect_equal(as.integer(table(dataKmeans$cluster)), 631)
  expect_equal(dim(dataKmeans$centers), c(1, 10))
})

test_that("Clustering with louvain works", {
  expect_false(is.null(dataLouvain))
  expect_equal(names(dataLouvain), c("cluster", "centers"))
  expect_equal(as.integer(table(dataLouvain$cluster)), c(178,146,162,145))
  expect_equal(dim(dataLouvain$centers), c(4, 10))
})
