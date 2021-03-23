x <- iris[, -5]
mem <- iris$Species

centers <- Canek:::CalculateCenters(x, mem)

test_that("CalculateCenters", {
  expect_is(centers, "matrix")
  expect_equal(dim(centers), c(3, 4))
  expect_equal(rownames(centers), levels(mem))
  expect_equal(colnames(centers), colnames(x))
})

mst <- Canek:::CalculateMST(centers)
edges1 <- igraph::as_edgelist(mst, names = TRUE)
edges2 <- igraph::as_edgelist(mst, names = FALSE)

test_that("CalculateMST", {
  expect_is(mst, "igraph")
  expect_is(edges1, "matrix")
  expect_is(edges2, "matrix")
  expect_equal(edges1[1, 1], "setosa")
  expect_equal(edges2[1, 1], 1)
})

