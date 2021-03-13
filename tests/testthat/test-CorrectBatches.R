context("test-CorrectBatches")

set.seed(0)

Batches <- SimBatches$batches
z <- CorrectBatches(Batches, doCosNorm = TRUE, clusterMethod = "kmeans")

test_that("CorrectBatches works", {
  expect_false(is.null(z))
  expect_equal(length(which(is.na(z))),0)
  expect_equal(dim(z), dim(cbind(Batches$B1,Batches$B2)))
  expect_true(ncol(z) == 1579 & nrow(z) == 500)
  expect_equal(z[1,1], 6.664089, tolerance = 0.0001)
})
