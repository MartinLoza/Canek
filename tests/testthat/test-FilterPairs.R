context("test-FilterPairs")

set.seed(0)

Batches <- SimBatches$batches
z <- CorrectBatches(Batches, queNumCelltypes = 2, pairsFilter = TRUE, doCosNorm = TRUE, clusterMethod = "kmeans", debug = TRUE)
pairs <- z$`B2/B1`$`Correction Data`$`MNN Pairs`

test_that("FilterPairs works", {
  expect_false(is.null(pairs))
  expect_equal(nrow(pairs), 5918)
  expect_equal(ncol(pairs), 2)
  expect_equal(length(which(is.na(pairs))), 0)
  expect_equal(length(unique(pairs[, 1])), 514)
  expect_equal(length(unique(pairs[, 2])), 674)
  expect_equal(length(which(pairs[,2] == 1)), 8)
  expect_equal(length(which(pairs[,2] == 689)), 0)
})
