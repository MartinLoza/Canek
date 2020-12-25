context("test-FilterPairs")

set.seed(0)

Batches <- SimBatches$batches
z <- CorrectBatch(Batches$B1, Batches$B2, queNumCelltypes = 2, pairsFilter = TRUE)
pairs <- z$`Correction Data`$`MNN Pairs`

test_that("FilterPairs works", {
  expect_false(is.null(pairs))
  expect_true((nrow(pairs) == 5918) & (ncol(pairs) == 2))
  expect_true(length(which(is.na(pairs))) == 0)
  expect_true(length(unique(pairs[,1])) == 674 & length(unique(pairs[,2])) == 514)
  expect_true(length(which(pairs[,1] == 1)) == 8)
  expect_true(length(which(pairs[,1] == 689)) == 0)
})
