context("test-CorrectBatch")

set.seed(0)

x <- SimBatches$batches[[1]]
y <- SimBatches$batches[[2]]
z <- Canek:::CorrectBatch(x, y, doCosNorm = TRUE, clusterMethod = "kmeans")

test_that("CorrectBatch works", {
  expect_false(is.null(z))
  expect_equal(length(which(is.na(z$`Corrected Query Batch`))), 0)
  expect_equal(length(z), 6)
  expect_equal(names(z), c("Reference Batch (B1)", "Query Batch (B2)", "Corrected Query Batch", "Correction Data", "Correction_Time", "debug"))
  expect_equal(z$`Corrected Query Batch`[1,1], 6.652662, tolerance = 0.0001)
  expect_equal(dim(z$`Corrected Query Batch`), dim(y) )
})
