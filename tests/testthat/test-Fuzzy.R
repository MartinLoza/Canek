set.seed(0)

Batches <- SimBatches$batches
z <- CorrectBatches(Batches, queNumCelltypes = 2, doCosNorm = TRUE, clusterMethod = "kmeans", debug = TRUE)
Fuzzy_Data <- z$`B2/B1`$`Correction Data`$`Fuzzy Data`

test_that("Fuzzy works", {
  expect_false(is.null(Fuzzy_Data))
  expect_equal(names(Fuzzy_Data), c("Fuzzy Memberships", "MST", "Fuzzied", "Edges Data"))
  expect_equal(length(which(is.na(Fuzzy_Data$`Fuzzy Memberships`))),0)
  expect_equal(length(Fuzzy_Data),4)
  expect_equal(min(rowSums(Fuzzy_Data$`Fuzzy Memberships`)), 1)
  expect_equal(max(rowSums(Fuzzy_Data$`Fuzzy Memberships`)), 1, tolerance = 1e-04)
  expect_equal(nrow(Fuzzy_Data$`Fuzzy Memberships`), ncol(z$`B2/B1`$`Query Batch (B2)`))
  expect_false(length(which(Fuzzy_Data$Fuzzied == TRUE)) != 0)
})
