context("test-Filter_Pairs")

set.seed(0)

Batches <- SimBatches$batches
z <- Correct_Batches(Batches, Num_Clusters = 3, Sampling = FALSE)
M2_Pairs_Data <- z$`B1/B2`$`Correction Data`$`Membership Data`$`Membership Correction Data`$`Membership 2`$`Pairs Selection Data`
M3_Pairs_Data <- z$`B1/B2`$`Correction Data`$`Membership Data`$`Membership Correction Data`$`Membership 3`$`Pairs Selection Data`

test_that("Filter_Pairs works", {
  expect_false(is.null(M2_Pairs_Data) || is.null(M3_Pairs_Data))
  expect_true( (length(M2_Pairs_Data) == 3) &&  (length(M3_Pairs_Data) == 3) )
  expect_equal( names(M2_Pairs_Data), c("Selected Pairs", "Clusters", "Selected Cluster") )
  expect_equal( names(M3_Pairs_Data), c("Selected Pairs", "Clusters", "Selected Cluster") )
  expect_equal(length(which(is.na(M2_Pairs_Data$`Selected Pairs`))),0)
  expect_equal(length(which(is.na(M3_Pairs_Data$`Selected Pairs`))),0)
  expect_equal(dim(M2_Pairs_Data$`Selected Pairs`), c(199,2))
  expect_equal(dim(M3_Pairs_Data$`Selected Pairs`), c(775,2))

  expect_true( (M2_Pairs_Data$`Selected Cluster` == 2) && (M3_Pairs_Data$`Selected Cluster` == 1) )

  expect_false(is.null(M2_Pairs_Data$Clusters) || is.null(M3_Pairs_Data$Clusters))
})
