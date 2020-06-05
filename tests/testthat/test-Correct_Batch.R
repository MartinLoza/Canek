context("test-Correct_Batch")

set.seed(0)

x <- SimBatches$batches[[1]]
y <- SimBatches$batches[[2]]
z <- Canek:::Correct_Batch(x,y)

test_that("Correct_Batch works", {
  expect_false(is.null(z))
  expect_equal(length(which(is.na(z$`Corrected Query Batch`))), 0)
  expect_equal(length(z), 4)
  expect_equal(names(z), c("Reference Batch (B1)", "Query Batch (B2)", "Corrected Query Batch", "Correction Data"))
  expect_equal(z$`Corrected Query Batch`[1,1], 6.681331, tolerance = 0.0001)
  expect_equal(dim(z$`Corrected Query Batch`), dim(y) )
})
