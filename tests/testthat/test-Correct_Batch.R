context("test-Correct_Batch")

set.seed(0)
x <- matrix(sample(100, 50*30, replace = TRUE), ncol = 50)
y <- matrix(sample(100, 50*30, replace = TRUE), ncol = 50)

z <- Canek:::Correct_Batch(x,y,Dimensions = 10)
test_that("Correct_Batch works", {
  expect_false(is.null(z))
  expect_equal(length(z),4)
  expect_equal(names(z), c("Reference Batch (B1)", "Query Batch (B2)", "Corrected Query Batch", "Correction Data"))
  expect_equal(z$`Corrected Query Batch`[1,1],82.175386, tolerance = 0.0001)
})
