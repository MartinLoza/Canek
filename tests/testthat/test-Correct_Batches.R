context("test-Correct_Batches")

set.seed(0)

Batches <- SimBatches
z <- Correct_Batches(Batches)

test_that("Correct_Batches works", {
  expect_false(is.null(z))
  expect_equal(length(which(is.na(z$`Batches Integrated`))),0)
  expect_equal(length(z),2)
  expect_equal(names(z), c("B1/B2","Batches Integrated"))
  expect_equal(z$`B1/B2`$`Corrected Query Batch`[1,1], 6.670866, tolerance = 0.0001)
  expect_equal( dim(z$`Batches Integrated`), dim(cbind(Batches$B1,Batches$B2) ))
})
