context("test-Estimation")

set.seed(0)

Batches <- SimBatches$batches
Pairs <- SimBatches$pairs
x <- Canek:::EkfBE(refBatch = Batches[[1]], queBatch = Batches[[2]], pairs = Pairs, sampling = TRUE)
y <- Canek:::MedianBE(refBatch = Batches[[1]], queBatch = Batches[[2]], pairs = Pairs)

test_that("EKF Method", {
  expect_false(is.null(x))
  expect_equal( names(x), c("Correction Vector", "Sampled Pairs"))

  expect_false(is.null(x$`Correction Vector`))
  expect_true(length(x$`Correction Vector`) == nrow(Batches$B1))
  expect_equal(length(which(is.finite(x$`Correction Vector`))), length(x$`Correction Vector`))
  expect_equal(x$`Correction Vector`[1], -0.01188622, tolerance = 1e-4 )

  expect_equal(ncol(x$`Sampled Pairs`), 2)
  expect_equal(nrow(x$`Sampled Pairs`), 819)
  expect_equal(nrow(x$`Sampled Pairs`)/(nrow(Pairs)*0.2), 1.0, tolerance = 1e-3)
})

test_that("Median Method", {
  expect_false(is.null(y))
  expect_equal(names(y), c("Correction Vector", "Sampled Pairs"))

  expect_false(is.null(y$`Correction Vector`))
  expect_true(length(y$`Correction Vector`) == nrow(Batches$B1))
  expect_equal(length(which(is.finite(y$`Correction Vector`))),length(y$`Correction Vector`))
  expect_equal(y$`Correction Vector`[1],0.00523153, tolerance = 1e-4 )

  expect_equal(y$`Sampled Pairs`, NULL)
})
