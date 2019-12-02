context("test-Estimation")

set.seed(0)

Batches <- SimBatches$Batches
Pairs <- SimBatches$Pairs
z <- Canek:::EKF_BE(B1 = Batches[[1]], B2 = Batches[[2]], Pairs = Pairs, Sampling = TRUE)

test_that("Estimation works", {
  expect_false(is.null(z))
  expect_equal( names(z), c("Sampled Pairs", "Correction Vector"))

  expect_false(is.null(z$`Correction Vector`))
  expect_true(length(z$`Correction Vector`) == nrow(Batches$B1))
  expect_equal( length( which( is.finite(z$`Correction Vector`))), length(z$`Correction Vector`))
  expect_equal(z$`Correction Vector`[1], -10.92e-3, tolerance = 1e-4 )

  expect_equal(ncol(z$`Sampled Pairs`), 2)
  expect_equal(nrow(z$`Sampled Pairs`), 819)
  expect_equal(nrow(z$`Sampled Pairs`)/(nrow(Pairs)*0.2), 1.0, tolerance = 1e-3)
})
