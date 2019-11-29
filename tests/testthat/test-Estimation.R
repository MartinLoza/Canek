context("test-Estimation")

set.seed(0)

Batches <- SimBatches
z <- Correct_Batches(Batches)
M1_Data <- z$`B1/B2`$`Correction Data`$`Membership Data`$`Membership Correction Data`$`Membership 1`
M2_Data <- z$`B1/B2`$`Correction Data`$`Membership Data`$`Membership Correction Data`$`Membership 2`
M3_Data <- z$`B1/B2`$`Correction Data`$`Membership Data`$`Membership Correction Data`$`Membership 3`

test_that("Estimation works", {
  expect_false(is.null(M1_Data) || is.null(M2_Data) || is.null(M3_Data) )
  expect_true( (length(M1_Data) == 4) &&  (length(M2_Data) == 4) &&  (length(M3_Data) == 4) )
  expect_equal( names(M1_Data), c("Cells Index", "Pairs Selection Data", "Sampled MNN Pairs", "Correction Vector"  ) )
  expect_equal( names(M2_Data), c("Cells Index", "Pairs Selection Data", "Sampled MNN Pairs", "Correction Vector"  ) )
  expect_equal( names(M3_Data), c("Cells Index", "Pairs Selection Data", "Sampled MNN Pairs", "Correction Vector"  ) )

  expect_false(is.null(M1_Data$`Correction Vector`) || is.null(M2_Data$`Correction Vector`) || is.null(M3_Data$`Correction Vector`) )
  expect_true( length(M1_Data$`Correction Vector`) == nrow(z$`Batches Integrated`) )
  expect_true( (length(M1_Data$`Correction Vector`) == length(M2_Data$`Correction Vector`)) &&
                 (length(M1_Data$`Correction Vector`) == length(M3_Data$`Correction Vector`)) )

  expect_equal( length( which( is.finite(M1_Data$`Correction Vector`) ) ), length(M1_Data$`Correction Vector`) )
  expect_equal( length( which( is.finite(M2_Data$`Correction Vector`) ) ), length(M2_Data$`Correction Vector`) )
  expect_equal( length( which( is.finite(M3_Data$`Correction Vector`) ) ), length(M3_Data$`Correction Vector`) )

  expect_equal( length(which(M1_Data$`Correction Vector` == 0)), length(M1_Data$`Correction Vector`) )
  expect_equal( M2_Data$`Correction Vector`[1], 48.49e-3, tolerance = 1e-4 )
  expect_equal( M3_Data$`Correction Vector`[1], -23.89e-3, tolerance = 1e-4 )


})
