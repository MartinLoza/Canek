x <- matrix(c(1,3,1,1,2,3), ncol = 2, byrow = TRUE)
y <- matrix(c(2,1,10,9,20,18,19,15), ncol = 2, byrow = TRUE)

mnn <- Canek:::Get_MNN_Pairs(t(x), t(y), kNN = 1)

test_that("Get_MNN_Pairs works", {
  expect_true(!is.null(mnn))
  expect_named(mnn, "Pairs")
  expect_equal(dim(mnn$Pairs), c(1, 2))
  expect_equal(colnames(mnn$Pairs), c("B2-Cells-Index", "B1-Cells-Index"))
  expect_equal(unname(mnn$Pairs[1, ]), c(1, 2))
})

mnn <- Canek:::Get_MNN_Pairs(t(x), t(y), kNN = 2)

test_that("Get_MNN_Pairs works", {
  expect_true(!is.null(mnn))
  expect_named(mnn, "Pairs")
  expect_equal(dim(mnn$Pairs), c(4, 2))
  expect_equal(colnames(mnn$Pairs), c("B2-Cells-Index", "B1-Cells-Index"))
  expect_equal(unname(mnn$Pairs[1, ]), c(1, 2))
  expect_equal(unname(mnn$Pairs[2, ]), c(1, 3))
  expect_equal(unname(mnn$Pairs[3, ]), c(2, 1))
  expect_equal(unname(mnn$Pairs[4, ]), c(2, 3))

})
