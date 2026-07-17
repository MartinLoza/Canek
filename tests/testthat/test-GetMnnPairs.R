x <- matrix(c(1,3,1,1,2,3), ncol = 2, byrow = TRUE)
y <- matrix(c(2,1,10,9,20,18,19,15), ncol = 2, byrow = TRUE)

mnn <- Canek:::GetMnnPairs(t(x), t(y), kNN = 1)

test_that("GetMnnPairs works", {
  expect_true(!is.null(mnn))
  expect_named(mnn, "Pairs")
  expect_equal(dim(mnn$Pairs), c(1,2))
  expect_equal(colnames(mnn$Pairs), c("queBatch-Cells-Index", "refBatch-Cells-Index"))
  expect_equal(unname(mnn$Pairs[1,]), c(1,2))
})

mnn <- Canek:::GetMnnPairs(t(x), t(y), kNN = 2)

test_that("Get_MNN_Pairs works", {
  expect_true(!is.null(mnn))
  expect_named(mnn, "Pairs")
  expect_equal(dim(mnn$Pairs), c(4, 2))
  expect_equal(colnames(mnn$Pairs), c("queBatch-Cells-Index", "refBatch-Cells-Index"))
  expect_equal(unname(mnn$Pairs[1,]), c(1,2))
  expect_equal(unname(mnn$Pairs[2,]), c(1,3))
  expect_equal(unname(mnn$Pairs[3,]), c(2,1))
  expect_equal(unname(mnn$Pairs[4,]), c(2,3))

})

test_that("FindMnnPairs matches an independent reference implementation on a larger random example", {
  set.seed(123)
  nRef <- 40
  nQue <- 35
  kNN <- 5

  B1_B2_NN <- cbind(rep(seq_len(nRef), each = kNN), sample(seq_len(nQue), nRef * kNN, replace = TRUE))
  B2_B1_NN <- cbind(rep(seq_len(nQue), each = kNN), sample(seq_len(nRef), nQue * kNN, replace = TRUE))
  colnames(B1_B2_NN) <- c("Batch-1", "Batch-2")
  colnames(B2_B1_NN) <- c("Batch-2", "Batch-1")

  # Independent reference: a pair is mutual when the same (ref, que)
  # combination appears in both NN tables, so a relational join on both
  # columns is exactly the definition of a mutual nearest neighbor. Uses
  # base R's merge() rather than which()/split() grouping, so it can't
  # share a bug with either FindMnnPairs implementation. merge() also
  # naturally reproduces the row multiplicity from duplicate NN entries,
  # matching the nested loop in FindMnnPairs.
  ReferenceMnnPairs <- function(B1_B2_NN, B2_B1_NN) {
    df1 <- data.frame(ref = B1_B2_NN[, 1], que = B1_B2_NN[, 2])
    df2 <- data.frame(que = B2_B1_NN[, 1], ref = B2_B1_NN[, 2])
    merged <- merge(df1, df2, by = c("ref", "que"))
    as.matrix(merged[, c("que", "ref")])
  }

  standardize <- function(m) {
    m <- unname(m)
    m[order(m[, 1], m[, 2]), ]
  }

  expected <- standardize(ReferenceMnnPairs(B1_B2_NN, B2_B1_NN))
  result <- standardize(Canek:::FindMnnPairs(B1_B2_NN = B1_B2_NN, B2_B1_NN = B2_B1_NN, B2_NCells = nQue)$Pairs)

  expect_equal(result, expected)
})
