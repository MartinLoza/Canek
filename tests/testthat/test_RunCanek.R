context("test-RunCanek")

# Create Seurat object.
x <- lapply(names(SimBatches$batches), function(batch) {
  Seurat::CreateSeuratObject(SimBatches$batches[[batch]], project = batch)
})
x <- merge(x[[1]], x[[2]])

# Create SingleCellExperiment object.
y <- Seurat::as.SingleCellExperiment(x)

# RunCanek.
x <- RunCanek(x, "orig.ident")
y <- RunCanek(y, "orig.ident")
z <- RunCanek(SimBatches$batches, debug = TRUE)

test_that("RunCanek works on Seurat objects", {
  expect_false(is.null(x))
  expect_is(x, "Seurat")
  expect_length(Seurat::Assays(x), 2)
  expect_equal(Seurat::Assays(x), c("RNA", "Canek"))
})

test_that("RunCanek works on SingleCellExperiment objects", {
  expect_false(is.null(y))
  expect_is(y, "SingleCellExperiment")
  expect_length(SummarizedExperiment::assays(y), 3)
  expect_equal(names(SummarizedExperiment::assays(y)), c("counts", "logcounts", "Canek"))
})

test_that("RunCanek works on lists", {
  expect_false(is.null(z))
  expect_is(z, "list")
  expect_length(z, 3)
  expect_equal(names(z), c("B2/B1", "Batches Integrated", "Total_Correction_Time"))
})
