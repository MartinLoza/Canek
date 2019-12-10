context("test-RunCanek")

x <- lapply(names(SimBatches$Batches), function(batch) {
  Seurat::CreateSeuratObject(SimBatches$Batches[[batch]], project = batch)
})
x <- merge(x[[1]], x[[2]])
x <- RunCanek(x, "orig.ident")

test_that("Correct_Batch works", {
  expect_false(is.null(x))
  expect_is(x, "Seurat")
  expect_equal(length(x@assays), 2)
  expect_equal(names(x@assays), c("RNA", "Canek"))
})
