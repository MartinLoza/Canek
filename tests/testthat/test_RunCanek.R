context("test-RunCanek")

x <- lapply(names(SimBatches$batches), function(batch) {
  Seurat::CreateSeuratObject(SimBatches$batches[[batch]], project = batch)
})
x <- merge(x[[1]], x[[2]])
x <- RunCanek(x, "orig.ident")

test_that("RunCanek works on Seurat objects", {
  expect_false(is.null(x))
  expect_is(x, "Seurat")
  expect_equal(Seurat::Assays(x), 2)
  expect_equal(Seurat::Assays(x), c("RNA", "Canek"))
})
  expect_equal(length(x@assays), 2)
  expect_equal(names(x@assays), c("RNA", "Canek"))
})
