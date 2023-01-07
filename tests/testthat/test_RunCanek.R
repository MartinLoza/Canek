# Create Seurat object.
x <- lapply(names(SimBatches$batches), function(batch) {
  Seurat::CreateSeuratObject(SimBatches$batches[[batch]], project = batch)
})
x <- merge(x[[1]], x[[2]])
cellnames <- colnames(x)

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
  expect_equal(colnames(x), cellnames)
})

test_that("RunCanek works on SingleCellExperiment objects", {
  expect_false(is.null(y))
  expect_is(y, "SingleCellExperiment")
  expect_length(SummarizedExperiment::assays(y), 3)
  expect_equal(names(SummarizedExperiment::assays(y)), c("counts", "logcounts", "Canek"))
  expect_equal(colnames(y), cellnames)
})

test_that("RunCanek works on lists", {
  expect_false(is.null(z))
  expect_is(z, "list")
  expect_length(z, 3)
  expect_equal(names(z), c("B2/B1", "Batches Integrated", "Total_Correction_Time"))
  expect_equal(colnames(z$`Batches Integrated`), cellnames)
  expect_error(CorrectBatches(list(B1 = SimBatches$batches$B1, B2 = SimBatches$batches$B1)))
})

x <- RunCanek(x, "orig.ident", integration.name = "CanekRNA")
y <- RunCanek(y, "orig.ident", integration.name = "CanekRNA")
test_that("Setting RunCanek integration.name argument works", {
  expect_true("CanekRNA" %in% names(x))
  expect_true("CanekRNA" %in% names(SummarizedExperiment::assays(y)))
})
