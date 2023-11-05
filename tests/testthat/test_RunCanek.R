# Get matrices for each batch.
m1 <- SimBatches$batches[[1]]
m2 <- SimBatches$batches[[2]]

# Fix column (cells) names.
colnames(m1) <- paste0("B1_", colnames(m1))
colnames(m2) <- paste0("B2_", colnames(m2))

# Initialize batch information.
b1 <- rep("B1", ncol(m1))
b2 <- rep("B2", ncol(m2))

# Combine batches.
b <- c(b1, b2)
m <- cbind(m1, m2)
cellnames <- colnames(m)

# Create Seurat object.
x <- Seurat::CreateSeuratObject(Seurat::as.sparse(m))
x$batch <- b

# Create SingleCellExperiment object.
y <- SingleCellExperiment::SingleCellExperiment(list(counts=m, logcounts=m))
y$batch <- b

# RunCanek.
x <- RunCanek(x, "batch", slot="counts")
y <- RunCanek(y, "batch")
z <- RunCanek(list(B1=m1, B2=m2), debug = TRUE)

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

x <- RunCanek(x, "batch", integration.name = "CanekRNA")
y <- RunCanek(y, "batch", integration.name = "CanekRNA")

test_that("Setting RunCanek integration.name argument works", {
  expect_true("CanekRNA" %in% names(x))
  expect_true("CanekRNA" %in% names(SummarizedExperiment::assays(y)))
})
