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
x <- RunCanek(x, "batch", slot="counts", correctEmbeddings = FALSE)
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

x <- RunCanek(x, "batch", integration.name = "CanekRNA", correctEmbeddings = FALSE)
y <- RunCanek(y, "batch", integration.name = "CanekRNA")

test_that("Setting RunCanek integration.name argument works", {
  expect_true("CanekRNA" %in% names(x))
  expect_true("CanekRNA" %in% names(SummarizedExperiment::assays(y)))
})

test_that("RunCanek defaults to correctEmbeddings and falls back to pcaDim=30 with a warning when no PCA reduction exists", {
  xe <- Seurat::CreateSeuratObject(Seurat::as.sparse(m))
  xe$batch <- b

  expect_warning(
    xe <- RunCanek(xe, "batch", slot = "counts"),
    "No existing 'pca' reduction"
  )

  expect_true("canek" %in% Seurat::Reductions(xe))
  expect_equal(Seurat::Assays(xe), "RNA")
  expect_equal(ncol(Seurat::Embeddings(xe, "canek")), 30)
})

test_that("RunCanek infers pcaDim from an existing PCA reduction", {
  xe <- Seurat::CreateSeuratObject(Seurat::as.sparse(m))
  xe$batch <- b
  xe <- Seurat::NormalizeData(xe, verbose = FALSE)
  xe <- Seurat::FindVariableFeatures(xe, nfeatures = 100, verbose = FALSE)
  xe <- Seurat::ScaleData(xe, verbose = FALSE)
  xe <- Seurat::RunPCA(xe, npcs = 15, verbose = FALSE)

  expect_warning(
    xe <- RunCanek(xe, "batch"),
    NA
  )

  expect_equal(ncol(Seurat::Embeddings(xe, "canek")), 15)
})

test_that("RunCanek respects an explicitly passed pcaDim over inference and the default", {
  xe <- Seurat::CreateSeuratObject(Seurat::as.sparse(m))
  xe$batch <- b

  expect_warning(
    xe <- RunCanek(xe, "batch", slot = "counts", pcaDim = 5),
    NA
  )

  expect_equal(ncol(Seurat::Embeddings(xe, "canek")), 5)
})

test_that("CorrectBatches sets maxLoop to 1 with a warning when correctEmbeddings = FALSE", {
  expect_warning(
    z <- CorrectBatches(list(B1 = m1, B2 = m2), maxLoop = 5),
    "maxLoop > 1 is only supported for correctEmbeddings = TRUE"
  )
  expect_equal(dim(z), dim(m))
})

test_that("CorrectBatches iterates the correction and records the correction magnitude when correctEmbeddings = TRUE", {
  z <- CorrectBatches(list(B1 = m1, B2 = m2), correctEmbeddings = TRUE, pcaDim = 10,
                       maxLoop = 4, loopTol = 1e-8, debug = TRUE)
  info <- z$`B2/B1`$debug$info

  expect_equal(info$loops, 4)
  expect_length(info$loopMagnitude, 4)
  # each pass should keep refining the correction, i.e. shrinking its magnitude
  expect_true(all(diff(info$loopMagnitude) < 0))
})

test_that("CorrectBatches stops iterating early once the correction magnitude stops improving", {
  z <- CorrectBatches(list(B1 = m1, B2 = m2), correctEmbeddings = TRUE, pcaDim = 10,
                       maxLoop = 10, loopTol = 0.4, debug = TRUE)
  info <- z$`B2/B1`$debug$info

  expect_true(info$loops < 10)
  expect_length(info$loopMagnitude, info$loops)
})

test_that("RunCanek forwards maxLoop/loopTol to CorrectBatches through ...", {
  xe <- Seurat::CreateSeuratObject(Seurat::as.sparse(m))
  xe$batch <- b
  xe <- Seurat::NormalizeData(xe, verbose = FALSE)
  xe <- Seurat::FindVariableFeatures(xe, nfeatures = 100, verbose = FALSE)
  xe <- Seurat::ScaleData(xe, verbose = FALSE)
  xe <- Seurat::RunPCA(xe, npcs = 10, verbose = FALSE)

  xe <- RunCanek(xe, "batch", maxLoop = 3, loopTol = 1e-8, debug = TRUE)
  info <- xe@tools$RunCanek[[1]]$debug$info

  expect_equal(info$loops, 3)
})
