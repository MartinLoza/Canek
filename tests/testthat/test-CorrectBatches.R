set.seed(0)

Batches <- SimBatches$batches
z <- CorrectBatches(Batches, doCosNorm = TRUE, clusterMethod = "kmeans")

test_that("CorrectBatches works", {
  expect_false(is.null(z))
  expect_equal(length(which(is.na(z))),0)
  expect_equal(dim(z), dim(cbind(Batches$B1,Batches$B2)))
  expect_true(ncol(z) == 1579 & nrow(z) == 500)
  expect_equal(z[1,1], 6.133852, tolerance = 0.0001)
})

test_that("CorrectBatches uses hierarchical selection for correctEmbeddings = TRUE with more than two batches", {
  set.seed(42)
  m1 <- SimBatches$batches[[1]] # 631 cells, real structure, paired with m2
  m2 <- SimBatches$batches[[2]] # 948 cells, largest -> becomes the reference
  colnames(m1) <- paste0("Sim_", colnames(m1))
  colnames(m2) <- paste0("Ref_", colnames(m2))

  # unrelated noise batch, sized between m1 and m2 (bigger than m1, smaller than m2) so a
  # size-only merge order would pick it over m1 despite it having no real relationship to m2
  nGenes <- nrow(m1)
  noise <- matrix(rpois(nGenes * 800, lambda = 3), nrow = nGenes)
  rownames(noise) <- rownames(m1)
  colnames(noise) <- paste0("Noise_", seq_len(800))

  batches <- list(Ref = m2, Sim = m1, Noise = noise)

  resHier <- suppressWarnings(CorrectBatches(batches, correctEmbeddings = TRUE, pcaDim = 15,
                                              hierarchical = TRUE, debug = TRUE))
  resSize <- suppressWarnings(CorrectBatches(batches, correctEmbeddings = TRUE, pcaDim = 15,
                                              hierarchical = FALSE, debug = TRUE))

  # hierarchical picks the batch with more MNN pairs with the reference (Sim), regardless of size
  expect_equal(names(resHier)[1], "Ref/Sim")
  # size-only ordering picks the larger remaining batch (Noise) regardless of similarity
  expect_equal(names(resSize)[1], "Ref/Noise")
})
