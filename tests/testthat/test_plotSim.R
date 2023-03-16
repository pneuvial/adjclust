library("adjclust")

context("Check plotSim plots for all types of input")

test_that("'plotSim' works for 'matrix'", {
  sim <- matrix(c(1.0, 0.1, 0.2, 0.3, 0.1, 1.0 ,0.4 ,0.5, 0.2, 0.4, 1.0, 0.6, 
                  0.3, 0.5, 0.6, 1.0), 
                nrow = 4)
  p <- plotSim(sim, "similarity", axis = TRUE, naxis = 2)
  expect_s3_class(p, "ggplot")
  
  fit1 <- adjClust(sim, "similarity")
  p <- plotSim(sim, "similarity", dendro = fit1)
  expect_s3_class(p, "ggplot")
  
  clust1 <- cutree_chac(fit1, k = 2)
  p <- plotSim(sim, "similarity", dendro = fit1, clustering = clust1)
  expect_s3_class(p, "ggplot")
})

test_that("'plotSim' works for 'dgCMatrix'", {
  sim <- Matrix::Matrix(
    c(0, 2:0, 0, 0, 0, 2:0, 0, 0, 0, 2:0, 2:0, 0, 2:0, 0, 0),
    5, 5)
  sim <- sim + t(sim)
  p <- plotSim(sim, "dissimilarity", axis = TRUE, naxis = 2)
  expect_s3_class(p, "ggplot")
  
  expect_message({ fit1 <- adjClust(1 - 2*sim, "similarity") },
                 "merges with non increasing heights", fixed = FALSE)
  fit1 <- correct(fit1)
  p <- plotSim(sim, "dissimilarity", dendro = fit1)
  expect_s3_class(p, "ggplot")
  
  clust1 <- cutree_chac(fit1, k = 2)
  p <- plotSim(sim, "dissimilarity", dendro = fit1, clustering = clust1)
  expect_s3_class(p, "ggplot")
})

test_that("'plotSim' works for 'dsCMatrix'", {
  sim <- Matrix::Matrix(toeplitz(c(10, 0, 1, 0, 3)), sparse = TRUE)
  p <- plotSim(sim, "similarity", axis = TRUE, naxis = 2)
  expect_s3_class(p, "ggplot")
  
  expect_message({ fit1 <- adjClust(sim, "similarity") },
                 "merges with non increasing heights", fixed = FALSE)
  fit1 <- correct(fit1)
  p <- plotSim(sim, "dissimilarity", dendro = fit1)
  expect_s3_class(p, "ggplot")
  
  clust1 <- cutree_chac(fit1, k = 2)
  p <- plotSim(sim, "dissimilarity", dendro = fit1, clustering = clust1)
  expect_s3_class(p, "ggplot")
})

test_that("'plotSim' works for 'dist'", {
  data("iris")
  dissim <- dist(iris[1:10, 1:4])^2
  fit0 <- hclust(dissim, method = "ward.D")
  # permute so as to have constrained HAC = HAC
  dissim <- as.dist(as.matrix(dissim)[fit0$order,fit0$order])
  expect_message(
    { p <- plotSim(dissim, axis = TRUE, naxis = 2) },
    "input class is 'dist' so 'type' is supposed to be 'dissimilarity'",
    fixed = FALSE)
  expect_s3_class(p, "ggplot")
  
  sim <- 1-as.matrix(dissim)/2
  fit2 <- adjClust(sim*2/9)
  p <- plotSim(dissim, "dissimilarity", dendro = fit2)
  expect_s3_class(p, "ggplot")
})

test_that("'plotSim' works for 'HTCexp'", {
  testthat::skip_if_not_installed("HiTC")
  load(system.file("extdata", "hic_imr90_40_XX.rda", package = "adjclust"))
  p <- plotSim(hic_imr90_40_XX, axis = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("'plotSim' works for 'snpStats'", {
  skip_if_not_installed("snpStats")

  data("ld.example", package = "snpStats")
  ceph.1mb[4, 286]@.Data[1, 1] <- as.raw(3) ## to avoid NaNs
  p <- plotSim(ceph.1mb)
  expect_s3_class(p, "ggplot")
  
  p <- plotSim(ceph.1mb, h = 100, stats = "D.prime", axis = TRUE)
  expect_s3_class(p, "ggplot")
})

