library("adjclust")

context("Check plotSim plots for all types of input")

test_that("'plotSim' works for 'matrix'", {
  sim <- matrix(c(1.0, 0.1, 0.2, 0.3, 0.1, 1.0 ,0.4 ,0.5, 0.2, 0.4, 1.0, 0.6, 
                  0.3, 0.5, 0.6, 1.0), 
                nrow = 4)
  expect_error(plotSim(sim, "similarity", xaxis = TRUE, naxis = 2), NA)
  
  fit1 <- adjClust(sim, "similarity")
  expect_error(plotSim(sim, "similarity", dendro = fit1), NA)
  
  clust1 <- cutree(fit1, k = 2)
  expect_error(
    plotSim(sim, "similarity", dendro = fit1, clustering = clust1),
    NA)
})

test_that("'plotSim' works for 'dgCMatrix'", {
  sim <- Matrix::Matrix(
    c(0, 2:0, 0, 0, 0, 2:0, 0, 0, 0, 2:0, 2:0, 0, 2:0, 0, 0),
    5, 5)
  sim <- sim + t(sim)
  expect_error(plotSim(sim, "dissimilarity", xaxis = TRUE, naxis = 2), NA)
  
  expect_message(
    fit1 <- adjClust(1 - 2*sim, "similarity"),
    "merges with non increasing heights")
  fit1 <- correct(fit1)
  expect_error(plotSim(sim, "dissimilarity", dendro = fit1), NA)
  
  clust1 <- cutree(fit1, k = 2)
  expect_error(
    plotSim(sim, "dissimilarity", dendro = fit1, clustering = clust1),
    NA)
})

test_that("'plotSim' works for 'dsCMatrix'", {
  sim <- Matrix::Matrix(toeplitz(c(10, 0, 1, 0, 3)), sparse = TRUE)
  expect_error(plotSim(sim, "similarity", xaxis = TRUE, naxis = 2), NA)
  
  expect_message(
    fit1 <- adjClust(sim, "similarity"),
    "merges with non increasing heights")
  fit1 <- correct(fit1)
  expect_error(plotSim(sim, "dissimilarity", dendro = fit1), NA)
  
  clust1 <- cutree(fit1, k = 2)
  expect_error(
    plotSim(sim, "dissimilarity", dendro = fit1, clustering = clust1),
    NA)
})

test_that("'plotSim' works for 'dist'", {
  data("iris")
  dissim <- dist(iris[1:10,1:4])^2
  fit0 <- hclust(dissim, method = "ward.D")
  # permute so as to have constrained HAC = HAC
  dissim <- as.dist(as.matrix(dissim)[fit0$order,fit0$order])
  expect_message(
    plotSim(dissim, xaxis = TRUE, naxis = 2),
    "'type' is supposed to be 'dissimilarity'")
  
  sim <- 1-as.matrix(dissim)/2
  fit2 <- adjClust(sim*2/9)
  expect_error(plotSim(dissim, "dissimilarity", dendro = fit2), NA)
})

test_that("'plotSim' works for 'HTCexp'", {
  testthat::skip_if_not_installed("HiTC")
  load(system.file("extdata", "hic_imr90_40_XX.rda", package = "adjclust"))
  expect_error(plotSim(hic_imr90_40_XX, xaxis = TRUE), NA)
})

test_that("'plotSim' works for 'snpStats'", {
  skip_if_not_installed("snpStats")

  data("ld.example", package = "snpStats")
  ceph.1mb[4,286]@.Data[1,1] <- as.raw(3) ## to avoid NaNs
  expect_error(plotSim(ceph.1mb), NA)
  
  expect_error(
    plotSim(ceph.1mb, h = 100, stats = "D.prime", xaxis = TRUE),
    NA)
})

