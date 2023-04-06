context("Check plotSim plots for all types of input")

test_that("'plotSim' works for 'matrix'", {
  sim <- matrix(c(1.0, 0.1, 0.2, 0.3, 
                  0.1, 1.0 ,0.4 ,0.5, 
                  0.2, 0.4, 1.0, 0.6, 
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

  p <- plotSim(sim, dendro = as.hclust(fit1), clustering = clust1)
  expect_s3_class(p, "ggplot")
  
  expect_error(plotSim(sim, clustering = 1))
  expect_error(plotSim(sim, clustering = c(1, 3, 3, 3)))
  expect_s3_class(p, "ggplot")
  
  expect_error(plotSim(sim, dendro = clust1))
  expect_error(plotSim(sim, dendro = clust1))

  expect_error(plotSim(sim, log = 1), "'log' must be logical")
  expect_error(plotSim(sim, legendName = 0), "'legendName' must be a string")
  expect_error(plotSim(sim, main = 1), "'main' must be a string")
  expect_error(plotSim(sim, axis = 1), "'axis' must be logical")
  expect_error(plotSim(sim, priorCount = "1"), "'priorCount' must be a single non-negative number!")
  expect_error(plotSim(sim, axis = TRUE, naxis = 1.2), "'naxis' must be a single value of type integer!")
  expect_warning(plotSim(sim, axis = TRUE, naxis = 5), "Reducing the number of ticks on x-axis to the number of objects.")
  expect_error(plotSim(sim, axis = TRUE, axistext = "1"), "'axistext' length must be equal to 'naxis'")
  expect_error(plotSim(sim, axis = TRUE, xlab = FALSE), "'xlab' must be a string!")

  p <- plotSim(sim, axis = TRUE, axistext = 1:nrow(sim))
  p <- plotSim(sim, log = FALSE)
  p <- plotSim(sim, main = "main title")
  p <- plotSim(sim, axis = FALSE, dendro = fit1)
  mk <- make_coords(1:4, c(1:3, NA), 1:2) 
})

test_that("'plotSim' works for 'dgCMatrix'", {
  sim <- Matrix::Matrix(
    c(0, 2:0, 0, 0, 0, 2:0, 0, 0, 0, 2:0, 2:0, 0, 2:0, 0, 0),
    5, 5)
  expect_warning(plotSim(sim, "dissimilarity"), 
                         "Input matrix was not symmetric. Plotting only the upper-triangular part of the matrix.")
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
  expect_error(plotSim(hic_imr90_40_XX, type = "dissimilarity"), 
                       "type 'dissimilarity' does not match 'HTCexp' data")
})

test_that("'plotSim' works for 'snpMatrix'", {
  skip_if_not_installed("snpStats")

  data("ld.example", package = "snpStats")
  ceph.1mb[4, 286]@.Data[1, 1] <- as.raw(3) ## to avoid NaNs
  p <- plotSim(ceph.1mb)
  expect_s3_class(p, "ggplot")
  
  p <- plotSim(ceph.1mb, h = 100, stats = "D.prime", axis = TRUE)
  expect_s3_class(p, "ggplot")

  msg <- "'h' should be numeric, larger than 0 and smaller than p."
  expect_error(plotSim(ceph.1mb, h = 0), msg)
  expect_error(plotSim(ceph.1mb, h = NA_character_), msg)
  expect_error(plotSim(ceph.1mb, h = ncol(ceph.1mb)), msg)
})

