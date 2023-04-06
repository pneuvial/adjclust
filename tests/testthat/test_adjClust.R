test_that("adjClust methods returns expected 'calls'", {
  sim <- matrix(
    c(1.0, 0.1, 0.2, 0.3,
      0.1, 1.0 ,0.4 ,0.5,
      0.2, 0.4, 1.0, 0.6,
      0.3, 0.5, 0.6, 1.0), nrow = 4)
  
  ## similarity, full width
  fit1 <- adjClust(sim, "similarity")
  lst <- as.list(fit1$call)
  expect_identical(lst[[1]], as.symbol("adjClust"))

  ## similarity, h < p-1
  fit2 <- adjClust(sim, "similarity", h = 2)
  lst <- as.list(fit2$call)
  expect_identical(lst[[1]], as.symbol("adjClust"))
  
  ## dissimilarity
  dist <- as.dist(sqrt(2-(2*sim)))
  
  ## dissimilarity, full width
  fit3 <- adjClust(dist, "dissimilarity")
  lst <- as.list(fit3$call)
  expect_identical(lst[[1]], as.symbol("adjClust"))
  
  ## dissimilarity, h < p-1
  fit4 <- adjClust(dist, "dissimilarity", h = 2)
  lst <- as.list(fit4$call)
  expect_identical(lst[[1]], as.symbol("adjClust"))
})

test_that("adjClust methods properly catches unexpected  'calls'", {
  mat <- matrix(NA_character_)
  expect_error(adjClust(mat), "Input matrix is not numeric")

  mat <- matrix(NA_real_)
  expect_error(adjClust(mat), "Missing values in the input are not allowed")
  
  mat <- matrix(1:2)
  expect_error(adjClust(mat), "Input matrix is not symmetric")

  mat <- matrix(rep(1, 4), 2, 2)
  expect_error(adjClust(mat, h = NA_character_), "Input band width 'h' must be numeric")
  expect_error(adjClust(mat, h = -1), "Input band width 'h' must be non negative")
  expect_error(adjClust(mat, h = 0.1), "Input band width 'h' must be an integer")
  expect_error(adjClust(mat, h = 2), "Input band width 'h' must be strictly less than dimensions of matrix")
  adjClust(mat, strictCheck = FALSE)
  
  # dsyMatrix/dgeMatrix
  mat <- matrix(rep(1, 4), 2, 2)
  smat <- as(as(as(mat, "dMatrix"), "symmetricMatrix"), "unpackedMatrix")
  smat[1, 2] <- 2 # automatic coercion to dgeMatrix
  expect_error(adjClust(smat), "Input matrix is not symmetric")
  
  
  # dgTMatrix
  mat <- matrix(rep(1, 4), 2, 2)
  smat <- as(as(as(mat, "dMatrix"), "symmetricMatrix"), "sparseMatrix")
  expect_error(adjClust(smat, type = "dissimilarity"), 
               "'type' can only be 'similarity' with sparse Matrix inputs")
  
  dmat <- as(as(as(smat, "dMatrix"), "symmetricMatrix"), "TsparseMatrix") 
  expect_error(adjClust(dmat, type = "dissimilarity"), 
               "'type' can only be 'similarity' with sparse Matrix inputs")
  dmat <-  as(mat, "dgTMatrix")
  dmat[1, 2] <- 0
  expect_error(adjClust(dmat), "Input matrix is not symmetric")
})

test_that("'matL' and 'matR' are consistent with C++ versions", {
  sim <- matrix(
    c(1.0, 0.1, 0.2, 0.3,
      0.1, 1.0 ,0.4 ,0.5,
      0.2, 0.4, 1.0, 0.6,
      0.3, 0.5, 0.6, 1.0), nrow = 4)
  ml <- matL(sim, h = 2)
  mr <- matR(sim, h = 2)
  mat <- as(sim, "dgCMatrix")
  expect_identical(matL_full(sim, h = 2), ml)
  expect_identical(matL(mat, h = 2), ml)
  expect_identical(matR_full(sim, h = 2), mr)
  expect_identical(matR(mat, h = 2), mr)
  expect_identical(matR_sparse(mat, h = 2), 
                   as(mr, "sparseMatrix"))
  expect_identical(matL_sparse(mat, h = 2), 
                   as(ml, "sparseMatrix"))
})

test_that("WCSS functions", {
  sim <- matrix(
    c(1.0, 0.1, 0.2, 0.3,
      0.1, 1.0 ,0.4 ,0.5,
      0.2, 0.4, 1.0, 0.6,
      0.3, 0.5, 0.6, 1.0), nrow = 4)
  mat <- as(sim, "dgCMatrix")

  clust <- rep(1, ncol(mat))
  wcss_single(mat, clust)
  
  clust_mat <- rbind(clust, clust)
  WCSS(mat, clust_mat)
})