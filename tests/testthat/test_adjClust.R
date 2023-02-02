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

  mat <- matrix(1:2)
  expect_error(adjClust(mat), "Input matrix is not symmetric")

  # dsyMatrix
  mat <- matrix(rep(1, 4), 2, 2)
  smat <- as(as(as(mat, "dMatrix"), "symmetricMatrix"), "unpackedMatrix")
  smat[1, 2] <- 2
  expect_error(adjClust(smat), "Input matrix is not symmetric")
  
  # dgTMatrix
  smat <- as(as(as(mat, "dMatrix"), "symmetricMatrix"), "sparseMatrix")
  expect_error(adjClust(smat, type = "dissimilarity"), 
               "'type' can only be 'similarity' with sparse Matrix inputs")
  
  dmat <-  as(mat, "dsTMatrix")
  expect_error(adjClust(dmat, type = "dissimilarity"), 
               "'type' can only be 'similarity' with sparse Matrix inputs")
  dmat <-  as(mat, "dgTMatrix")
  dmat[1, 2] <- 0
  expect_error(adjClust(dmat), "Input matrix is not symmetric")
})