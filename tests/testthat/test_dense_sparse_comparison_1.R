library("adjclust")
library("Matrix")

context("Comparison between the results of adjClust with sparse and dense matrices")

test_that("test that adjClust gives identical results for sparse matrices and dense matrices when h<p-1", {

  m <- matrix(c(1,0,0,0,0,0.1,1,0,0,0,0.5,0.2,1,0,0,0.8,0.6,0.3,1,0,0,0.9,0.7,0.4,1), nrow=5)
  dgc <- as(m, "dgCMatrix")
  dsc <- as(forceSymmetric(m), "dsCMatrix")        
  p <- nrow(m)
  
  m1 <- as.matrix(forceSymmetric(m))
  
  fit1 <- adjClust(m1, "similarity", p-2)
  fit2 <- adjClust(dgc, "similarity", p-2)
  fit3 <- adjClust(dsc, "similarity", p-2)
  
  expect_equal(fit1$merge, fit2$merge)
  expect_equal(fit1$height, fit2$height)
  
  expect_equal(fit1$merge, fit3$merge)
  expect_equal(fit1$height, fit3$height)
  
})
