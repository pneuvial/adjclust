library("adjclust")
library("Matrix")

context("Comparison between the results of sparseBand function with dgCMatrix and corresponding dsCMatrix")

test_that("sparseBand function gives identical results for dgCMatrix and corresponding dsCMatrix", {
  
  m <- matrix(c(1,0,0,0,0,0.1,1,0,0,0,0.5,0.2,1,0,0,0.8,0.6,0.3,1,0,0,0.9,0.7,0.4,1), nrow=5)
  dgc <- as(m, "dgCMatrix")
  dsc <- as(forceSymmetric(m), "dsCMatrix")        
  p <- nrow(m)
  h <- 3
  
  res1 <- adjclust:::sparseBand(dgc@x, dgc@p, dgc@i, p, h)
  res2 <- adjclust:::sparseBand(dsc@x, dsc@p, dsc@i, p, h)
  
  expect_equal(res1, res2)

})