library("adjclust")

context("Correctness of the band.R function")

test_that("band.R function is correctly extracting data around diagonal on Rsquared linkage disequilibrium data", {
  
  data("R2.100", package="adjclust")
  h <- 100
  
  m1 <- as.matrix(R2.100)
  m1 <- m1 + t(m1)  ## upper triangular matrix to symmetric matrix
  
  x1 <- adjclust:::band(r,h)
    
  expect_equal(x1, R2.100@x)

})