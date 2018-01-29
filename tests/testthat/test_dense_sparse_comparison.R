library("adjclust")

context("Comparison between the results of adjClust with sparse and dense matrices")

mat <- matrix(c(1,0,0,0,0,0.1,1,0,0,0,0.5,0.2,1,0,0,0.8,0.6,0.3,1,0,0,0.9,0.7,0.4,1), nrow = 5)
smat1 <- as(mat, "dgCMatrix")
mat <- Matrix::forceSymmetric(mat)
smat2 <- as(mat, "dsCMatrix")
mat <- as(mat, "matrix")
p <- nrow(mat)

test_that("test that adjClust gives identical results for sparse and dense matrices when h < p-1", {
  fit1 <- adjClust(mat, h = 2)
  fit2 <- adjClust(smat1, h = 2)
  fit3 <- adjClust(smat2, h = 2)
  
  expect_equal(fit1$merge, fit2$merge)
  expect_equal(fit1$height, fit2$height)
  
  expect_equal(fit1$merge, fit3$merge)
  expect_equal(fit1$height, fit3$height)
})

test_that("test that adjClust gives identical results for sparse and dense matrices when h is p-1", {
  fit1 <- adjClust(mat)
  fit2 <- adjClust(smat1)
  fit3 <- adjClust(smat2)
  
  expect_equal(fit1$merge, fit2$merge)
  expect_equal(fit1$height, fit2$height)
  
  expect_equal(fit1$merge, fit3$merge)
  expect_equal(fit1$height, fit3$height)
})