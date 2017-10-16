library("adjclust")
library("HiTC")

context("Consistency of the results of 'hicClust' across various input formats")

test_that("'hicClust' gives identical results regardless of data input format", {

  #case1: Input as HiTC::HTCexp object
  data("hic_imr90_40_XX", package="adjclust")
  fit1 <- hicClust(hic_imr90_40_XX)
  
  #case2: Input as Matrix::dsCMatrix contact map
  mat <- intdata(hic_imr90_40_XX) 
  fit2 <- hicClust(mat)
  
  V1 <- mat@Dimnames[[1]][mat@i+1]          #loci1names
  V2 <- rep(mat@Dimnames[[2]], diff(mat@p)) #loci2names
  V3 <- mat@x
  
  content <- cbind(as.numeric(V1), as.numeric(V2), as.numeric(V3))
  
  tf <- tempfile(fileext = ".txt")
  write.table(content, tf, sep = " ", col.names = FALSE, row.names = FALSE)
  fit3 <- hicClust(tf, sep = " ")  
  
  expect_equal(fit1$merge, fit2$merge)
  expect_equal(fit1$height, fit2$height)  
  
  expect_equal(fit2$merge, fit3$merge)
  expect_equal(fit2$height, fit3$height)  
})
