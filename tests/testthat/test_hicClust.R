library("adjclust")
library("HiTC")

context("Checking the consistency of the results of hicclust function across various input formats")

test_that("hicClust function identical results for the same data in all three inputs formats", {

  #case1: Input as HiTC::HTCexp object
  data("hic_imr90_40$chrXchrX", package="adjclust")
  h <- 3881
  fit1 <- hicClust(obj, h)
  
  #case2: Input as Matrix::dsCMatrix contact map
  mat <- intdata(obj) 
  fit2 <- hicClust(mat, h)
  
  #case3: Input as text file
  V3 <- mat@x
  V1 <- mat@Dimnames[[1]][mat@i+1]          #loci1names
  V2 <- rep(mat@Dimnames[[2]], diff(mat@p)) #loci2names
  
  l1 <- unique(c(as.numeric(V1), as.numeric(V2)))  
  
  l2 <- setdiff(c(72037:75918), l1) 
  V1 <- append(V1, l2)
  V2 <- append(V2, l2)
  V3 <- append(V3, rep(0, length(l2)))
  
  content <- cbind(as.numeric(V1), as.numeric(V2), as.numeric(V3))
  
  tf <- tempfile(fileext = ".txt")
  write.table(content, tf, sep = " ", col.names = FALSE, row.names = FALSE)
  fit3 <- hicClust(tf, h, sep = " ")  
  
  expect_equal(fit1$merge, fit2$merge)
  expect_equal(fit1$height, fit2$height)  
  
  expect_equal(fit2$merge, fit3$merge)
  expect_equal(fit2$height, fit3$height)  
  
})
