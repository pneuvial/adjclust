library("adjclust")

context("Consistency of the results of 'hicClust' across various input formats")

test_that("'hicClust' gives identical results regardless of data input format", {
  testthat::skip_if_not_installed("HiTC")
  #case1: Input as HiTC::HTCexp object
  load(system.file("extdata", "hic_imr90_40_XX.rda", package = "adjclust"))
  
  fit1 <- hicClust(hic_imr90_40_XX)

  #case2: Input as Matrix::dsCMatrix contact map
  mat <- HiTC::intdata(hic_imr90_40_XX) 
  
  fit2 <- hicClust(mat)
  
  V1 <- mat@Dimnames[[1]][mat@i+1]          #loci1names
  V2 <- rep(mat@Dimnames[[2]], diff(mat@p)) #loci2names
  V3 <- mat@x
  
  content <- data.frame(V1, V2, V3)
  
  tf <- tempfile(fileext = ".txt")
  write.table(content, tf, sep = " ", col.names = FALSE, row.names = FALSE)
  fit3 <- hicClust(tf, sep = " ")
  df <- read.table(tf, header = FALSE, sep = " ")
  fit4 <- hicClust(df)
  
  expect_equal(fit1$merge, fit2$merge)
  expect_equal(fit1$height, fit2$height)  
  
  expect_equal(fit1$merge, fit3$merge)
  expect_equal(fit1$height, fit3$height)  
  
  expect_equal(fit1$merge, fit4$merge)
  expect_equal(fit1$height, fit4$height)  
})
