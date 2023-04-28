context("Consistency of the results of 'hicClust' across various input formats")

test_that("'hicClust' gives identical results regardless of data input format", {
  testthat::skip_if_not_installed("HiTC")
  #case1: Input as HiTC::HTCexp object
  load(system.file("extdata", "hic_imr90_40_XX.rda", package = "adjclust"))
  
  fit1 <- hicClust(hic_imr90_40_XX)
  fit1_log <- hicClust(hic_imr90_40_XX, log = TRUE)
  expect_error(hicClust(hic_imr90_40_XX, h = "1"), "h should be numeric")
  
  #case2: Input as Matrix::dsCMatrix contact map
  mat <- HiTC::intdata(hic_imr90_40_XX) 
  
  fit2 <- hicClust(mat)
  fit2_log <- hicClust(mat, log = TRUE)

  V1 <- mat@Dimnames[[1]][mat@i+1]          #loci1names
  V2 <- rep(mat@Dimnames[[2]], diff(mat@p)) #loci2names
  V3 <- mat@x
  
  content <- data.frame(V1, V2, V3)
  fit_df <- hicClust(as.matrix(content)) 
  
  tf <- tempfile(fileext = ".txt")
  write.table(content, tf, sep = "\t", col.names = FALSE, row.names = FALSE)
  fit3 <- hicClust(tf)
  fit3b <- hicClust(tf, sep = "\t")
  expect_identical(fit3[-5], fit3b[-5])  # '5' <=> 'call', which does differ

  fit3_log <- hicClust(tf, log = TRUE)
  
  expect_error(hicClust("non existing file name"), 
               "Input of type 'character' should be a valid file.")
  
  df <- read.table(tf, header = FALSE, sep = "\t")
  fit4 <- hicClust(df)
  fit4_log <- hicClust(df, log = TRUE)
  
  # case 5: Input as matrix
  fit5 <- hicClust(as.matrix(mat))
  fit5_log <- hicClust(as.matrix(mat), log = TRUE)
  
  # case 6: Input as Matrix
  fit6 <- hicClust(Matrix(mat))
  fit6_log <- hicClust(Matrix(mat), log = TRUE)
  
  # case 7: Input as dgCMatrix
  smat <- as(as(as(mat, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  fit7 <- hicClust(smat)
  fit7_log <- hicClust(smat, log = TRUE)
  
  # case 8: Input as dgeMatrix
  smat <- as(as(as(mat, "dMatrix"), "generalMatrix"), "unpackedMatrix")
  fit8 <- hicClust(smat)
  fit8_log <- hicClust(smat, log = TRUE)

  expect_equal(fit1$merge, fit2$merge)
  expect_equal(fit1$height, fit2$height)  
  
  expect_equal(fit1$merge, fit3$merge)
  expect_equal(fit1$height, fit3$height)  
  
  expect_equal(fit1$merge, fit4$merge)
  expect_equal(fit1$height, fit4$height)  

  # test that hicClust methods returns expected 'calls'
  expect_identical(as.list(fit1$call)[[1]], as.symbol("hicClust"))
  expect_identical(as.list(fit2$call)[[1]], as.symbol("hicClust"))
  expect_identical(as.list(fit3$call)[[1]], as.symbol("hicClust"))
  expect_identical(as.list(fit4$call)[[1]], as.symbol("hicClust"))
})
