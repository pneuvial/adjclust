library("adjclust")

context("Comparison between the results of the 'rioja' and 'adjclust' packages")

check_rioja <- function() {
  if (!require("rioja")) {
    skip("rioja package not available")
  }
}

test_that("rioja and adjClustBand_heap with full band give idenctical results on simmatrix.rda", {
  check_rioja()
  
  data("simmatrix", package="adjclust")
  sim <- simmatrix
  rownames(sim) <- colnames(sim) #to avoid dissimilarity error  
  p <- nrow(sim)
      
  dis_sq <- 2 - (2*sim)
  dis <- sqrt(dis_sq)
  
  dis_sq <- as.dist(dis_sq)
  dis <- as.dist(dis)
  
  fit1 <- adjclust:::adjClustBand_heap(sim, "similarity", p-1)
  fit2 <- chclust(dis_sq, method="coniss")
  fit3 <- adjclust:::adjClustBand_heap(dis, "dissimilarity", p-1)
  
  expect_equal(fit1$merge, fit2$merge)
  expect_equal(fit1$height, fit2$height, tolerance=0.000001)
  
  expect_equal(fit3$merge, fit2$merge)
  expect_equal(fit3$height, fit2$height, tolerance=0.000001)
})
