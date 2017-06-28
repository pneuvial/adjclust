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
  p <- nrow(sim)
      
  dist <- 2 - (2*sim)
  dist <- as.dist(dist)
  
  fit1 <- adjclust:::adjClustBand_heap(sim, p-1)
  fit2 <- chclust(dist,method="coniss")
  fit3 <- adjclust:::adjClustBand_heap(dist, p-1)
  
  expect_equal(fit1$merge, fit2$merge)
  expect_equal(fit3$merge, fit2$merge)
  
})