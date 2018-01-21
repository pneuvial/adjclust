library("adjclust")

context("Comparison between the results of the 'rioja' and 'adjclust' packages")

check_rioja <- function() {
  if (!require("rioja")) {
    skip("rioja package not available")
  }
}

test_that("rioja and adjClust with full band give identical results on toy data", {
  check_rioja()
  
  data("iris")
  sim <- cor(t(iris[, 1:4]))
  p <- nrow(sim)
      
  dis_sq <- 2 - (2*sim)
  dis <- sqrt(dis_sq)
  
  dis_sq <- as.dist(dis_sq)
  dis <- as.dist(dis)
  
  fit1 <- adjClust(sim, "similarity", p-1)
  fit2 <- chclust(dis_sq, method="coniss")
  fit3 <- adjClust(dis, "dissimilarity", p-1)
  fit4 <- adjClust(sim, "similarity")
  
  expect_equal(fit1$merge, fit2$merge)
  expect_equal(cumsum(fit1$height), fit2$height, tolerance = 0.00001)
  
  expect_equal(fit3$merge, fit2$merge)
  expect_equal(cumsum(fit3$height), fit2$height, tolerance = 0.00001)
  
  expect_equal(fit4$merge, fit2$merge)
  expect_equal(cumsum(fit4$height), fit2$height, tolerance = 0.00001)
})
