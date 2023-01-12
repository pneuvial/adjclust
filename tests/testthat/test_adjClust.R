test_that("adjClust methods returns expected 'calls'", {
  sim <- matrix(
    c(1.0, 0.1, 0.2, 0.3,
      0.1, 1.0 ,0.4 ,0.5,
      0.2, 0.4, 1.0, 0.6,
      0.3, 0.5, 0.6, 1.0), nrow = 4)
  
  ## similarity, full width
  fit1 <- adjClust(sim, "similarity")
  lst <- as.list(fit1$call)
  expect_identical(lst[[1]], as.symbol("adjClust"))

  ## similarity, h < p-1
  fit2 <- adjClust(sim, "similarity", h = 2)
  lst <- as.list(fit2$call)
  expect_identical(lst[[1]], as.symbol("adjClust"))
  
  ## dissimilarity
  dist <- as.dist(sqrt(2-(2*sim)))
  
  ## dissimilarity, full width
  fit3 <- adjClust(dist, "dissimilarity")
  lst <- as.list(fit3$call)
  expect_identical(lst[[1]], as.symbol("adjClust"))
  
  ## dissimilarity, h < p-1
  fit4 <- adjClust(dist, "dissimilarity", h = 2)
  lst <- as.list(fit4$call)
  expect_identical(lst[[1]], as.symbol("adjClust"))
})
