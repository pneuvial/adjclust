context("Test cuttree in various situations (decreasing merges or not, k and/or
        h given.")

test_that("'cuttree_chac' must ignore 'h' when reversals are present.", {
  data("iris")
  dissim <- dist(iris[ ,1:4])^2
  sim <- 1 - as.matrix(dissim)/2
  fit <- adjClust(sim)
  fit2 <- correct(fit)
  
  clust1 <- cutree_chac(fit, k = 4)
  expect_error(cutree_chac(fit, h = 0.02))
  clust2 <- cutree_chac(fit, k = 4, h = 50)
  expect_equal(clust1, clust2)
  
  clust3 <- cutree_chac(fit2, k = 4)
  expect_equal(clust1, clust3)
  height4four <- (fit2$height[nrow(sim) - 3] + fit2$height[nrow(sim) - 4])/2
  clust4 <- cutree_chac(fit2, h = height4four)
  expect_equal(clust1, clust4)
  clust5 <- cutree_chac(fit2, k = 4, h = 50)
  expect_equal(clust1, clust5)
})
