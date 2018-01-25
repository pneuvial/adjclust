library("adjclust")
library("Matrix")

context("Correctness of handling general similarity matrices")

test_that("Output of 'checkCondition' is correct", {
  data("iris")
  dissim <- dist(iris[1:10,1:4])^2
  sim <- 1-as.matrix(dissim)/2
  sim2 <- sim - diag(0.9, ncol(sim))
  # normalized similarity
  expect_null(adjclust::checkCondition(sim))
  # unormalized similarity
  expect_gt(adjclust::checkCondition(sim2), 0)
})

test_that("Results of 'adjclust' are shifted by lambda when similarity is shifted by lambda", {
  fit2 <- adjClust(sim)
  fit3 <- adjClust(sim + diag(rep(3, ncol(sim))))
  
  expect_equal(fit2$height, fit3$height - 3, tolerance = 0.00001)
  expect_equal(fit2$merge, fit3$merge)
})

test_that("Expect a message when an unnormalized matrix is shifted", {
  expect_message(fit3 <- adjClust(sim2), "added")
})

test_that("Results of the algorithm are shifted by lambda when similarity is unnormalized and heights are positive", {  
  expect_equal(fit2$height, fit3$height - adjclust::checkCondition(sim2) + 0.9, tolerance = 0.00001)
  expect_equal(fit2$merge, fit3$merge)
  expect_equal(sum(fit3$height < 0), 0)
})