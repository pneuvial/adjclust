library("adjclust")

context("Check that the sum of heights is the dataset (pseudo) inertia")

test_that("'adjClust' returns an object for which the sum of heights is the 
          dataset (pseudo) inertia", {
  data("iris")
  dissim <- dist(iris[ ,1:4])^2
  sim <- 1-as.matrix(dissim)/2
  fit <- adjClust(sim)
  
  expect_equal(sum(fit$height), sum(dissim)/nrow(sim), tolerance = 0.00001)
  
  inertia_dendo <- as.hclust(plot(fit, mode = "within-disp"))
  expect_equal(inertia_dendo$height[nrow(sim)-1], sum(dissim)/nrow(sim), 
               tolerance = 0.00001)
  
  dispersion_dendo <- as.hclust(plot(fit, mode = "total-disp"))
  expect_equal(dispersion_dendo$height[nrow(sim)-1], sum(dissim)/nrow(sim), 
               tolerance = 0.00001)
})
