library("adjclust")

context("Equivalence between similarity and dissimilarity implementations")

data("iris")
dissim <- as.matrix(dist(iris[1:10,1:4]))
sim <- 12-dissim^2/2
fit1 <- adjClust(sim)

test_that("Case of a dissimilarity of type 'matrix'", {
  fit2 <- adjClust(dissim, type = "dissimilarity")
  
  expect_equal(fit1$height, fit2$height, tolerance = 0.00001)
  expect_equal(fit1$merge, fit2$merge)
})

test_that("Case of a dissimilarity of type 'dist'", {
  dissim <- dist(iris[1:10,1:4])
  expect_message(fit2 <- adjClust(dissim), "type")
  
  expect_equal(fit1$height, fit2$height, tolerance = 0.00001)
  expect_equal(fit1$merge, fit2$merge)
})