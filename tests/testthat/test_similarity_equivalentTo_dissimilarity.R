context("Equivalence between similarity and dissimilarity implementations")

#Sys.setenv("OMP_THREAD_LIMIT" = 2)
data("iris")
dissim <- as.matrix(dist(iris[1:10,1:4]))
sim <- 12-dissim^2/2
fit1 <- adjClust(sim)

test_that("Case of a dissimilarity of type 'matrix'", {
  #Sys.setenv("OMP_THREAD_LIMIT" = 2)
  fit2 <- adjClust(dissim, type = "dissimilarity")
  
  expect_equal(fit1$height, fit2$height, tolerance = 0.00001)
  expect_equal(fit1$merge, fit2$merge)
})

test_that("Case of a dissimilarity of type 'dist'", {
  #Sys.setenv("OMP_THREAD_LIMIT" = 2)
  dissim <- dist(iris[1:10,1:4])
  expect_message(fit2 <- adjClust(dissim), "type")
  
  expect_equal(fit1$height, fit2$height, tolerance = 0.00001)
  expect_equal(fit1$merge, fit2$merge)
})