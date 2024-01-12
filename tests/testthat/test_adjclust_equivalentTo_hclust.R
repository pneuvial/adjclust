context("Comparison between the results of the 'hclust' and 'adjclust' when 
        optimal merge is always an adjacent merge")

test_that("'hclust' and 'adjClust' give identical results on toy data when the
          best merges are always adjacent merges", {
  #Sys.setenv("OMP_THREAD_LIMIT" = 2)
  data("iris")
  dissim <- dist(iris[1:10,1:4])^2  ## Note the "^2"
  fit0 <- hclust(dissim, method = "ward.D")
  # permute so as to have constrained HAC = HAC
  dissim <- as.dist(as.matrix(dissim)[fit0$order,fit0$order])
  fit1 <- hclust(dissim/9, method = "ward.D")
  
  sim <- 1-as.matrix(dissim)/2
  fit2 <- adjClust(sim*2/9)
  
  expect_equal(fit1$height, fit2$height, tolerance = 0.00001)
  expect_equal(fit1$merge, fit2$merge)

  ## simpler and equivalent:
  fit1 <- hclust(dissim, method = "ward.D")
  
  sim <- 2-as.matrix(dissim)  ## why are there 2:s on the diagonal?
  fit2 <- adjClust(sim)
  
  expect_equal(fit1$height, fit2$height, tolerance = 0.00001)
  expect_equal(fit1$merge, fit2$merge)
})
