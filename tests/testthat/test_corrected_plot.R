library("adjclust")

context("Check that the corrected plots have increasing heights")

test_that("'adjClust' returns a dendrogram with increasing heights for 
          'mode=corrected'", {
  data("iris")
  dissim <- dist(iris[ ,1:4])^2
  sim <- 1-as.matrix(dissim)/2
  fit <- adjClust(sim)
  corrected_dendo <- as.hclust(plot(fit, mode = "corrected"))
  expect_equal(sum(diff(corrected_dendo$height) < 0), 0)
})
