test_that("Methods of class 'chac'", {
  data("iris")
  dissim <- dist(iris[, 1:4])^2
  sim <- 1-as.matrix(dissim)/2
  fit <- adjClust(sim)
  adjclust::select(fit, type = "capushe")
  adjclust::select(fit, type = "bs", graph = TRUE)
  #adjclust::select(fit, graph = TRUE) # error plot DDSE (uses base "plot" ?)
  
  class(fit)
  print(fit)
  head(fit)
  summary(fit)
  fit2 <- correct(fit)
  expect_error(plot(fit2, mode = "corrected"),
               "Already corrected 'chac' object. 'mode' must be set to 'standard'")
  p <- plot(fit2)
  p <- plot(fit2, nodeLabel = TRUE)
  p <- plot(fit2, leaflab = "none")
  p <- plot(fit2, leaflab = "perpendicular")
  diagnose(fit)
  diagnose(fit, graph = TRUE)
  diagnose(fit, verbose = TRUE)
  
  fit_h <- hclust(dissim)
  expect_error(cutree_chac(fit_h), 
               "'tree' must be of class 'chac'")
})