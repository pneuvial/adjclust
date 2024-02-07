test_that("Methods of class 'chac'", {
  data("iris")
  dissim <- dist(iris[, 1:4])^2
  sim <- 1 - as.matrix(dissim)/2
  fit <- adjClust(sim)
  
  fit2 <- correct(fit)
  expect_error(plot(fit2, mode = "corrected"),
               "Already corrected 'chac' object. 'mode' must be set to 'standard'")
  p <- plot(fit2)
  p <- plot(fit2, nodeLabel = TRUE)
  p <- plot(fit2, nodeLabel = TRUE, leaflab = "none")
  p <- plot(fit2, nodeLabel = TRUE, leaflab = "perpendicular")
  p <- plot(fit2, nodeLabel = TRUE, leaflab = "perpendicular", horiz = TRUE)
  p <- plot(fit2, nodeLabel = TRUE, leaflab = "textlike")
  attr(fit2, "edgetext") <- "test text"  # does not work
  p <- plot(fit2, nodeLabel = TRUE)
  p <- plot(fit2, nodeLabel = TRUE, leaflab = "textlike")
  
  fit_h <- hclust(dissim)
  expect_error(cutree_chac(fit_h), "'tree' must be of class 'chac'")
})