context("Test outputs of diagnose and correct.")

test_that("'diagnose' and 'correct' must return a warning or a message when no reversals are found.", {
  #Sys.setenv("OMP_THREAD_LIMIT" = 2)
  data("iris")
  dissim <- dist(iris[ ,1:4])^2
  sim <- 1-as.matrix(dissim)/2
  fit <- adjClust(sim)
  fit2 <- correct(fit)
  
  expect_equal(sum(diff(fit2$height) < 0), 0)
  expect_warning(correct(fit2))
  expect_message(diagnose(fit2))
  expect_null(suppressWarnings(correct(fit2)))
  expect_null(diagnose(fit2))
  
  corrected_dendro <- as.hclust(plot(fit, mode = "corrected"))
  expect_equal(corrected_dendro$height, fit2$height)
})