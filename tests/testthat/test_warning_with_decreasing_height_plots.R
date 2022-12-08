context("Check that the messages or warnings are produced for decreasing 
        heights")

test_that("'adjClust' returns a note when decreasing heights are produced and
          warnings when such results are plotted with 'mode=standard' and
          'mode=average-disp'", {
  data("iris")
  dissim <- dist(iris[ ,1:4])^2
  sim <- 1-as.matrix(dissim)/2
  expect_message(fit <- adjClust(sim), "merges with non increasing heights")
  
  expect_warning(plot(fit))
  expect_warning(plot(fit, mode = "average-disp"))
})
