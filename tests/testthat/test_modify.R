library("adjclust")

context("Correctness of handling general similarity matrices")

data("iris")
dissim <- dist(iris[1:10,1:4])^2
sim <- 1-as.matrix(dissim)/2
sim2 <- sim - diag(0.9, ncol(sim))

fit2 <- adjClust(sim)
fit3 <- adjClust(sim + diag(rep(3, ncol(sim))))

test_that("Results of 'adjclust' are shifted by lambda when similarity is shifted by lambda", {
  expect_equal(fit2$height, fit3$height - 3, tolerance = 0.00001)
  expect_equal(fit2$merge, fit3$merge)
})

test_that("Results of the algorithm are shifted by lambda when similarity is unnormalized and heights are positive", {
  expect_message(fit3 <- adjClust(sim2), "added")

  tmp <- sweep(-2*sim2, 1, diag(sim2), "+")
  tmp <- sweep(tmp, 2, diag(sim2), "+")
  tmp <- tmp[upper.tri(tmp)]
  expect_equal(fit2$height, fit3$height + min(tmp) * 1.01 + 0.9, tolerance = 0.00001)
  expect_equal(fit2$merge, fit3$merge)
  expect_equal(sum(fit3$height < 0), 0)
})