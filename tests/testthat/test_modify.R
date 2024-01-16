context("Correctness of handling general similarity matrices")

#Sys.setenv("OMP_THREAD_LIMIT" = 2)
data("iris")
dissim <- dist(iris[1:10,1:4])^2
sim <- 1-as.matrix(dissim)/2
sim2 <- sim - diag(0.9, ncol(sim))

fit <- adjClust(sim)
fit2 <- adjClust(sim + diag(rep(3, ncol(sim))))

test_that("Results of 'adjclust' are shifted by lambda when similarity is shifted by lambda", {
  #Sys.setenv("OMP_THREAD_LIMIT" = 2)
  expect_equal(fit$height, fit2$height - 3, tolerance = 0.00001)
  expect_equal(fit$merge, fit2$merge)
  expect_equal(fit$correction, 0)
})

test_that("Results of the algorithm are shifted by lambda when similarity is unnormalized and heights are positive", {
  #Sys.setenv("OMP_THREAD_LIMIT" = 2)
  expect_message(fit3 <- adjClust(sim2), "added")
  expect_message(fit4 <- adjClust(sim2), fit3$correction)

  tmp <- sweep(-2*sim2, 1, diag(sim2), "+")
  tmp <- sweep(tmp, 2, diag(sim2), "+")
  tmp <- tmp[upper.tri(tmp)]
  expect_equal(fit$height, fit3$height + min(tmp) * 1.01 + 0.9, tolerance = 0.00001)
  expect_equal(fit$merge, fit3$merge)
  expect_equal(sum(fit3$height < 0), 0)
  expect_gt(fit3$correction, 0)
})

test_that("A message is displayed when 'select' is used on results obtained from preprocessed matrices", {
  #Sys.setenv("OMP_THREAD_LIMIT" = 2)
  suppressMessages({fit3 <- adjClust(sim2)})
  expect_message(adjclust::select(fit3, type = "bstick"), "might be spurious")
})