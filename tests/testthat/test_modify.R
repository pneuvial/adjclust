library("adjclust")
library("Matrix")

context("Correctness of modify function")

test_that("Output of the modify function satisfies the conditions :For all 1 <= i,j <=p, (s[i,i] == 1) and (s[i,j] <= 0.5*(s[i,i]+s[j,j])) )", {
  
  m <- matrix(c(1,0.1,0.5,0.8,0,0.1,1,0.2,0.6,0.9,0.5,0.2,1,0.3,0.7,0.8,0.6,0.3,1,0.4,0,0.9,0.7,0.4,1), nrow=5)
  p <- nrow(m)
  h <- p-1
  
  #Case 1 :Input matrix satisfies both conditions (s[i,i] == 1) and (s[i,j] <= 0.5*(s[i,i]+s[j,j])) for all 1 <= i,j <=p
  m1 <- m
  res1 <- adjclust:::modify(m, p, h)
  
  expect_equal(adjclust:::condnCheck(res1), TRUE)
  expect_equal(diag(res1), rep(1, p))  
    
  #Case 2 :Input matrix satisfies condition (s[i,i] == 1) and but does NOT satisfy condition (s[i,j] <= 0.5*(s[i,i]+s[j,j])) for all 1 <= i,j <=p
  m2 <- m
  m2[2,3] <- 2.5
  res2 <- adjclust:::modify(m2, p, h)

  expect_equal(adjclust:::condnCheck(res2), TRUE)
  expect_equal(diag(res2), rep(1, p))  
    
  #Case 3 :Input matrix satisfies condition (s[i,j] <= 0.5*(s[i,i]+s[j,j])) and but does NOT satisfy condition (s[i,i] == 1) for all 1 <= i,j <=p
  m3 <- m
  m3[2,2] <- 2.2
  res3 <- adjclust:::modify(m2, p, h)
  
  expect_equal(adjclust:::condnCheck(res3), TRUE)
  expect_equal(diag(res3), rep(1, p))

  #Case 4 :Input matrix does NOT satisfies both conditions (s[i,j] <= 0.5*(s[i,i]+s[j,j])) and (s[i,i] == 1) for all 1 <= i,j <=p
  m4 <- m
  m4[2,3] <- 2.5
  m4[2,2] <- 2.2
  res4 <- adjclust:::modify(m2, p, h)
  
  expect_equal(adjclust:::condnCheck(res4), TRUE)
  expect_equal(diag(res4), rep(1, p))
  
})