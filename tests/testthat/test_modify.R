library("adjclust")
library("Matrix")

context("Correctness of modify function")

test_that("Output of the modify function satisfies the conditions: For all 1 <= i,j <=p, (s[i,i] == 1) and (s[i,j] <= 0.5*(s[i,i]+s[j,j])) )", {
  
  m <- matrix(c(1,0.1,0.5,0.8,0,0.1,1,0.2,0.6,0.9,0.5,0.2,1,0.3,0.7,0.8,0.6,0.3,1,0.4,0,0.9,0.7,0.4,1), nrow=5)
  p <- nrow(m)
  h <- p-1
  
  #Case 1: Input matrix satisfies both conditions (s[i,i] == 1) and (s[i,j] <= 0.5*(s[i,i]+s[j,j])) for all 1 <= i,j <=p
  m1 <- m
  expect_true(adjclust:::condnCheck(m))
  res1 <- adjclust:::modify(m, p, h)
  expect_identical(res1, m)
  expect_true(adjclust:::condnCheck(res1))
  expect_equal(diag(res1), rep(1, p))  
    
  #Case 2: Input matrix satisfies condition (s[i,i] == 1) but does NOT satisfy condition (s[i,j] <= 0.5*(s[i,i]+s[j,j])) for all 1 <= i,j <=p
  m2 <- m
  m2[2,3] <- 2.5
  
  expect_false(adjclust:::condnCheck(m2))
  res2 <- adjclust:::modify(m2, p, h)
  expect_true(adjclust:::condnCheck(res2))
  expect_equal(diag(res2), rep(1, p))  
    
  #Case 3: Input matrix satisfies condition (s[i,j] <= 0.5*(s[i,i]+s[j,j]))  but does NOT satisfy condition (s[i,i] == 1) for all 1 <= i,j <=p
  m3 <- m
  m3[2,2] <- 2.2
  expect_true(adjclust:::condnCheck(m3))  ## only checks for dominant diagonal
  res3 <- adjclust:::modify(m2, p, h)
  expect_true(adjclust:::condnCheck(res3))
  expect_equal(diag(res3), rep(1, p))

  #Case 4: Input matrix does NOT satisfy any of the conditions (s[i,j] <= 0.5*(s[i,i]+s[j,j])) and (s[i,i] == 1) for all 1 <= i,j <=p
  m4 <- m
  m4[2,3] <- 2.5
  m4[2,2] <- 2.2
  expect_false(adjclust:::condnCheck(m4))
  res4 <- adjclust:::modify(m2, p, h)
  expect_true(adjclust:::condnCheck(res4))
  expect_equal(diag(res4), rep(1, p))
})
