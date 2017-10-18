library("adjclust")
library("Matrix")

context("Correctness of modifySparse function")

test_that("Output of the modifySparse function satisfies the conditions :For all 1 <= i,j <=p, (s[i,i] == 1) and (s[i,j] <= 0.5*(s[i,i]+s[j,j])) )", {
  
  m <- matrix(c(1,0.1,0.5,0.8,0,0.1,1,0.2,0.6,0.9,0.5,0.2,1,0.3,0.7,0.8,0.6,0.3,1,0.4,0,0.9,0.7,0.4,1), nrow=5)
  sp <- as(m, "dgCMatrix")
  p <- nrow(m)
  h <- p-1
  sb <- adjclust:::sparseBand(sp@x, sp@p, sp@i, p, h)

  
  #Case 1 :Input matrix satisfies both conditions (s[i,i] == 1) and (s[i,j] <= 0.5*(s[i,i]+s[j,j])) for all 1 <= i,j <=p
  sb1 <- sb
  res1 <- adjclust:::modifySparse(sb, p, h)
  
  expect_equal(adjclust:::sparseCondnCheck(res1, p, h), TRUE)
  expect_equal(res1[ ,1], rep(1, p))  
  
  #Case 2 :Input matrix satisfies condition (s[i,i] == 1) and but does NOT satisfy condition (s[i,j] <= 0.5*(s[i,i]+s[j,j])) for all 1 <= i,j <=p
  sb2 <- sb
  sb2[3,2] <- 2.5
  res2 <- adjclust:::modifySparse(sb2, p, h)
  
  expect_equal(adjclust:::sparseCondnCheck(res2, p, h), TRUE)
  expect_equal(res2[ ,1], rep(1, p))  
  
  #Case 3 :Input matrix satisfies condition (s[i,j] <= 0.5*(s[i,i]+s[j,j])) and but does NOT satisfy condition (s[i,i] == 1) for all 1 <= i,j <=p
  sb3 <- sb
  sb3[2,1] <- 2.2
  res3 <- adjclust:::modifySparse(sb3, p, h)

  expect_equal(adjclust:::sparseCondnCheck(res3, p, h), TRUE)
  expect_equal(res3[ ,1], rep(1, p))   
  
  #Case 4 :Input matrix does NOT satisfies both conditions (s[i,j] <= 0.5*(s[i,i]+s[j,j])) and (s[i,i] == 1) for all 1 <= i,j <=p
  sb4 <- sb
  sb4[3,2] <- 2.5
  sb4[4,1] <- 2.2
  res4 <- adjclust:::modifySparse(sb4, p, h)
  
  expect_equal(adjclust:::sparseCondnCheck(res4, p, h), TRUE)
  expect_equal(res3[ ,1], rep(1, p))
  
})