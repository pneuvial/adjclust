## x : the LD matrix
.toMatLeft <- function(x, p, h) {
  ## x <- slot(t(x), "x")
  x1 <- utils::head(x, (p-1)*h-h*(h-1))
  m1 <- matrix(x1, ncol=h, byrow=TRUE)
  x2 <- utils::tail(x, h*(h-1)/2)
  m2t <- matrix(nrow=h, ncol=h-1)
  idxs <- seq(along=x2)
  for (jj in 1:(h-1)) {
    idxsJJ <- utils::head(idxs, h-jj)
    idxs <- utils::tail(idxs, -(h-jj))
    m2t[1:(h-jj), jj] <- x2[idxsJJ]
  }
  mat <- rbind(m1, t(m2t), NA)
  mat[is.na(mat)] <- 0
  return(mat)
}

## x : the LD matrix
.toMatRight <- function(x, p, h){
  ## x <- slot(x, "x")
  x1 <- utils::tail(x, (p-1)*h-h*(h-1))
  m1 <- matrix(x1, ncol=h, byrow=TRUE)
  x2 <- utils::head(x, h*(h-1)/2)
  m2t <- matrix(nrow=h-1, ncol=h)
  idxs <- seq(along=x2)
  for (jj in (h-1):1) {
    idxsJJ <- utils::tail(idxs, jj)
    idxs <- utils::head(idxs, -jj)
    m2t[jj, (h-jj+1):h] <- x2[idxsJJ]
  }
  mat <- rbind(NA, m2t, m1)
  mat[is.na(mat)] <- 0
  return(mat)
}

.rotate <- function(mat) { t(mat[nrow(mat):1,,drop=FALSE]) }
