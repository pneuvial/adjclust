matL <- function(mat, h) {
  p <- ncol(mat)
  x <- rep(1:p, each = h+1)
  y <- x + rep(0:h, p)
  coord <- cbind(x, y)
  out <- rep(0, p * (h+1))
  out[y <= p] <- mat[coord[y <= p, ]]
  out <- matrix(out, ncol = h+1, byrow = TRUE)
  out[ ,-1] <- 2*out[ ,-1]
  return(out)
}

matR <- function(mat, h) {
  if (class(mat) == "dgCMatrix")
    mat <- t(mat)
  p <- ncol(mat)
  x <- rep(p:1, each = h+1)
  y <- x - rep(0:h, p)
  coord <- cbind(x, y)
  out <- rep(0, p * (h+1))
  out[y >= 1] <- mat[coord[y >= 1, ]]
  out <- matrix(out, ncol = h+1, byrow = TRUE)
  out[ ,-1] <- 2*out[ ,-1]
  return(out)
}

checkCondition <- function(mat) {# used only with symmetric matrices to check condition
  tmp <- sweep(-2*mat, 1, diag(mat), "+")
  tmp <- sweep(tmp, 2, diag(mat), "+")
  tmp <- tmp[upper.tri(tmp)]
  if (any(tmp < 0)) {
    return(-min(tmp) * 1.01)
  } else return(NULL)
}
