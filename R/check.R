modify <- function(m) {

  if (!(is.matrix(m) || ("dist" %in% class(m))))
    stop("Input must be a similarity matrix or a dist object")

  if (is.matrix(m)) {
    if (!(nrow(m) == ncol(m)))
      stop("Input matrix is not a square matrix")
    if (!is.numeric(m))
      stop("Input matrix is not numeric")
    if (!(sum(m == t(m)) == (nrow(m)^2)))
      stop("Input matrix is not symmetric")
  }

  if (any(is.na(m)))
    stop("Missing values in the input")
  
  if ("dist" %in% class(m)) {
  d <- as.matrix(m)
  mat <- 1 - 0.5*(d^2)
  } else {
    if (all(diag(m) == 1))
    {
      mat <- m
    } else {
        if (!(condnCheck(m))) {
          m1 <- sweep(m, 1, diag(m))
          l <- max(m1)  
          if (l<0) l <- 0
          diag(m) <- diag(m) + l + 0.1   #Taking epsilon = 0.1 
        }
       mat <- makeDiagOne(m)
     }
  }
  return(mat)
}

condnCheck <- function(m) {
  res <- .Call("CCondnCheck", m)
  return(res)
}

makeDiagOne <- function(m) {
  res <- .Call("CMakeDiagOne", m)
  return(res)
}