band <- function(mat, h) {
  x <- .Call("DiagBand", mat, as.integer(h))
  return(x)
}