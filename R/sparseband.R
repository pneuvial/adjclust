sparseBand <- function(xslot, pslot, islot, p, h) {
  res <- .Call("CSparseBand", xslot, pslot, islot, as.integer(p), as.integer(h))
  return(res)  
}