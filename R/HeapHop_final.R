HeapHop_final <- function(mat, p, h, NbCLasses){
  p <- ncol(mat)
  
  x <- .Call("OneDiagBand",mat, as.integer(h))
  
  resP <- adjclust:::HeapHop(x, p, h, 1)
  
  gains <- t(resP)[,3]
  res <- t(resP)[,1:2]
  height <- cumsum(gains)
  traceW <- cbind(c(1:(p-1)),rev(height))
  tree <- list(traceW   = traceW,
                gains   = gains,
                merge   = res,
                height  = height,
                seqdist = height,
                order   = 1:p,
                labels  = paste("",1:p))
  class(tree) <- "hclust"
  return(tree)
}