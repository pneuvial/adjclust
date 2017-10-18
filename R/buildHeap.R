#' @useDynLib adjclust percDown
buildHeap <- function(heap, D, l) {
  for (ii in floor(l/2) : 1) {
    heap <- .Call("percDown", heap, D, as.integer(l), as.integer(ii), PACKAGE = "adjclust")
  }
  return(heap)
}
