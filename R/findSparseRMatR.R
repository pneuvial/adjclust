findSparseRMatR <- function(mat, p, h) {
    res <- .Call("CSparseRmatR", mat, as.integer(p), as.integer(h))
    return(res)
}