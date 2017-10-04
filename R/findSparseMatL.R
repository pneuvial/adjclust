findSparseMatL <- function(mat, p, h) {
    res <- .Call("CSparseMatL" , mat, as.integer(p), as.integer(h))
    return(res)
}