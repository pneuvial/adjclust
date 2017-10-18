findMatL <- function(mat, p, h) {
    res <- .Call("CmatL", mat, as.integer(p), as.integer(h))
    return(res)
}