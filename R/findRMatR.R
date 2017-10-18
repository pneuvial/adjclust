findRMatR <- function(mat, p, h) {
    res <- .Call("CRmatR", mat, as.integer(p), as.integer(h))
    return(res)
}