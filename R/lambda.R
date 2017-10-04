findLambda <- function(mat, p, h) {
    res <- .Call("CFindLambda", mat, as.integer(p), as.integer(h))
    return(res)  
}