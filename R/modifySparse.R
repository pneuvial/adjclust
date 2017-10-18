modifySparse <- function(m, p, h) {
    
    if (!(sparseCondnCheck(m, p, h))) {
        if (h < p-1) {
            stop("Some elements s[i,j] do not satisfy the condition s[i,j] <= 0.5*(s[i,i]+s[j,j]) and h<p-1; this case is not handled by the current implementation")
        } else {
            l <- findLambda(m, p, h)
            m[ ,1] <- m[ ,1] + (1.1*l) 
        }
    }
    
    if (!(all(m[ ,1] == 1))) {
        if (h < p-1) {
            stop("Some diagonal elements are not 1 and h<p-1; this case is not handled by the current implementation")
        } else {
            m <- makeSparseDiagOne(m, p, h)
        }
    }
    return(m)
}

sparseCondnCheck <- function(m, p, h) {
    res <- .Call("CSparseCondnCheck", m, as.integer(p), as.integer(h))
    return(res)
}

makeSparseDiagOne <- function(m, p, h) {
    res <- .Call("CMakeSparseDiagOne", m, as.integer(p), as.integer(h))
    return(res)
}