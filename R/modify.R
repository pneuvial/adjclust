modify <- function(m, p, h) {
    if (!(condnCheck(m))) {
        if (h < p-1) {
            stop("Some elements s[i,j] do not satisfy the condition s[i,j] <= 0.5*(s[i,i]+s[j,j]) and h<p-1; this case is not handled by the current implementation")
        } else {
            m1 <- sweep(m, 1, diag(m))
            l <- max(m1)  
            diag(m) <- diag(m) + l + (0.1*l)   #Taking epsilon = 0.1*l
        }
    }
    
    if (!(all(diag(m) == 1))) {
        if (h < p-1) {
            stop("Some diagonal elements are not 1 and h<p-1; this case is not handled by the current implementation")
        } else {
            m <- makeDiagOne(m)
        }
    }
    return(m)
}

condnCheck <- function(m) {
    res <- .Call("CCondnCheck", m)
    return(res)
}

makeDiagOne <- function(m) {
    res <- .Call("CMakeDiagOne", m)
    return(res)
}
