#' @export
transpose <- function(x, p, h) {
    len <- length(x)
    stopifnot(len==(p-1)*h-h*(h-1)/2)
    tx <- numeric(len)
    kk <- 0
    for (ii in 1:(p-1)){
        for (jj in (ii+1):min(p,(ii+h))){ ## NB: jj loops over less than h terms!
            kk <- kk+1
            if (jj>h) {
                idx <- (h*(h-1))/2 + (jj-h-1)*h +h+1-(jj-ii)
            } else {
                idx <- ((jj-2)*(jj-1))/2 + ii
            }
            tx[kk] <- x[idx]
        }
    }
    return(tx)
}
