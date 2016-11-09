#' Adjacency-constrained clustering of a band diagonal similarity matrix
#'
#' @importFrom matrixStats rowCumsums
#' @importFrom matrixStats colCumsums
#' @export
#'
#' @examples
#'
#' data("ld_ceph", package="adjclust")
#' x <- R2.100
#' x <- Dprime.100
#'
#' h <- 100
#' p <- 603
#' res <- adjClustBand(x, p, h)
#' resP <- adjClustBand(x, p, h, flavor="pseudoMatrix")

adjClustBand <- function(x, p, h, flavor=c("heap", "pseudoMatrix"), minNbBlocks=1, trace.time=FALSE, verbose=FALSE) {
    flavor <- match.arg(flavor)
    if (flavor=="heap") {
        if (minNbBlocks>1) {
            stop("flavor 'heap' not implemented (properly) when minNbMblocks>1")
        }
        res <- adjClustBand_heap(x, p, h, blMin=1, trace.time=trace.time, verbose=verbose)
    } else if (flavor=="pseudoMatrix") {
        resP <- HeapHop(x, p, h, minNbBlocks)
        mg <- t(resP[1:2, ])
        gains0 <- resP[3, ]  ## seems to be broken?
        gains <- 1:(p-minNbBlocks)  ## FIXME: terrible hack here
        res <- list(
            traceW=NULL,
            gains=gains,
            gains0=gains0, ## for debugging purposes
            merge=mg,
            height=gains,
            seqdist=gains,
            order=1:p,
            labels=as.character(1:p),
            method="adhclust-pseudoMatrix",
            call=NULL,
            distMethod="Ward")
        class(res) <- "hclust"
    }
    return(res)
}
