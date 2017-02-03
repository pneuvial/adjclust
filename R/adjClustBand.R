#' Adjacency-constrained clustering of a band diagonal similarity matrix
#'
#' @importFrom matrixStats rowCumsums
#' @importFrom matrixStats colCumsums
#' @useDynLib adjclust
#'
#'
#' @export
#'
#' @param x A vector of length \eqn{p*h - h*(h+1)/2} similarities
#' @param p An integer, the number of items to be clustered
#' @param h An integer, the size of the diagonal band to be considered
#' @param flavor A character vector, the type of algorithm to be performed. Defaults to 'crayons' as described in Dehman (2015)
#' @param minNbBlocks An integer, the number of blocks at which clustering should be stopped (to save time). Defaults to 1, which corresponds to perform the clustering until all items are merged.
#' @param verbose A logical value: should extra information be displayed?
#'
#' @references Dehman A. (2015). "Spatial Clustering of Linkage Disequilibrium blocks for Genome-Wide Association Studies". PhD thesis. \url{https://tel.archives-ouvertes.fr/tel-01288568/}
#'
#' @examples
#'
#' data("R2.100", package="adjclust")
#' x <- R2.100
#'
#' h <- 100
#' p <- 603
#' res <- adjClustBand(x, p, h)
#' resK <- adjClustBand(x, p, h, flavor="PseudoMatrix")

adjClustBand <- function(x, p, h, flavor=c("crayons", "PseudoMatrix"), minNbBlocks=1, verbose=FALSE) {
    flavor <- match.arg(flavor)
    if (flavor=="crayons") {
        if (minNbBlocks>1) {
            stop("flavor 'crayons' not implemented (properly) when minNbMblocks>1")
        }
        res <- adjClustBand_heap(x, p, h, blMin=1, verbose=verbose)
    } else if (flavor=="PseudoMatrix") {
        resP <- HeapHop(x, p, h, minNbBlocks)
        mg <- t(resP[1:2, ])
        gains <- cumsum(resP[3, ])  ## seems to be broken?
        ord <- 1:(ncol(resP)+1)
        res <- list(
            traceW=NULL,
            gains=gains,
            merge=mg,
            height=gains,
            seqdist=gains,
            order=ord,
            labels=as.character(1:p),
            method="adhclust-PseudoMatrix",
            call=NULL,
            distMethod="Ward")
        class(res) <- "hclust"
    }
    return(res)
}
