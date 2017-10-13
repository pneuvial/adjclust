#' Methods for class 'chac'
#'
#' @name chac
#' @param x,object an object of class 'chac'
#' @param y not used
#' @param ... for \code{\link{plot}}, arguments passed to the function 
#' \code{\link[stats]{plot.dendrogram}}. Default values for \code{type} and
#' \code{leaflab} are respectively set to \code{"triangle"} and \code{"none"}
NULL

#' @rdname chac
#' @aliases as.hclust.chac
#' @importFrom stats as.hclust
#' @export
as.hclust.chac <- function(x, ...) {
    res <- x
    class(res) <- "hclust"
    return(res)
}

#' @rdname chac
#' @aliases print.chac
#' @export
print.chac <- function(x, ...) {
    x <- as.hclust(x)
    print(x)
}

#' @rdname chac
#' @aliases head.chac
#' @importFrom utils head
#' @export
head.chac <- function(x, ...) {
    sapply(x, head)
}

#' @rdname chac
#' @aliases summary.chac
#' @export
summary.chac <- function(object, ...) {
    print(object)
}

#' @rdname chac
#' @aliases plot.chac
#' @export
#' @importFrom graphics plot
#' @importFrom stats as.dendrogram
plot.chac <- function(x, y, ...) {
    args <- list(...)
    args$x <- as.dendrogram(as.hclust(x))
    if (is.null(args$type)) args$type <- "triangle"
    if (is.null(args$leaflab)) args$leaflab <- "none"
    do.call(plot, args)
}
