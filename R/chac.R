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
    if (is.null(args$ylim)) args$ylim <- range(x$height)
    args$x <- as.dendrogram(as.hclust(x))
    if (is.null(args$type)) args$type <- "triangle"
    if (is.null(args$leaflab)) args$leaflab <- "none"
    do.call(plot, args)
}


#' @rdname chac
#' @aliases diagnose
#' @export
diagnose <- function(x, ...) {
  UseMethod("diagnose")
}

#' @rdname chac
#' @aliases diagnose.chac
#' @param graph (logical) whether the diagnostic plot has to be displayed or 
#' not. Default to \code{TRUE}
#' @details \code{\link{diagnose}} invisibly exports a data frame with the 
#' numbers of non increasing merges described by the labels of the clusters 
#' being merged at this step and at the previous one, as well as the 
#' corresponding merge heights.
#' @export
diagnose.chac <- function(x, graph = TRUE) {
  diff_heights <- diff(x$height)
  if (any(diff_heights < 0)) {
    cat(sum(diff_heights < 0), "merges with non increasing heights:", "\n")
    
    # extract non increasing merges with the previous merge
    where_decrease <- which(diff_heights < 0)
    res <- sapply(where_decrease, function(adec) {
      out <- c(adec+1, x$merge[adec+1, ], x$height[adec+1], x$merge[adec, ],
               x$height[adec])
      names(out) <- c("number", "x1", "x2", "height", "px1", "px2", "pheight")
      return(out)
    })
    res <- data.frame(t(res))
    print(head(res))
    if (nrow(res) > 6) cat("...", nrow(res) - 6, "remaining rows...", "\n")
    
    # if requested, plot
    if (graph) {
      plot(x$height, type = "b", pch = "+", axes = FALSE, xlab = "Merge number",
           ylab = "height of merge")
      
      axis(1, at = where_decrease + 1, cex = 0.7)
                            
      axis(2, tick = FALSE)
      axis(2, at = x$height[where_decrease + 1], labels = FALSE, 
           col.ticks = "red")
      
      points(where_decrease + 1, x$height[where_decrease + 1], col = "red",
             pch = "+")
    }
    
    invisible(res)
  } else {
    print("All merges have increasing weights.")
    return(NULL)
  }
}
