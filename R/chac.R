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
#' @param mode type of dendrogram to plot (see Details). Default to 
#' \code{"standard"}
#' @details When \code{\link{plot.chac}} is called with 
#' \code{mode = "standard"}, the standard dendrogram is plotted, even though,
#' due to contingency constrains, some branches are reversed (decreasing
#' merges). When \code{\link{plot.chac}} is called with 
#' \code{mode = "corrected"}, a correction is applied to original heights so as
#' to have only non decreasing merges). It does not change the result of the 
#' clustering, only the look of the dendrogram for easier interpretation.\cr\cr
#' Other modes are provided that correspond to different alternatives
#' described in Grimm (1987): \itemize{
#' \item in \code{mode = "within-disp"}, heights correspond to whithin-cluster
#' dispersion, \emph{i.e.}, for a corresponding cluster, its height is 
#' \deqn{I(C) = \sum_{i \in C} d(i,g_C)} where \eqn{d} is the dissimilarity 
#' used to cluster objects and \eqn{g_C} is the center of gravity of cluster
#' \eqn{C}. In this case, heights are always non decreasing;
#' \item in \code{mode = "total-disp"}, heights correspond to the sum of all
#' whithin-cluster dispersions. It is obtained from \code{mode = "standard"} by
#' the cumulartive sum of its heights. In this case, heights are always
#' non decreasing;
#' \item in \code{mode = "average-disp"}, heights correspond to the 
#' whithin-cluster dispersion divided by the cluster size. In this case, there 
#' is no guaranty that the heights are non decreasing. When reversals are 
#' detected, a warning is printed to advice the user to change the mode of the
#' representation.
#' }\cr
#' Grimm (1987) indicates that heights as provided by 
#' \code{mode = "within-disp"} are highly dependant on cluster sizes and that 
#' the most advisable representation is the one provided by 
#' \code{mode = "total-disp"}.
#' @references { Grimm, E.C. (1987) CONISS: a fortran 77 program for
#' stratigraphically constrained analysis by the method of incremental sum of
#' squares. \emph{Computer & Geosciences}, \strong{13}(1), 13-35. }
#' @return The function \code{plot.chac} displays the dendrogram and 
#' additionally invisibly returns an object of class 
#' \code{\link[stats]{dendrogram}} with heights as specified by the user through
#' the option \code{mode}.
#' @export
#' @importFrom graphics plot segments
#' @importFrom stats as.dendrogram cutree
plot.chac <- function(x, y, ..., 
                      mode = c("standard", "corrected", "total-disp", 
                               "within-disp", "average-disp")) {
  mode <- match.arg(mode)
  args <- list(...)
  if (is.null(args$type)) args$type <- "triangle"
  if (is.null(args$leaflab)) args$leaflab <- "none"
  
  if (mode == "standard") {
    if (any(diff(x$height) < 0)) 
      warning(paste0("\nDetected reversals in dendrogram: ",
                     "mode = 'corrected', 'within-disp' or 'total-disp' might be more relevant."))
    if (is.null(args$ylim)) args$ylim <- range(x$height)
  } else if (mode == "corrected") {
    res_diagnose <- diagnose(x, graph = FALSE, verbose = FALSE)
    to_add <- data.frame(res_diagnose$number,
                         add = res_diagnose$pheight - res_diagnose$height)
    to_add <- apply(to_add, 1, function(acol) {
      c(rep(0, acol[1] - 1), rep(acol[2], length(x$height) - acol[1] + 1))
    })
    to_add <- rowSums(to_add)
    x$height <- x$height + to_add
    ## note: remaining decreasing gains due to numerical approximations
  } else if (mode == "total-disp") {
    x$height <- cumsum(x$height)
  } else if (mode == "within-disp" | mode == "average-disp") {
    to_correct <- which((x$merge[ ,1] > 0) | (x$merge[ ,2] > 0))
    for (ind in to_correct) {
      clusters <- x$merge[ind, ]
      clusters <- clusters[clusters > 0]
      x$height[ind] <- x$height[ind] + sum(x$height[clusters])
    }
    if (mode == "average-disp") {
      # search for current cluster size
      out <- sapply((length(x$height) + 1):1, function(num) {
        tmp <- table(table(cutree(as.hclust(x), k = num)))
        res <- rep(0, length(x$height) + 1)
        res[as.numeric(names(tmp))] <- tmp
        return(res)
      })
      out <- apply(out, 1, diff)
      sizes <- unlist(apply(out, 1, function(acol) which(acol > 0)))
      sizes <- as.numeric(sizes)
      x$height <- x$height / sizes
      if (any(diff(x$height) < 0)) 
        warning(paste0("\nDetected reversals in dendrogram: ",
                       "mode = 'corrected', 'within-disp' or 'total-disp' might be more relevant."))
      if (is.null(args$ylim)) args$ylim <- range(x$height)
    }
    args$ylab <- "within-cluster dispersion"
  }
  
  args$x <- as.dendrogram(as.hclust(x))
  do.call(plot, args)
  
  # for "mode='corrected'", show the corrections
  if (mode == "corrected") {
    x0 <- rep(0, nrow(res_diagnose))
    x1 <- rep(0, nrow(res_diagnose))
    y1 <- x$height[res_diagnose$number]
    y0 <- y1 - (res_diagnose$pheight - res_diagnose$height)
    segments(x0, y0, x1, y1, col = "darkred", lwd = 2)
    segments(rep(-0.25, length(x0)), y0, rep(0.25, length(x0)), col = "darkred")
    segments(rep(-0.25, length(x0)), y1, rep(0.25, length(x0)), col = "darkred")
  }
  
  invisible(as.dendrogram(as.hclust(x)))
}

#' @export
diagnose <- function(x, ...) {
  UseMethod("diagnose")
}

#' @rdname chac
#' @aliases diagnose.chac
#' @aliases diagnose
#' @param graph (logical) whether the diagnostic plot has to be displayed or 
#' not. Default to \code{TRUE}
#' @param verbose (logical) whether to print a summary of the result or not.
#' Default to \code{TRUE}
#' @return \code{\link{diagnose}} invisibly exports a data frame with the 
#' numbers of decreasing merges described by the labels of the clusters being
#' merged at this step and at the previous one, as well as the corresponding
#' merge heights.
#' @importFrom graphics axis points
#' @export
diagnose.chac <- function(x, graph = TRUE, verbose = TRUE) {
  diff_heights <- diff(x$height)
  if (any(diff_heights < 0)) {
    if (verbose)
      cat(sum(diff_heights < 0), "merges with decreasing heights:", "\n")
    
    # extract decreasing merges with the previous merge
    where_decrease <- which(diff_heights < 0)
    res <- sapply(where_decrease, function(adec) {
      out <- c(adec+1, x$merge[adec+1, ], x$height[adec+1], x$merge[adec, ],
               x$height[adec])
      names(out) <- c("number", "x1", "x2", "height", "px1", "px2", "pheight")
      return(out)
    })
    res <- data.frame(t(res))
    if (verbose) {
      print(head(res))
      if (nrow(res) > 6) cat("...", nrow(res) - 6, "remaining rows...", "\n")
    }
    
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
    print("All merges have non decreasing heights.")
    invisible(NULL)
  }
}
