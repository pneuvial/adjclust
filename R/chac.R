#' Class chac
#' 
#' S3 class for Constrained Hierarchical Agglomerative Clustering results
#' 
#' 
#
#' Methods for class 'chac'
#'
#' @name chac
#' @param x,object,tree an object of class 'chac'
#' @param y not used
#' @param ... for \code{\link{plot}}, arguments passed to the function 
#' \code{\link{plot.dendrogram}}. Default values for \code{type} and
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
#' \item in \code{mode = "within-disp"}, heights correspond to within-cluster
#' dispersion, \emph{i.e.}, for a corresponding cluster, its height is 
#' \deqn{I(C) = \sum_{i \in C} d(i,g_C)} where \eqn{d} is the dissimilarity 
#' used to cluster objects and \eqn{g_C} is the center of gravity of cluster
#' \eqn{C}. In this case, heights are always non decreasing;
#' \item in \code{mode = "total-disp"}, heights correspond to the total
#' within-cluster dispersion. It is obtained from \code{mode = "standard"} by
#' the cumulative sum of its heights. In this case, heights are always
#' non decreasing;
#' \item in \code{mode = "average-disp"}, heights correspond to the 
#' within-cluster dispersion divided by the cluster size. In this case, there 
#' is no guaranty that the heights are non decreasing. When reversals are 
#' detected, a warning is printed to advice the user to change the mode of the
#' representation.}
#' Grimm (1987) indicates that heights as provided by 
#' \code{mode = "within-disp"} are highly dependent on cluster sizes and that 
#' the most advisable representation is the one provided by 
#' \code{mode = "total-disp"}. Further details are provided in the vignette 
#' "Notes on CHAC implementation in adjclust".
#' @param nodeLabel (logical) whether the order of merging has to be displayed
#' or not. \code{nodeLabel=TRUE} prints orders of fusion at corresponding
#' nodes. Default to \code{FALSE}
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
                               "within-disp", "average-disp"), 
                      nodeLabel = FALSE) {
  mode <- match.arg(mode)
  args <- list(...)
  if (is.null(args$type)) args$type <- "triangle"
  if (is.null(args$leaflab)) args$leaflab <- "none"
  
  if ((x$method == "adjClust-corrected") & (mode != "standard")) {
    stop("Already corrected 'chac' object. 'mode' must be set to 'standard'")
  }
  
  if (mode == "standard") {
    if (any(diff(x$height) < 0)) 
      warning(paste0("\nDetected reversals in dendrogram: ",
                     "mode = 'corrected', 'within-disp' or 'total-disp' might be more relevant."))
    if (is.null(args$ylim)) args$ylim <- range(x$height)
  } else if (mode == "corrected") {
    res_diagnose <- diagnose(x, graph = FALSE, verbose = FALSE)
    x <- correct(x)
    args$ylab <- "corrected"
  } else if (mode == "total-disp") {
    x$height <- cumsum(x$height)
    args$ylab <- "total dispersion"
  } else if (mode == "within-disp" | mode == "average-disp") {
    args$ylab <- "within-cluster dispersion"
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
      args$ylab <- "average dispersion"
    }
  }
  if (nodeLabel) {
	  args$x <- alt.as.dendrogram(as.hclust(x))
	  do.call(alt.plot, args)
  } else {
	  args$x <- as.dendrogram(as.hclust(x))
	  do.call(plot, args)
  }
  
  # for "mode='corrected'", show the corrections
  if (mode == "corrected") {
    x0 <- rep(0, nrow(res_diagnose))
    x1 <- rep(0, nrow(res_diagnose))
    y1 <- x$height[res_diagnose$number]
    y0 <- y1 - (res_diagnose$pheight - res_diagnose$height) * 1.0001
    segments(x0, y0, x1, y1, col = "darkred", lwd = 2)
    segments(rep(-0.25, length(x0)), y0, rep(0.25, length(x0)), col = "darkred")
    segments(rep(-0.25, length(x0)), y1, rep(0.25, length(x0)), col = "darkred")
  }
  
  invisible(as.dendrogram(as.hclust(x)))
}

#' @rdname chac
#' @aliases diagnose
#' @aliases diagnose.chac
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
diagnose <- function(x, graph = TRUE, verbose = TRUE) {
  UseMethod("diagnose")
}

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
    message("All merges have non decreasing heights.")
    invisible(NULL)
  }
}

#' @rdname chac
#' @aliases correct
#' @aliases correct.chac
#' @return The function \code{\link{correct}} returns a \code{chac} objects with
#' modified heights so as they are increasing. The new heights are calculated in
#' an way identical to the option \code{mode = "corrected"} of the function
#' \code{plot.chac} (see Details). In addition, the \code{chac} object has its
#' field \code{method} modified from \code{adjClust} to 
#' \code{adjClust-modified}.
#' @export
correct <- function(x) {
  UseMethod("correct")
}

#' @export
correct.chac <- function(x) {
  if (any(diff(x$height) < 0)) {
    res_diagnose <- diagnose(x, graph = FALSE, verbose = FALSE)
    to_add <- data.frame(res_diagnose$number,
                         add = res_diagnose$pheight - res_diagnose$height)
    to_add$add <- 1.0001 * to_add$add

    value <- rep(0, length(x$height))
    value[to_add[, 1]] <- to_add[,2]
    x$height <- x$height + cumsum(value)
    
    x$method <- "adjClust-corrected"
    return(x) 
  } else {
    warning("No reversal. Returned nothing.")
    invisible(NULL)
  }
}

#' @rdname chac
#' @aliases cuttree_chac
#' @param k an integer scalar or vector with the desired number of groups
#' @param h numeric scalar or vector with heights where the tree should be cut.
#' Only available when the heights are increasing
#' @return The function \code{\link{cutree_chac}} returns the clustering with 
#' \code{k} groups or with the groups obtained by cutting the tree at height
#' \code{h}. If the heights are not increasing, the cutting of the tree is based
#' on the corrected heights as provided by the function \code{correct}.
#' @export
#' 
cutree_chac <- function(tree, k = NULL, h = NULL) {
  if (class(tree) != "chac")
    stop("'tree' must be of class 'chac'")
  
  if (any(diff(tree$height) < 0)) {
    if (is.null(k)) {
      stop("With decreasing heights, 'k' must be provided.")
    } else {
      tree <- correct(tree)
    }
  }
  tree <- as.hclust(tree)
  res <- cutree(tree, k = k, h = h)
  return(res)
}

#' @name select
#' @aliases select.chac
#' @title Clustering selection
#' @description Clustering selection from a chac object with the slope heuristic
#' or the broken stick heuristic
#' @param x an object of class 'chac'
#' @param type model selection approach between slope heuristic 
#' (\code{"capushe"}) and broken stick approach (\code{"bstick"})
#' @param k.max maximum number of clusters that can be selected. Default to 
#' \code{NULL}, in which case it is set to 
#' \eqn{\min(\max(100, \frac{n}{\log(n)}), \frac{n}{2})} where \eqn{n} is the
#' number of objects to be clustered for capushe and to \eqn{n} for the broken
#' stick model
#' @param graph logical. Whether the diagnostic plot for the capushe selection
#' is displayed or not. Default to \code{FALSE}
#' @param pct minimum percentage of points for the plateau selection in 
#' capushe selection. See \code{\link[capushe]{DDSE}} for further details
#' @return The function returns the clustering selected by the slope heuristic,
#' as implemented in the R package \code{capushe}.
#' @importFrom capushe DDSE
#' @importFrom capushe Djump
#' @importFrom graphics lines
#' @references Baudry, J.P., Maugis, C. and Michel, B. (2012) Slope heuristics: 
#' overview and implementation. \emph{Statistics and Computing}, \strong{22}(2),
#' 355-470.
#' MacArthur, R.H. (1957) On the relative abundance of bird species. 
#' \emph{Proceedings of the National Academy of Sciences}, \strong{43}, 293-295.
#' @examples \dontrun{if (require("HiTC", quietly = TRUE)) {
#'   load(system.file("extdata", "hic_imr90_40_XX.rda", package = "adjclust"))
#'   res <- hicClust(hic_imr90_40_XX, log = TRUE)
#'   selected.capushe <- select(res)
#'   table(selected.capushe)
#'   selected.bs <- select(res, type = "bstick")
#'   table(selected.bs)
#' }}
#' 
#' res <- adjClust(dist(iris[ ,1:4]))
#' select.clust <- select(res, "bs")
#' table(select.clust)
#' 
#' @export

select <- function(x, type = c("capushe", "bstick"), k.max = NULL, 
                   graph = FALSE, pct = 0.15) {
  UseMethod("select")
}

#' @export
select.chac <- function(x, type = c("capushe", "bstick"), k.max = NULL, 
                        graph = FALSE, pct = 0.15) {
  type <- match.arg(type)
  n <- length(x$labels)

  if (type == "capushe") {
    if (is.null(k.max)) {
      k.max <- round(min(max(100, n/log(n)), n/2))
    }
    
    in_capushe <- data.frame(name = 1:k.max, 
                             pen.shape = lchoose(n - 1, 0:(k.max-1)),
                             complexity = 1:k.max, 
                             contrast = cumsum(x$height)[(n-1):(n-k.max)])
    
    KC <- try(DDSE(in_capushe, pct = pct), silent = TRUE)
    
    if (class(KC) == "try-error")
      KC <- Djump(in_capushe)
    
    if (graph)
      plot(KC)
    
    res <- cutree_chac(x, k = as.integer(KC@model))
  } else if (type == "bstick") {
    if (is.null(k.max)) k.max <- n
    disp <- rev(x$height)
    tot.disp <- sum(disp)
    disp <- abs(disp)
    bs <- tot.disp * rev(cumsum(1/((n-1):1))/(n-1))
    cdt <- which(disp[1:(k.max-1)] <= bs[1:(k.max-1)])
    if (length(cdt) > 0) {
      KC <- min(cdt)
    } else {
      KC <- k.max
    }
    
    if (graph) {
      plot(bs, type = "l", xlab = "number of clusters", ylab = "broken stick")
      lines(seq_along(bs)[1:KC], disp[1:KC], type = "b", pch = "+", 
            col = "darkgreen")
      lines(seq_along(bs)[KC:length(bs)], disp[KC:length(bs)], type = "b", 
            pch = "+", col = "red")
    }
    
    res <- cutree_chac(x, k = KC)
  }
  return(res)
}