#' Plot (dis)similarity matrix
#' 
#' Heatmap of the (dis)similarity matrix
#' 
#' This function produces a heatmap for the used (dis)similarity matrix that 
#' can be used as a diagnostic plot to check the consistency between the 
#' obtained clustering and the original (dis)similarity
#' 
#' @param mat matrix to plot. It can be of class \code{'matrix'}, 
#' \code{'dgCMatrix'}, \code{'dsCMatrix'}, \code{'dist'}, \code{'HTCexp'},
#' \code{'snpMatrix'}.
#' @param type input matrix type. Can be either \code{"similarity"} or 
#' \code{"dissimilarity"} (kernels are supposed to be of type 
#' \code{"similarity"}).
#' @param clustering vector of length the number of rows (columns) of the 
#' matrix that contains a contiguity constrained clustering (as provided by
#' \code{\link{select}} for instance). If supplied the clustering is 
#' superimposed over the heatmap.
#' @param palette color palette. Default to \code{\link[grDevices]{heat.colors}}
#' @param breaks number of breaks used to set colors from the palette. Those
#' are based on the quantiles of the matrix entries and for skewed distributions
#' the actual number used to set the palette can be lower than \code{breaks}.
#' @param log logical. Should the breaks be based on log-scaled values of the
#' matrix entries. Default to \code{TRUE}.
#' @param h if \code{mat} is of class \code{"snpMatrix"}, band parameter used to
#' compute the linkage desiquilibrium (see \code{\link[snpStats]{ld}}).
#' @param stats if \code{mat} is of class \code{"snpMatrix"}, type of linkage
#' desiquilibrium measure (see \code{\link[snpStats]{ld}}).
#' @param main graphic title.
#' @param col.clust color for the borders of the clusters (if \code{clustering}
#' is provided).
#' @param lwd.clust line width for the borders of the clusters (if 
#' \code{clustering} is provided).
#' @importFrom grDevices heat.colors
#' @importFrom graphics image
#' @importFrom stats quantile
#' @examples
#' # input as HiTC::HTCexp object
#' if (require("HiTC", quietly = TRUE)) {
#'   load(system.file("extdata", "hic_imr90_40_XX.rda", package = "adjclust"))
#'   plotSim(hic_imr90_40_XX)
#'   
#'   # with a constrained clustering
#'   res <- hicClust(hic_imr90_40_XX, log = TRUE)
#'   selected.capushe <- select(res)
#'   plotSim(hic_imr90_40_XX, clustering = selected.capushe)
#' }
#' plotSim(dist(iris[ ,1:4]), log = FALSE)
#' @seealso \code{\link{select}}
#' @export

plotSim <- function(mat, type = c("similarity", "dissimilarity"),
                    clustering = NULL, palette = heat.colors, breaks = 10, 
                    log = TRUE, h = p - 1, stats = c("R.squared", "D.prime"),
                    main = NULL, col.clust = "darkblue", lwd.clust = 2) {
  # checks
  type <- match.arg(type)
  stats <- match.arg(stats)

  CLASS <- c("matrix", "dgCMatrix", "dsCMatrix", "dist", "HTCexp", 
             "snpMatrix")
  classcheck <- pmatch(class(mat), CLASS)
  if (is.na(classcheck)) {
    stop("Input matrix class not supported")
  }
  if (classcheck == -1) {
    stop("Ambiguous matrix class")
  }
  inclass <- CLASS[classcheck]
  
  if (inclass == "dist") {
    mat <- as.matrix(mat)
    if (type != "dissimilarity") {
      type <- "dissimilarity"
      message("type 'dissimilarity' used for objects of class 'dist'")
    }
  }
  
  if (!is.null(clustering)) {
    if (sum(!(diff(clustering) %in% c(0,1))))
      stop("'clustering' is not a contiguity constrained clustering.")
  }
  
  # special cases and preprocessing
  if (inclass == "HTCexp") {
    if (!requireNamespace("HiTC")) {
      stop("Package 'HiTC' not available. 'HTCexp' input cannot be used.")
    }
    mat <- HiTC::intdata(mat)
    if (type != "similarity")
      stop("type 'dissimilarity' does not match 'HTCexp' data")
  }
  if (inclass == "SnpMatrix") {
    if (!requireNamespace("snpStats")) {
      stop("Package 'snpStats' not available. 'SnpMatrix' input cannot be used.")
    }
    if (h >= p) {
      stop("h should be strictly less than p")
    }
    mat <- snpStats::ld(mat, stats = stats, depth = h)
    mat[mat > 1] <- 1  ## fix numerical aberrations
    mat[mat < 0] <- 0  ## fix numerical aberrations
    diag(mat) <- rep(1, nrow(mat))  ## by default the diagonal is 0 after 'snpStats::ld'
    #x <- round(x, digits = 10) ## ensure ascending compatibility but removed for sanity
    if (any(is.na(x))) {
      ww <- which(is.na(as.matrix(x)), arr.ind = TRUE)
      warning("Clustering could not be performed due to missing value(s) or NaN(s) in LD estimates. Returning these indices")
      return(ww)
    }
  }
  
  # extract dimension and initialize
  p <- nrow(mat)
  z <- matrix(0, ncol = p, nrow = 2*p)
  
  if (type == "dissimilarity") {
    mat <- max(mat) - mat
    mat <- as(mat, "matrix")
  }
  
  # extract data
  if (class(mat) == "matrix") {
    all_coords <- expand.grid(1:p, 1:p)
    all_coords <- all_coords[all_coords[ ,1] <= all_coords[ ,2], ]
    x <- all_coords[ ,1] + all_coords[ ,2]
    y <- (all_coords[ ,2] - all_coords[ ,1]) + 1
    values <- mat[upper.tri(mat, diag = TRUE)]
  } else if (class(mat) == "dsCMatrix") {
    mat <- as(mat, "TsparseMatrix")
    x <- (mat@i + mat@j) + 2
    y <- (mat@j - mat@i) + 1
    values <- mat@x
  }
  
  pos_val <- values[values != 0]
  if (log) {
    bvalues <- quantile(log(pos_val + 1), probs = (0:(breaks))/(breaks))
    bvalues <- exp(bvalues) - 1
  } else {
    bvalues <- quantile(pos_val, probs = (0:(breaks))/(breaks))
  }
  bvalues <- c(0, unique(bvalues))
  all_colors <- c("white", rev(palette(length(bvalues) - 2)))
  z[cbind(x, y)] <- values
  z[cbind(x-1, y)] <- values
  
  # plot
  image(1:(2*p+1), 1:p, z, breaks = bvalues, col = all_colors, 
        useRaster = TRUE, axes = FALSE, ylab = "", xlab = "", main = main)
  
  if (!is.null(clustering)) {
    starts <- sapply(unique(clustering), function(aclust) {
      min(which(clustering == aclust))
    })
    ends <- sapply(unique(clustering), function(aclust) {
      max(which(clustering == aclust))
    })
    x0 <- c(2*starts - 1.5, 2*ends + 0.5)
    y0 <- rep(0.5, 2*length(starts))
    x1 <- rep(starts + ends - 0.5, 2)
    y1 <- rep(ends - starts + 2, 2)
    segments(x0, y0, x1, y1, lwd = lwd.clust, col = col.clust)
  }
  
  invisible(NULL)
}
