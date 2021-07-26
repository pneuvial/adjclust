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
#' @param dendro \code{\link{chac}} object as provided, e.g., by the function
#' \code{\link{adjClust}} (or any of the other wrappers).
#' @param palette color palette. Default to \code{\link{heat.colors}}
#' @param breaks number of breaks used to set colors from the palette. Those
#' are based on the quantiles of the matrix entries and for skewed distributions
#' the actual number used to set the palette can be lower than \code{breaks}.
#' @param log logical. Should the breaks be based on log-scaled values of the
#' matrix entries. Default to \code{TRUE}.
#' @param h if \code{mat} is of class \code{"snpMatrix"}, band parameter used to
#' compute the linkage disequilibrium (see \code{\link[snpStats]{ld}}).
#' @param stats if \code{mat} is of class \code{"snpMatrix"}, type of linkage
#' disequilibrium measure (see \code{\link[snpStats]{ld}}).
#' @param main graphic title.
#' @param col.clust color for the borders of the clusters (if \code{clustering}
#' is provided).
#' @param lwd.clust line width for the borders of the clusters (if 
#' \code{clustering} is provided).
#' @param xaxis logical. Should a x-axis be displayed? Default to \code{FALSE}
#' @param naxis number of breaks to display on the x-axis. For 
#' \code{HTCexp} objects, the axis is displayed in terms of Mpb and for the 
#' other types of input, it is displayed in terms of bin number. Default to 
#' \code{10}.
#' @importFrom grDevices heat.colors
#' @importFrom graphics image
#' @importFrom stats quantile
#' @importFrom stats rect.hclust
#' @importFrom graphics par
#' @examples
#' # input as HiTC::HTCexp object
#' \dontrun{
#' if (require("HiTC", quietly = TRUE)) {
#'   load(system.file("extdata", "hic_imr90_40_XX.rda", package = "adjclust"))
#'   plotSim(hic_imr90_40_XX)
#'   
#'   # with a constrained clustering
#'   res <- hicClust(hic_imr90_40_XX, log = TRUE)
#'   selected.capushe <- select(res)
#'   plotSim(hic_imr90_40_XX, clustering = selected.capushe, xaxis = TRUE)
#'   plotSim(hic_imr90_40_XX, clustering = selected.capushe, dendro = res)
#' }}
#' 
#' plotSim(dist(iris[ ,1:4]), log = FALSE)
#' @seealso \code{\link{select}}, \code{\link{adjClust}}
#' @export

plotSim <- function(mat, type = c("similarity", "dissimilarity"),
                    clustering = NULL, dendro = NULL, palette = heat.colors, 
                    breaks = 10, log = TRUE, h = NULL, 
                    stats = c("R.squared", "D.prime"), main = NULL, 
                    col.clust = "darkblue", lwd.clust = 2, xaxis = FALSE, 
                    naxis = 10) {
  UseMethod("plotSim")
}

extract_values <- function(mat) {
  UseMethod("extract_values")
}

extract_values.default <- function(mat) {
  p <- ncol(mat)
  all_coords <- expand.grid(1:p, 1:p)
  all_coords <- all_coords[all_coords[ ,1] <= all_coords[ ,2], ]
  x <- all_coords[ ,1] + all_coords[ ,2]
  y <- (all_coords[ ,2] - all_coords[ ,1]) + 1

  values <- mat[upper.tri(mat, diag = TRUE)]
  
  return(list("p" = p, "x" = x, "y" = y, "values" = values))
}

extract_values.dgCMatrix <- function(mat) {
  out <- extract_values_sparse(mat)
  return(out)
}

extract_values.dsCMatrix <- function(mat) {
  out <- extract_values_sparse(mat)
  return(out)
}

extract_values_sparse <- function(mat) {
  p <- ncol(mat)
  mat <- as(mat, "TsparseMatrix")
  x <- (mat@i + mat@j) + 2
  y <- (mat@j - mat@i) + 1
  values <- mat@x
  
  return(list("p" = p, "x" = x, "y" = y, "values" = values))
}

run_xaxis <- function(mat, dendro, naxis) {
  UseMethod("run_xaxis")
}

run_xaxis.default <- function(mat, dendro, naxis) {
  p <- ncol(mat)
  step <- floor((p-1) / naxis)
  all_breaks <- seq(1, p, by = step)
  
  all_pos <- 1.5 + (all_breaks - 1) / (p - 1) * (2*p - 1)
  if (is.null(dendro)) {
    axis_pos <- NA
  } else axis_pos <- 1.1*p
  axis(1, at = all_pos, labels = all_breaks, pos = axis_pos)
  
  invisible(NULL)
}

run_xaxis.HTCexp <- function(mat, dendro, naxis) {
  p <- ncol(mat)
  hic_ranges <- mat@xgi
  
  requireNamespace("BiocGenerics")
  startend <- c(min(BiocGenerics::start(hic_ranges@ranges)), 
                max(BiocGenerics::end(hic_ranges@ranges)))
  startend <- startend / 10^6
  step <- floor((floor(startend[2] - ceiling(startend[1]))) / naxis)
  all_breaks <- seq(ceiling(startend[1]), floor(startend[2]), by = step)
  
  where_1 <- (BiocGenerics::start(hic_ranges@ranges)[1] +
                BiocGenerics::start(hic_ranges@ranges)[1]) / 2 / 10^6
  where_p <- (BiocGenerics::start(hic_ranges@ranges)[p] +
                BiocGenerics::start(hic_ranges@ranges)[p]) / 2 / 10^6
  all_pos <- 1.5 + (all_breaks - where_1) / (where_p - where_1) * (2*p - 1)
  if (is.null(dendro)) {
    axis_pos <- NA
  } else axis_pos <- 1.1*p
  
  axis(1, at = all_pos, labels = paste0(all_breaks, "Mb"), pos = axis_pos)
  
  invisible(NULL) 
}

#' @export
plotSim.matrix <- function(mat, type = c("similarity", "dissimilarity"),
                           clustering = NULL, dendro = NULL, 
                           palette = heat.colors, breaks = 10, log = TRUE, 
                           h = NULL, stats = c("R.squared", "D.prime"), 
                           main = NULL, col.clust = "darkblue", lwd.clust = 2,
                           xaxis = FALSE, naxis = 10) {
  
  plotSim_generic(mat, type, clustering, dendro, palette, breaks, log, h, main, 
                  col.clust, lwd.clust, xaxis, naxis)
  
  invisible(NULL)
}

#' @export
plotSim.dgCMatrix <- function(mat, type = c("similarity", "dissimilarity"),
                              clustering = NULL, dendro = NULL, 
                              palette = heat.colors, breaks = 10, log = TRUE, 
                              h = NULL, stats = c("R.squared", "D.prime"), 
                              main = NULL, col.clust = "darkblue", 
                              lwd.clust = 2, xaxis = FALSE, naxis = 10) {
  
  plotSim_generic(mat, type, clustering, dendro, palette, breaks, log, h, main, 
                  col.clust, lwd.clust, xaxis, naxis)
  
  invisible(NULL)
}

#' @export
plotSim.dsCMatrix <- function(mat, type = c("similarity", "dissimilarity"),
                              clustering = NULL, dendro = NULL, 
                              palette = heat.colors, breaks = 10, log = TRUE, 
                              h = NULL, stats = c("R.squared", "D.prime"), 
                              main = NULL, col.clust = "darkblue", 
                              lwd.clust = 2, xaxis = FALSE, naxis = 10) {
  
  plotSim_generic(mat, type, clustering, dendro, palette, breaks, log, h, main, 
                  col.clust, lwd.clust, xaxis, naxis)
  
  invisible(NULL)
}

#' @export
plotSim.dist <- function(mat, type = c("similarity", "dissimilarity"),
                         clustering = NULL, dendro = NULL, 
                         palette = heat.colors, breaks = 10, log = TRUE, 
                         h = NULL, stats = c("R.squared", "D.prime"), 
                         main = NULL, col.clust = "darkblue", lwd.clust = 2, 
                         xaxis = FALSE, naxis = 10) {
  type <- match.arg(type)
  if (type != "dissimilarity") {
    message("Note: input class is 'dist' so 'type' is supposed to be 'dissimilarity'")
    type <- "dissimilarity"
  }
  
  mat <- as.matrix(mat)
  
  plotSim_generic(mat, type, clustering, dendro, palette, breaks, log, h, main, 
                  col.clust, lwd.clust, xaxis, naxis)
  
  invisible(NULL)
}

#' @export
plotSim.HTCexp <- function(mat, type = c("similarity", "dissimilarity"),
                           clustering = NULL, dendro = NULL, 
                           palette = heat.colors, breaks = 10, log = TRUE, 
                           h = NULL, stats = c("R.squared", "D.prime"), 
                           main = NULL, col.clust = "darkblue", lwd.clust = 2, 
                           xaxis = FALSE, naxis = 10) {
  if (!requireNamespace("HiTC")) {
    stop("Package 'HiTC' not available. 'HTCexp' input cannot be used.")
  }
  
  type <- match.arg(type)
  if (type != "similarity")
    stop("type 'dissimilarity' does not match 'HTCexp' data")
  
  mat <- HiTC::intdata(mat)
  plotSim_generic(mat, type, clustering, dendro, palette, breaks, log, h, main, 
                  col.clust, lwd.clust, xaxis, naxis)
  
  invisible(NULL)
}

plotSim.SnpMatrix <- function(mat, type = c("similarity", "dissimilarity"),
                              clustering = NULL, dendro = NULL, 
                              palette = heat.colors, breaks = 10, log = TRUE, 
                              h = NULL, stats = c("R.squared", "D.prime"), 
                              main = NULL, col.clust = "darkblue", 
                              lwd.clust = 2, xaxis = FALSE, naxis = 10) {
  if (!requireNamespace("snpStats")) {
    stop("Package 'snpStats' not available. 'SnpMatrix' input cannot be used.")
  }
  stats <- match.arg(stats)
  if (is.null(h)) h <- ncol(mat) - 1
  if (h >= ncol(mat)) {
    stop("h should be strictly less than p")
  }
  mat <- snpStats::ld(mat, stats = stats, depth = h)
  mat[mat > 1] <- 1  ## fix numerical aberrations
  mat[mat < 0] <- 0  ## fix numerical aberrations
  diag(mat) <- rep(1, nrow(mat))  ## by default the diagonal is 0 after 'snpStats::ld'
  if (any(is.na(mat))) {
    ww <- which(is.na(as.matrix(mat)), arr.ind = TRUE)
    warning("Clustering could not be performed due to missing value(s) or NaN(s) in LD estimates. Returning these indices")
    return(ww)
  }
  
  plotSim_generic(mat, type, clustering, dendro, palette, breaks, log, h, main, 
                  col.clust, lwd.clust, xaxis, naxis)
  
  invisible(NULL)
}

plotSim_generic <- function(mat, type = c("similarity", "dissimilarity"), 
                            clustering, dendro, palette, breaks, log, h, 
                            main, col.clust, lwd.clust, xaxis, naxis) {
  type <- match.arg(type)
  if (type == "dissimilarity") {
    mat <- max(mat) - mat
  }
  
  all_data <- extract_values(mat)
  if (!is.null(clustering)) {
    if (sum(!(diff(clustering) %in% c(0,1))))
      stop("'clustering' is not a contiguity constrained clustering.")
  }
  
  if (!is.null(dendro)) {
    if (class(dendro) != "chac")
      stop("'dendro' is not a 'chac' object.")
  }
  
  # initialize
  if (!is.null(clustering)) {
    if (length(clustering) != all_data$p)
      stop("'clustering' must have the same length that the similarity matrix.")
  }
  z <- matrix(0, ncol = all_data$p, nrow = 2 * all_data$p)
  
  pos_val <- all_data$values[all_data$values != 0]
  if (log) {
    bvalues <- quantile(log(pos_val + 1), probs = (0:(breaks))/(breaks))
    bvalues <- exp(bvalues) - 1
  } else {
    bvalues <- quantile(pos_val, probs = (0:(breaks))/(breaks))
  }
  bvalues <- c(0, unique(bvalues))
  all_colors <- c("white", rev(palette(length(bvalues) - 2)))
  z[cbind(all_data$x, all_data$y)] <- all_data$values
  z[cbind(all_data$x - 1, all_data$y)] <- all_data$values
  
  # plot matrix with or without dendro
  if (is.null(dendro)) {
    image(1:(2 * all_data$p), 1:all_data$p, z, breaks = bvalues, 
          col = all_colors, useRaster = TRUE, axes = FALSE, ylab = "", 
          xlab = "", main = main)
  } else {
    # dendro
    par(mfrow = c(2, 1))
    par(mar = c(0, 1, 2, 1), xaxs = "i")
    suppressWarnings(plot(dendro, center = TRUE, main = main, axes = FALSE))
    if (!is.null(clustering)) {
      rect.hclust(dendro, k = length(unique(clustering)))
    }
    
    # matrix plot
    if (xaxis) {
      ylim <- c(1, all_data$p * 1.1)
    } else {
      ylim <- c(1, all_data$p)
    }
    par(mar = c(0, 1, 0, 1))
    image(1:(2 * all_data$p), (1:all_data$p), z[ ,all_data$p:1], 
          breaks = bvalues, col = all_colors, useRaster = TRUE, axes = FALSE, 
          ylab = "", xlab = "", ylim = ylim)
  }
  
  # plot clustering
  if (!is.null(clustering)) {
    starts <- sapply(unique(clustering), function(aclust) {
      min(which(clustering == aclust))
    })
    ends <- sapply(unique(clustering), function(aclust) {
      max(which(clustering == aclust))
    })
    x0 <- c(2*starts - 1.5, 2*ends + 0.5)
    y0 <- rep(0.5, 2*length(starts))
    if (!is.null(dendro)) y0 <- (all_data$p + 1) - y0
    x1 <- rep(starts + ends - 0.5, 2)
    y1 <- rep(ends - starts + 2, 2)
    if (!is.null(dendro)) y1 <- (all_data$p + 1) - y1
    segments(x0, y0, x1, y1, lwd = lwd.clust, col = col.clust)
  }
  
  if (xaxis)
    run_xaxis(mat, dendro, naxis)
  
  invisible(NULL)
}
