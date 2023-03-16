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
    if (!inherits(dendro, "chac")) stop("'dendro' is not a 'chac' object.")
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

#### ggplot2 version (work in progress)
# HERE: WORK IN PROGRESS #####
# WORKING ON ISSUE 40 ####

take_bottomleft <- function(amat) amat[lower.tri(amat, diag = TRUE)]

make_coords <- function(indi, indj, values) {
  if (any(is.na(values))) {
    selected <- !is.na(values)
    indi <- indi[selected]
    indj <- indj[selected]
  }
  
  # defining polygon borders (square matrix) and rotate them
  coords <- cbind(rep(seq_along(indi), 4),
                  c(indi - 1, indi, indi, indi - 1),
                  c(indj, indj, indj - 1, indj - 1))
  theta <- 3*pi/4
  rotation_matrix <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)),
                            ncol = 2)
  coords <- cbind(coords[, 1], rep(indi, 4), rep(indj, 4), 
                  t(tcrossprod(rotation_matrix, as.matrix(coords[, 2:3]))))
  coords <- data.frame(coords, rep(values, 4))
  names(coords) <- c("id", "i", "j", "x", "y", "IF")
  coords$x <- - coords$x
  
  # remove auto-IF bottom polygon border
  bottom_poly_rm <- (coords$i == coords$j) & (coords$y < 0)
  coords <- coords[!bottom_poly_rm, ]
  
  return(coords)
}

rescale_coords <- function(x, y, ymax_i, ymax_f, xmin_i, xmax_i, xmax_f) {
  # ymax <- max(fake_coords$y)
  y <- y / ymax_i * ymax_f
  y <- - y
  
  x <- (xmax_f^2 - 1) / xmax_f / (xmax_i - xmin_i) * x
  x <- x + (xmax_i * (xmax_f+1) + xmin_i * (-2 * xmax_f^2 - xmax_f + 1)) / 
    (2 * xmax_f * (xmax_i - xmin_i))
  
  out <- data.frame("x" = x, "y" = y)
  return(out)
}

poly_coords_sparse <- function(mat) {
  indi <- mat@j + 1
  indj <- mat@i + 1
  values <- mat@x
  selected <- indi >= indj
  indif <- indi[selected]
  indjf <- indj[selected]
  valuesf <- values[selected]
  
  coords <- make_coords(indif, indjf, valuesf)
  
  return(coords)
}

poly_coords <- function(mat) {
  UseMethod("poly_coords")
}

poly_coords.default <- function(mat) {
  # extracting coordinates in the matrix (genomic) and IF
  p <- ncol(mat)
  indi <- row(mat)
  indi <- take_bottomleft(indi)
  indj <- col(mat)
  indj <- take_bottomleft(indj)
  values <- take_bottomleft(mat)
  
  coords <- make_coords(indi, indj, values)
  
  return(coords)
}

poly_coords.dsCMatrix <- function(mat) {
  p <- ncol(mat)
  mat <- as(mat, "TsparseMatrix")
  coords <- poly_coords_sparse(mat)
  
  return(coords)
}

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
#' @param log logical. Should the breaks be based on log-scaled values of the
#' matrix entries. Default to \code{TRUE}.
#' @param legendName character. Title of the legend. Default to 
#' \code{"intensity"}.
#' @param main character. Title of the plot. Default to \code{NULL} (no title).
#' @param priorCount numeric. Average count to be added to each entry of the
#' matrix to avoid taking log of zero. Used only if \code{log = TRUE}.
#' @import ggplot2
#' @examples \dontrun{
#' clustering <- rep(1:3, each = 50)
#' dist_data <- as.matrix(dist(iris[, 1:4]))
#' dendro_iris <- adjClust(dist_data, type = "dissimilarity")
#' ggPlotSim(dist_data, type = "dissimilarity", dendro = dendro_iris,
#'           axis = TRUE)
#' ggPlotSim(dist_data, type = "dissimilarity", dendro = dendro_iris,
#'           clustering = clustering)
#' ggPlotSim(dist_data, type = "dissimilarity", legendName = "IF", axis = TRUE, 
#'           clustering = clustering)
#' p <- ggPlotSim(dist(iris[, 1:4]), type = "dissimilarity", log = FALSE, 
#'                clustering = clustering, cluster_col = "blue")
#' # custom palette
#' p + scale_fill_gradient(low = "yellow", high = "red")
#' # dsCMatrix
#' m <- Matrix(c(0, 0, 2, 0, 3, 0, 2, 0, 0), ncol = 3)
#' res <- adjClust(m)
#' ggPlotSim(m, axis = TRUE)
#' ggPlotSim(m, dendro = res)
#' # dgCMatrix
#' m <- as(m, "generalMatrix")
#' ggPlotSim(m)
#' m <- as.dist(m)
#' if (require("HiTC", quietly = TRUE)) {
#'   load(system.file("extdata", "hic_imr90_40_XX.rda", package = "adjclust"))
#'   res <- hicClust(hic_imr90_40_XX, log = TRUE)
#'   ggPlotSim(hic_imr90_40_XX, axis = TRUE)
#' }
#' if (requireNamespace("snpStats", quietly = TRUE)) {
#'   data(testdata, package = "snpStats")
#'   ggPlotSim(Autosomes[1:200, 1:5], h = 3, stats = "R.squared")
#' }

#' @seealso \code{\link{select}}, \code{\link{adjClust}}
#' @export

ggPlotSim <- function(mat, type = c("similarity", "dissimilarity"),
                      clustering = NULL, dendro = NULL, log = TRUE, 
                      legendName = "intensity", main = NULL, priorCount = 0.5, 
                      stats = c("R.squared", "D.prime"), h = NULL,
                      axis = FALSE, naxis = 10, axistext = NULL, 
                      xlab = "objects", cluster_col = "darkred") {
  UseMethod("ggPlotSim")
}

#' @export
ggPlotSim.dsCMatrix <- function(mat, type = c("similarity", "dissimilarity"),
                                clustering = NULL, dendro = NULL, log = TRUE, 
                                legendName = "intensity", main = NULL, 
                                priorCount = 0.5, 
                                stats = c("R.squared", "D.prime"), h = NULL,
                                axis = FALSE, naxis = 10, axistext = NULL, 
                                xlab = "objects", cluster_col = "darkred") {
  p <- ggPlotSim.default(mat, type, clustering, dendro, log, legendName, main, 
                         priorCount, axis = axis, naxis = naxis, 
                         axistext = axistext, xlab = xlab, 
                         cluster_col = cluster_col)

  return(p)
}

#' @export
ggPlotSim.dgCMatrix <- function(mat, type = c("similarity", "dissimilarity"),
                                clustering = NULL, dendro = NULL, log = TRUE, 
                                legendName = "intensity", main = NULL, 
                                priorCount = 0.5, 
                                stats = c("R.squared", "D.prime"), h = NULL,
                                axis = FALSE, naxis = 10, axistext = NULL,
                                xlab = "objects", cluster_col = "darkred") {
  p <- ncol(mat)
  if (!isSymmetric(mat)) 
    warning(paste("Input matrix was not symmetric. Plotting only the",
                  "upper-triangular part of the matrix."))
  
  mat <- forceSymmetric(mat)
  
  p <- ggPlotSim.dsCMatrix(mat, type, clustering, dendro, log, legendName, main, 
                           priorCount, axis = axis, naxis = naxis, 
                           axistext = axistext, xlab = xlab, 
                           cluster_col = cluster_col)
  
  return(p)
}

#' @export
ggPlotSim.dist <- function(mat, type = c("similarity", "dissimilarity"),
                           clustering = NULL, dendro = NULL, log = TRUE, 
                           legendName = "intensity", main = NULL, 
                           priorCount = 0.5, stats = c("R.squared", "D.prime"),
                           h = NULL, axis = FALSE, naxis = 10, 
                           axistext = NULL, xlab = "objects", 
                           cluster_col = "darkred") {
  
  type <- match.arg(type)
  if (type != "dissimilarity") {
    message(paste("Note: input class is 'dist' so 'type' is supposed to be",
                  "'dissimilarity'."))
    type <- "dissimilarity"
  }
  
  mat <- as.matrix(mat)
  
  p <- ggPlotSim.default(mat, type, clustering, dendro, log, legendName, main, 
                         priorCount, axis = axis, naxis = naxis, 
                         axistext = axistext, xlab = xlab, 
                         cluster_col = cluster_col)
  
  return(p)
}

#' @export
ggPlotSim.HTCexp <- function(mat, type = c("similarity", "dissimilarity"),
                             clustering = NULL, dendro = NULL, log = TRUE, 
                             legendName = "IF", main = NULL, priorCount = 0.5, 
                             stats = c("R.squared", "D.prime"), h = NULL,
                             axis = FALSE, naxis = 10, axistext = NULL, 
                             xlab = "bins", cluster_col = "darkred") {
  type <- match.arg(type)
  if (!requireNamespace("HiTC")) 
    stop("Package 'HiTC' not available. 'HTCexp' input cannot be used.")
  if (type != "similarity") 
    stop("type 'dissimilarity' does not match 'HTCexp' data")
    
  mat <- mat@intdata
  p <- ggPlotSim(mat, type, clustering, dendro, log, legendName, main, 
                 priorCount, axis = axis, naxis = naxis, axistext = axistext, 
                 xlab = xlab, cluster_col = cluster_col)
  
  return(p)
}

#' @export
ggPlotSim.SnpMatrix <- function(mat, type = c("similarity", "dissimilarity"),
                                clustering = NULL, dendro = NULL, log = TRUE, 
                                legendName = "correlation", main = NULL, 
                                priorCount = 0.5, 
                                stats = c("R.squared", "D.prime"), h = NULL,
                                axis = FALSE, naxis = 10, axistext = NULL,
                                xlab = "SNP index", cluster_col = "darkred") {
  if (!requireNamespace("snpStats")) 
    stop("Package 'snpStats' not available. 'SnpMatrix' input cannot be used.")
  
  stats <- match.arg(stats)
  if (is.null(h)) h <- ncol(mat) - 1
  if (h >= ncol(mat)) stop("h should be strictly less than p")
  
  mat <- snpStats::ld(mat, stats = stats, depth = h)
  mat[mat > 1] <- 1  ## fix numerical aberrations
  mat[mat < 0] <- 0  ## fix numerical aberrations
  diag(mat) <- rep(1, nrow(mat))  ## by default the diagonal is 0 after 'snpStats::ld'
  
  p <- ggPlotSim(mat, type, clustering, dendro, log, legendName, main, 
                 priorCount, axis, naxis = naxis, axistext = axistext, 
                 xlab = xlab, cluster_col = cluster_col)
  
  return(p)
}
 
#' @export
ggPlotSim.default <- function(mat, type = c("similarity", "dissimilarity"),
                              clustering = NULL, dendro = NULL, log = TRUE, 
                              legendName = "intensity", main = NULL, 
                              priorCount = 0.5, 
                              stats = c("R.squared", "D.prime"), h = NULL,
                              axis = FALSE, naxis = 10, axistext = NULL,
                              xlab = "objects", cluster_col = "darkred") {
  # Input checks ####
  d <- nrow(mat)
  type <- match.arg(type)
  if (!is.null(clustering)) {
    clusters <- unique(clustering)
    if ((length(clustering) != d) || (sort(clusters) != 1:max(clusters))) {
      stop(paste("'clustering' must be a vector of numeric values between 1",
                 "and K (nb of clusters). Please fix input."))
    }
  }
  if (!is.null(dendro)) {
    if (inherits(dendro, "hclust")) {
      dd <- dendro
    } else dd <- try(as.hclust(dendro), silent = TRUE)
    if (inherits(dd, "try-error")) {
      stop(paste("'dendro' can not be converted to class 'hclust'. Please,",
                 "provide a proper dendrogram."))
    }
  }
  if (!is.logical(log)) stop("'log' must be logical!")
  if (!is.character(legendName)) stop("'legendName' must be a string!")
  if (!is.null(main) && !is.character(main)) stop("'main' must be a string!")
  if (log) {
    if (!is.numeric(priorCount) || length(priorCount) > 1 || priorCount < 0) {
      stop(paste("'priorCount' must be a single non-negative number!"))
    }
  }
  if (!is.logical(axis)) stop("'axis' must be logical!")
  if (length(naxis) > 1 || !is.integer(naxis)) {
    stop(paste("'naxis' must be a single value of type integer!"))
  }
  if (!is.null(axistext)) {
    if (length(axistext) != naxis) {
      stop("'axistext' length must be equal to 'naxis'.")
    }
    if (!is.character(axistext)) axistext <- as.character(axistext)
  }
  if (!is.character(xlab)) stop("'xlab' must be a string!")

  # Coordinate computation ####
  if (type == "dissimilarity") mat <- max(mat) - mat
  
  coordinates <- poly_coords(mat)
  fake_coords <- make_coords(c(1, d, d), c(1, d, 1), rep(0, 3))
  
  if (!is.null(dendro)) {
    ymax_i <- max(fake_coords$y)
    ymax_f <- max(dd$height)
    xmin_i <- min(fake_coords$x)
    xmax_i <- max(fake_coords$x)
    coordinates[, c("x", "y")] <- rescale_coords(coordinates$x, coordinates$y,
                                                 ymax_i, ymax_f, xmin_i, xmax_i,
                                                 d)
    fake_coords[, c("x", "y")] <- rescale_coords(fake_coords$x, fake_coords$y,
                                                 ymax_i, ymax_f, xmin_i, xmax_i,
                                                 d)

    dd <- as.dendrogram(dd)
    dd <- dd %>% set("labels", value = rep(NA, d)) %>% as.ggdend()
    p <- ggplot(dd) + theme_void() +
      geom_polygon(data = fake_coords, aes(x = x, y = y), fill = "lightgrey")
  } else {
    p <- ggplot() + theme_void() +
      geom_polygon(data = fake_coords, aes(x = x, y = y), fill = "lightgrey")
  }
  
  # Plot ####
  if (log) {
    if (legendName != "") {
      legendName <- ifelse(priorCount == 0,
                           paste0("log(", legendName, ")"),
                           paste0("log(", legendName, " + ", priorCount, ")"))
    }
    p <- p + geom_polygon(data = coordinates, 
                          aes(x = x, y = y, group = id, 
                              fill = log(IF + priorCount))) + 
      scale_fill_viridis_b(name = legendName)
  } else {
    p <- p + geom_polygon(data = coordinates, 
                          aes(x = x, y = y, group = id, fill = IF)) + 
      scale_fill_viridis_b(name = legendName)
  }
  
  # Additional option ####
  if (!is.null(clustering)) {
    all_clusters <- unique(clustering)
    nb_clust <- length(all_clusters)
    start_clusters <- sapply(all_clusters, 
                             function(cc) min(which(clustering == cc)))
    end_clusters <- sapply(all_clusters, 
                           function(cc) max(which(clustering == cc)))
    indi <- c(start_clusters, end_clusters, end_clusters)
    indj <- c(start_clusters, start_clusters, end_clusters)
    cluster_coords <- make_coords(indi, indj, rep(0, length(indi)))
    cluster_coords <- cluster_coords[order(cluster_coords$id), ]
    to_keep <- c(3*(1:nb_clust), 
                 3*nb_clust + 4*(1:nb_clust) - 1,
                 7*nb_clust + 3*(1:nb_clust) - 2)
    cluster_coords <- cluster_coords[to_keep, ]
    if (!is.null(dendro)) {
      cluster_coords[, c("x", "y")] <- rescale_coords(cluster_coords$x, 
                                                      cluster_coords$y,
                                                      ymax_i, ymax_f, xmin_i, 
                                                      xmax_i, d)
    }
    cluster_coords$cluster <- rep(1:nb_clust, 3)
    p <- p + geom_path(data = cluster_coords, 
                       aes(x = x, y = y, group = cluster),
                       size = 1, colour = cluster_col)
  }
  
  if (!is.null(main)) p <- p + ggtitle(main)
  
  if (axis) {
    displayed_bins <- floor(seq(1, d, length.out = naxis))
    displayed_x <- make_coords(displayed_bins, displayed_bins, rep(0, naxis))
    displayed_x <- displayed_x$x[(naxis + 1):(2 * naxis)]
    if (!is.null(dendro)) {
      # note: y has no meaning here!
      displayed_x <- rescale_coords(displayed_x, displayed_x, ymax_i, ymax_f, 
                                    xmin_i, xmax_i, d)$x
    }
    if (is.null(axistext)) axistext <- displayed_bins
    p <- p + theme(axis.title.x = element_text(), axis.text.x = element_text(),
                   axis.ticks.x = element_line(), 
                   axis.ticks.length.x = unit(0.25, "cm")) +
      scale_x_continuous(name = xlab, breaks = displayed_x, labels = axistext)
    if (!is.null(dendro)) {
      p <- p + theme(axis.line.x = element_line(), 
                     panel.grid.major.x = element_line(colour = "lightgrey"))
    }
  }
  
  return(p)
}
