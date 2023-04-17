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
#' @param clustering vector of clusters to display on the matrix (if not 
#' \code{NULL}). If \code{clustering} is provided, it must be a numeric vector
#' of length identical to the matrix size with clusters identified as 
#' consecutive integers 1, 2, 3, ...
#' @param dendro dendrogram provided as an \code{hclust} object, or as another
#' type of object that can be converted to \code{hclust} with \code{as.hclust}.
#' @param k number of clusters to display. Used only when \code{dendro} is not
#' \code{NULL} and \code{clustering} is \code{NULL}. The clustering is then 
#' deduced from the dendrogram by a standard cut.
#' @param log logical. Should the breaks be based on log-scaled values of the
#' matrix entries. Default to \code{TRUE}.
#' @param legendName character. Title of the legend. Default to 
#' \code{"intensity"}.
#' @param main character. Title of the plot. Default to \code{NULL} (no title).
#' @param priorCount numeric. Average count to be added to each entry of the
#' matrix to avoid taking log of zero. Used only if \code{log == TRUE} and 
#' default to 0.5.
#' @param stats input SNP correlation type. Used when \code{mat} is of type 
#' \code{SnpMatrix}.
#' @param h positive integer. Threshold distance for SNP correlation 
#' computation. Used when \code{mat} is of type \code{SnpMatrix}. 
#' @param axis logical. Should x-axis be displayed on the plot? Default to 
#' \code{FALSE}.
#' @param naxis integer. If \code{axis == TRUE}, number of ticks to display on
#' the x-axis.
#' @param axistext character vector. If \code{axis == TRUE}, labels to display
#' of the x-axis (its length has to be equal to \code{naxis}).
#' @param xlab character. If \code{axis == TRUE}, x-axis title.
#' @param cluster_col colour for the cluster line if \code{clustering} is not
#' \code{NULL}.
#' @param mode type of dendrogram to plot (see \code{\link{plot.chac}}). Default
#' to \code{"standard"}.
#' @import ggplot2
#' @importFrom dendextend set as.ggdend %>%
#' @importFrom rlang .data
#' @examples \dontrun{
#' clustering <- rep(1:3, each = 50)
#' dist_data <- as.matrix(dist(iris[, 1:4]))
#' dendro_iris <- adjClust(dist_data, type = "dissimilarity")
#' plotSim(dist_data, type = "dissimilarity", dendro = dendro_iris, axis = TRUE)
#' plotSim(dist_data, type = "dissimilarity", dendro = dendro_iris,
#'         clustering = clustering)
#' plotSim(dist_data, type = "dissimilarity", dendro = dendro_iris, axis = TRUE,
#'         k = 3)
#' plotSim(dist_data, type = "dissimilarity", legendName = "IF", axis = TRUE, 
#'         clustering = clustering)
#' p <- plotSim(dist(iris[, 1:4]), type = "dissimilarity", log = FALSE, 
#'              clustering = clustering, cluster_col = "blue")
#' # custom palette
#' p + scale_fill_gradient(low = "yellow", high = "red")
#' # dsCMatrix
#' m <- Matrix(c(0, 0, 2, 0, 3, 0, 2, 0, 0), ncol = 3)
#' res <- adjClust(m)
#' plotSim(m, axis = TRUE)
#' plotSim(m, dendro = res)
#' # dgCMatrix
#' m <- as(m, "generalMatrix")
#' plotSim(m)
#' m <- as.dist(m)
#' if (require("HiTC", quietly = TRUE)) {
#'   load(system.file("extdata", "hic_imr90_40_XX.rda", package = "adjclust"))
#'   res <- hicClust(hic_imr90_40_XX, log = TRUE)
#'   plotSim(hic_imr90_40_XX, axis = TRUE)
#' }
#' if (requireNamespace("snpStats", quietly = TRUE)) {
#'   data(testdata, package = "snpStats")
#'   plotSim(Autosomes[1:200, 1:5], h = 3, stats = "R.squared", axis = TRUE,
#'           axistext = c("A", "B", "C", "D", "E"))
#' }
#' }

#' @seealso \code{\link{select}}, \code{\link{adjClust}}
#' @export

plotSim <- function(mat, type = c("similarity", "dissimilarity"),
                    clustering = NULL, dendro = NULL, k = NULL, log = TRUE, 
                    legendName = "intensity", main = NULL, priorCount = 0.5, 
                    stats = c("R.squared", "D.prime"), h = NULL, axis = FALSE,
                    naxis = min(10, nrow(mat)), axistext = NULL, 
                    xlab = "objects", cluster_col = "darkred",
                    mode = c("standard", "corrected", "total-disp", 
                             "within-disp", "average-disp")) {
  UseMethod("plotSim")
}

#' @export
plotSim.dsCMatrix <- function(mat, type = c("similarity", "dissimilarity"),
                              clustering = NULL, dendro = NULL, k = NULL, 
                              log = TRUE, legendName = "intensity", main = NULL, 
                              priorCount = 0.5, 
                              stats = c("R.squared", "D.prime"), h = NULL,
                              axis = FALSE, naxis = min(10, nrow(mat)), 
                              axistext = NULL,xlab = "objects", 
                              cluster_col = "darkred",
                              mode = c("standard", "corrected", "total-disp", 
                                       "within-disp", "average-disp")) {
  p <- plotSim.default(mat, type, clustering, dendro, k, log, legendName, main, 
                       priorCount, axis = axis, naxis = naxis, 
                       axistext = axistext, xlab = xlab, 
                       cluster_col = cluster_col, mode = mode)

  return(p)
}

#' @export
plotSim.dgCMatrix <- function(mat, type = c("similarity", "dissimilarity"),
                              clustering = NULL, dendro = NULL, k = NULL, 
                              log = TRUE, legendName = "intensity", 
                              main = NULL, priorCount = 0.5, 
                              stats = c("R.squared", "D.prime"), h = NULL,
                              axis = FALSE, naxis = min(10, nrow(mat)), 
                              axistext = NULL, xlab = "objects", 
                              cluster_col = "darkred",
                              mode = c("standard", "corrected", "total-disp", 
                                       "within-disp", "average-disp")) {
  if (!isSymmetric(mat)) 
    warning(paste("Input matrix was not symmetric. Plotting only the",
                  "upper-triangular part of the matrix."))
  
  mat <- forceSymmetric(mat)
  
  p <- plotSim.default(mat, type, clustering, dendro, k, log, legendName, main, 
                       priorCount, axis = axis, naxis = naxis, 
                       axistext = axistext, xlab = xlab, 
                       cluster_col = cluster_col, mode = mode)
  
  return(p)
}

#' @export
plotSim.dist <- function(mat, type = c("similarity", "dissimilarity"),
                         clustering = NULL, dendro = NULL, k= NULL, log = TRUE, 
                         legendName = "intensity", main = NULL, 
                         priorCount = 0.5, stats = c("R.squared", "D.prime"),
                         h = NULL, axis = FALSE, naxis = min(10, nrow(mat)), 
                         axistext = NULL, xlab = "objects", 
                         cluster_col = "darkred",
                         mode = c("standard", "corrected", "total-disp", 
                                  "within-disp", "average-disp")) {
  
  type <- match.arg(type)
  if (type != "dissimilarity") {
    message(paste("Note: input class is 'dist' so 'type' is supposed to be",
                  "'dissimilarity'."))
    type <- "dissimilarity"
  }
  
  mat <- as.matrix(mat)
  
  p <- plotSim.default(mat, type, clustering, dendro, k, log, legendName, main, 
                       priorCount, axis = axis, naxis = naxis, 
                       axistext = axistext, xlab = xlab, 
                       cluster_col = cluster_col, mode = mode)
  
  return(p)
}

#' @export
plotSim.HTCexp <- function(mat, type = c("similarity", "dissimilarity"),
                           clustering = NULL, dendro = NULL, k = NULL, 
                           log = TRUE, legendName = "IF", main = NULL, 
                           priorCount = 0.5, stats = c("R.squared", "D.prime"), 
                           h = NULL, axis = FALSE, naxis = min(10, nrow(mat)),
                           axistext = NULL, xlab = "bins", 
                           cluster_col = "darkred",
                           mode = c("standard", "corrected", "total-disp", 
                                    "within-disp", "average-disp")) {
  type <- match.arg(type)
  if (!requireNamespace("HiTC")) 
    stop("Package 'HiTC' not available. 'HTCexp' input cannot be used.")
  if (type != "similarity") 
    stop("type 'dissimilarity' does not match 'HTCexp' data")
    
  mat <- mat@intdata
  p <- plotSim(mat, type, clustering, dendro, k, log, legendName, main, 
               priorCount, axis = axis, naxis = naxis, axistext = axistext, 
               xlab = xlab, cluster_col = cluster_col, mode = mode)
  
  return(p)
}

#' @export
plotSim.SnpMatrix <- function(mat, type = c("similarity", "dissimilarity"),
                              clustering = NULL, dendro = NULL, k = NULL,
                              log = TRUE, legendName = "correlation", 
                              main = NULL, priorCount = 0.5, 
                              stats = c("R.squared", "D.prime"), h = NULL,
                              axis = FALSE, naxis = min(10, nrow(mat)), 
                              axistext = NULL, xlab = "SNP index", 
                              cluster_col = "darkred",
                              mode = c("standard", "corrected", "total-disp", 
                                       "within-disp", "average-disp")) {
  if (!requireNamespace("snpStats")) 
    stop("Package 'snpStats' not available. 'SnpMatrix' input cannot be used.")
  
  stats <- match.arg(stats)
  if (is.null(h)) h <- ncol(mat) - 1
  if (!is.numeric(h) || h <= 0 || h >= ncol(mat)) 
    stop("'h' should be numeric, larger than 0 and smaller than p.")
  
  mat <- snpStats::ld(mat, stats = stats, depth = h)
  mat[mat > 1] <- 1  ## fix numerical aberrations
  mat[mat < 0] <- 0  ## fix numerical aberrations
  diag(mat) <- rep(1, nrow(mat))  ## by default the diagonal is 0 after 'snpStats::ld'
  mat <- forceSymmetric(mat)
  
  p <- plotSim.default(mat, type, clustering, dendro, k, log, legendName, main, 
                       priorCount, h = h, axis = axis, naxis = naxis, 
                       axistext = axistext, xlab = xlab, 
                       cluster_col = cluster_col, mode = mode)
  
  return(p)
}
 
#' @export
plotSim.default <- function(mat, type = c("similarity", "dissimilarity"),
                            clustering = NULL, dendro = NULL, k = NULL, 
                            log = TRUE, legendName = "intensity", main = NULL, 
                            priorCount = 0.5, stats = c("R.squared", "D.prime"), 
                            h = NULL, axis = FALSE, naxis = min(10, nrow(mat)),
                            axistext = NULL, xlab = "objects", 
                            cluster_col = "darkred", 
                            mode = c("standard", "corrected", "total-disp", 
                                     "within-disp", "average-disp")) {
  # Input checks ####
  d <- nrow(mat)
  type <- match.arg(type)
  if (!is.null(dendro)) {
    mode <- match.arg(mode)
    if (inherits(dendro, "hclust")) {
      dd <- dendro
    } else dd <- try(as.hclust(dendro), silent = TRUE)
    dendro <- dendro_for_mode(dendro, mode)
    if (inherits(dd, "try-error")) {
      stop(paste("'dendro' can not be converted to class 'hclust'. Please ",
                 "provide a proper dendrogram."))
    }
    if (length(dd$order) != d)
      stop(paste("'dendro' must be a dendrogram with the same number of",
                 "objects as the input matrix. Please fix input."))
  }
  if (!is.null(k)) {
    if (!is.null(clustering))
      stop("Either 'clustering' or 'k' can be not NULL, not both.")
    
    if (is.null(dendro))
      stop("'dendro' must be provided when 'k' is not NULL.")
    
    if (k != round(k, 0) | k <= 1 | k >= length(dd$order))
      stop(paste("'k' must be an integer larger than 1 and smaller than the",
                 "number of objects in the input matrix. Please fix input."))
    clustering <- cutree(dd, k = k)
  }
  if (!is.null(clustering)) {
    clusters <- sort(unique(clustering))
    if ((length(clustering) != d) || !all.equal(clusters, 1:max(clusters))) {
      stop(paste("'clustering' must be a vector of size nrow(mat) containing",
                 "all the integers between 1 and K (nb of clusters).",
                 "Please fix input."))
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
  if (axis) {
    if (length(naxis) > 1 || round(naxis) != naxis) {
      stop(paste("'naxis' must be a single value of type integer!"))
    }
    if (naxis > d) {
      warning(paste("Reducing the number of ticks on x-axis to the number of",
                    "objects."))
      naxis <- d
    }
    if (!is.null(axistext)) {
      if (length(axistext) != naxis) {
        stop("'axistext' length must be equal to 'naxis'.")
      }
      if (!is.character(axistext)) axistext <- as.character(axistext)
    }
    if (!is.character(xlab)) stop("'xlab' must be a string!")
  }

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
  }
  if (!is.null(h) && h < d - 1) {
    fake_coords <- make_coords(c(1, d, d, h+1), c(1, d, d-h, 1), rep(0, 4))
  }
  if (!is.null(dendro)) {
    fake_coords[, c("x", "y")] <- rescale_coords(fake_coords$x, fake_coords$y,
                                                 ymax_i, ymax_f, xmin_i, xmax_i,
                                                 d)

    dd <- as.dendrogram(dd)
    dd <- dd %>% set("labels", value = rep(NA, d)) %>% as.ggdend()
    p <- ggplot(dd) + theme_void() +
      geom_polygon(data = fake_coords, aes(x = .data$x, y = .data$y), 
                   fill = "lightgrey")
  } else {
    p <- ggplot() + theme_void() +
      geom_polygon(data = fake_coords, aes(x = .data$x, y = .data$y), 
                   fill = "lightgrey")
  }
  
  # Plot ####
  if (log) {
    if (legendName != "") {
      legendName <- ifelse(priorCount == 0,
                           paste0("log(", legendName, ")"),
                           paste0("log(", legendName, " + ", priorCount, ")"))
    }
    p <- p + geom_polygon(data = coordinates, 
                          aes(x = .data$x, y = .data$y, group = .data$id, 
                              fill = log(.data$IF + priorCount))) + 
      scale_fill_viridis_b(name = legendName)
  } else {
    p <- p + geom_polygon(data = coordinates, 
                          aes(x = .data$x, y = .data$y, group = .data$id, 
                              fill = .data$IF)) + 
      scale_fill_viridis_b(name = legendName)
  }
  
  # Additional options ####
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
                       aes(x = .data$x, y = .data$y, group = .data$cluster),
                       linewidth = 1, colour = cluster_col)
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

# Auxiliary functions ####

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