#' Adjacency-constrained Clustering of Hi-C contact maps
#' 
#' Adjacency-constrained hierarchical agglomerative clustering of Hi-C contact
#' maps
#' 
#' Adjacency-constrained hierarchical agglomerative clustering (HAC) is HAC in
#' which each observation is associated to a position, and the clustering is 
#' constrained so as only adjacent clusters are merged. Genomic regions (loci)
#' are clustered according to information provided by high-throughput
#' conformation capture data (Hi-C).
#' 
#' @param x either: 1. A pxp contact sparse or dense matrix (classes matrix,
#' Matrix, dscMatrix, dgTMatrix, dgCMatrix, dgeMatrix). Its entries are the 
#' number of counts of physical interactions observed between all pairs of loci. 
#' 2. An object of class HiTC::HTCexp. The corresponding Hi-C data is stored as 
#' a Matrix::dsCMatrix object in the intdata slot. 3. A text file path with one 
#' line per pair of loci for which an interaction has been observed (in the 
#' format: locus1<tab>locus2<tab>signal) or a matrix or data frame with similar
#' data (3 columns).
#'   
#' @param h band width. If not provided, \code{h} is set to default value `p-1`.
#' 
#' @param log logical. Whether to log-transform the count data. Default to 
#' \code{FALSE}.
#'   
#' @param \dots further arguments to be passed to \code{\link{read.table}} 
#'   function when \code{x} is a text file name. If not provided, the text file 
#'   is supposed to be separated by tabulations, with no header.
#'   
#' @return An object of class \code{\link{chac}}.
#'   
#' @seealso \code{\link{adjClust}}
#'   
#' @references Ambroise C., Dehman A., Neuvial P., Rigaill G., and Vialaneix N
#'   (2019). \emph{Adjacency-constrained hierarchical clustering of a band
#'   similarity matrix with application to genomics}, Algorithms for Molecular
#'   Biology 14(22)"
#'
#' @references Servant N. \emph{et al} (2012). \emph{HiTC : Exploration of 
#'   High-Throughput 'C' experiments. Bioinformatics}.
#'   
#'   
#' @examples
#' # input as HiTC::HTCexp object
#' \dontrun{
#' if (require("HiTC", quietly = TRUE)) {
#'   load(system.file("extdata", "hic_imr90_40_XX.rda", package = "adjclust"))
#'   res1 <- hicClust(hic_imr90_40_XX)
#' }
#' }
#' 
#' # input as Matrix::dsCMatrix contact map
#' \dontrun{
#' mat <- HiTC::intdata(hic_imr90_40_XX) 
#' res2 <- hicClust(mat)
#' }
#' 
#' # input as text file
#' \dontrun{
#' res3 <- hicClust(system.file("extdata", "sample.txt", package = "adjclust"))
#' }
#' 
#' @export
#' 
#' @importFrom utils read.table
#' @importFrom methods new

hicClust <- function(x, h = NULL, log = FALSE, ...) {
  UseMethod("hicClust")
}

#' @export
hicClust.data.frame <- function(x, h = NULL, log = FALSE, ...) { # bin pair list
  lis <- sort(unique(c(x[,1], x[,2])))
  p <- length(lis)
  rowindx <- match(x[,1], lis)
  colindx <- match(x[,2], lis)
  identical <- rowindx == colindx
  x[ ,3] <- as.numeric(x[ ,3])
  
  mat <- new("dgTMatrix", Dim = c(as.integer(p), as.integer(p)),
             i = c(rowindx, colindx[!identical]) - 1L, 
             j = c(colindx, rowindx[!identical]) - 1L, 
             x = c(x[ ,3], x[!identical,3]))
  if (log) mat@x <- log(mat@x + 1)
  
  res <- run.hicclust(mat, h = h)
  x <- sys.call()
  res$call <- update_call(x, "hicClust")
  return(res)
}

#' @export
hicClust.character <- function(x, h = NULL, log = FALSE, ...) { # file with bin pair list
  if (!file.exists(x)) {
    stop("Input of type 'character' should be a valid file.")
  }
  inoptions <- list(...)
  inoptions$file <- x
  if (is.null(inoptions$sep)) {
    inoptions$sep <- "\t"
  }
  if (is.null(inoptions$header)) {
    inoptions$header <- FALSE
  }
  if (is.null(inoptions$stringsAsFactors)) {
    inoptions$stringsAsFactors <- FALSE
  }
  df <- do.call("read.table", inoptions) 
  res <- hicClust.data.frame(df, h = h, log = log)
  x <- sys.call()
  res$call <- update_call(x, "hicClust")
  
  return(res)
}

#' @export
hicClust.matrix <- function(x, h = NULL, log = FALSE, ...) {
  if ((nrow(x) != ncol(x)) & (ncol(x) == 3)) {
    # bin pair list
    x <- as.data.frame(x)
    res <- hicClust.data.frame(x, h = h, log = log)
  } else {
    # full hic matrix
    if (log) x <- log(x + 1)
    res <- run.hicclust(x, h = h)
  }
  x <- sys.call()
  res$call <- update_call(x, "hicClust")
  return(res)
}

#' @export
hicClust.Matrix <- function(x, h = NULL, log = FALSE, ...) {
  if (log) x <- log(x + 1)
  res <- run.hicclust(x, h = h)
  x <- sys.call()
  res$call <- update_call(x, "hicClust")
  return(res)
}

#' @export
hicClust.dgCMatrix <- function(x, h = NULL, log = FALSE, ...) {
  if (log) x@x <- log(x@x + 1)
  res <- run.hicclust(x, h = h)
  x <- sys.call()
  res$call <- update_call(x, "hicClust")
  return(res)
}

#' @export
hicClust.dsCMatrix <- function(x, h = NULL, log = FALSE, ...) {
  if (log) x@x <- log(x@x + 1)
  res <- run.hicclust(x, h = h)
  x <- sys.call()
  res$call <- update_call(x, "hicClust")
  return(res)
}

#' @export
hicClust.dgeMatrix <- function(x, h = NULL, log = FALSE, ...) {
  if (log) x <- log(x + 1)
  res <- run.hicclust(x, h = h)
  x <- sys.call()
  res$call <- update_call(x, "hicClust")
  return(res)
}

#' @export
hicClust.HTCexp <- function(x, h = NULL, log = FALSE, ...) {
  if (!requireNamespace("HiTC", quietly = TRUE))
    stop("Package 'HiTC' not available. This function cannot be used with 'HTCexp' data.") # nocov
  x <- HiTC::intdata(x)
  res <- hicClust(x, h = h, log = log) # sparse or dense version
  x <- sys.call()
  res$call <- update_call(x, "hicClust")
  return(res)
}

run.hicclust <- function(x, h) {
  if (is.null(h)) h <- nrow(x) - 1
  if (!is.numeric(h))
    stop("h should be numeric")
  
  res <- adjClust(x, type = "similarity", h = h)
  res$method <- "hicClust"
  x <- sys.call()
  res$call <- update_call(x, "hicClust")
  return(res)
}
