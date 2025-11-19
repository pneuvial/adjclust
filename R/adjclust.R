#' @useDynLib adjclust, .registration = TRUE
NULL

#' Adjacency-constrained Clustering
#'
#' Adjacency-constrained hierarchical agglomerative clustering
#'
#' Adjacency-constrained hierarchical agglomerative clustering (HAC) is HAC in
#' which each observation is associated to a position, and the clustering is
#' constrained so as only adjacent clusters are merged. These methods are useful
#' in various application fields, including ecology (Quaternary data) and
#' bioinformatics (e.g., in Genome-Wide Association Studies (GWAS)).
#'
#' This function is a fast implementation of the method that takes advantage of
#' sparse similarity matrices (i.e., that have 0 entries outside of a diagonal
#' band of width \code{h}). The method is fully described in (Dehman, 2015) and
#' based on a kernel version of the algorithm. The different options for the
#' implementation are available in the package vignette entitled 
#' \href{https://pneuvial.github.io/adjclust/articles/notesCHAC.html}{"Notes on CHAC implementation in adjclust"}.
#'
#' @param mat A similarity matrix or a dist object. Most sparse formats from
#'   \code{\link[Matrix]{sparseMatrix}} are allowed
#' @param type Type of matrix : similarity or dissimilarity. Defaults to
#'   \code{"similarity"}
#' @param h band width. It is assumed that the similarity between two items is 0
#'   when these items are at a distance of more than band width h. Default value
#'   is \code{ncol(mat)-1}
#' @param strictCheck Logical (default to \code{TRUE}) to systematically check
#'   default of positivity in input similarities. Can be disabled to avoid
#'   computationally expensive checks when the number of features is large.
#' @param nthreads Integer (default to \code{1L}). Number of threads use for
#' matrix pre-computations.
#'
#' @returns An object of class \code{\link{chac}} which describes the tree
#'   produced by the clustering process. The object is a list with the same
#'   elements as an object of class \code{\link[stats]{hclust}} (\code{merge},
#'   \code{height}, \code{order}, \code{labels}, \code{call}, \code{method},
#'   \code{dist.method}), and two extra elements: 
#'     \item{\code{mat}}{: (the data on which the clustering has been performed, 
#'     possibly after the pre-transformations described in the vignette entitled
#'     \href{https://pneuvial.github.io/adjclust/articles/notesCHAC.html#notes-on-relations-between-similarity-and-dissimilarity-implementation}{"Notes on CHAC implementation in adjclust"}}.
#'     \item{\code{correction}}{: the value of the correction for non positive
#'     definite similarity matrices (also described in the same vignette). If
#'     \code{correction == 0}, it means that the initial data were not 
#'     pre-transformed.}
#'
#' @seealso \code{\link{snpClust}} to cluster SNPs based on linkage
#'   disequilibrium
#' @seealso \code{\link{hicClust}} to cluster Hi-C data
#'
#' @references Murtagh F., and Legendre P. (2014). Ward's hierarchical
#'   agglomerative clustering method: which algorithms implement Ward's 
#'   criterion? \emph{Journal of Classification}, \strong{31}, 274-295. 
#'   DOI: \doi{10.1007/s00357-014-9161-z}.
#' 
#' @references Dehman A. (2015). \emph{Spatial Clustering of Linkage
#'   Disequilibrium Blocks for Genome-Wide Association Studies}, PhD thesis,
#'   Universite Paris Saclay, France.
#'
#' @references Ambroise C., Dehman A., Neuvial P., Rigaill G., and Vialaneix N
#'   (2019). Adjacency-constrained hierarchical clustering of a band similarity 
#'   matrix with application to genomics. \emph{Algorithms for Molecular
#'   Biology}, \strong{14}(22).
#'   DOI: \doi{10.1007/s11222-018-9806-6}.
#'
#' @references Randriamihamison N., Vialaneix N., and Neuvial P. (2020).
#'   Applicability and interpretability of Ward's hierarchical agglomerative
#'   clustering with or without contiguity constraints. \emph{Journal of 
#'   Classification}, \strong{38}, 1-27.
#'   DOI: \doi{10.1007/s00357-020-09377-y}.
#' 
#' @note When performed on a distance matrix \eqn{d} with the option 
#' \code{type = "dissimilarity"}, \code{adjclust} is identical to using the 
#' option \code{"ward.D"} on \eqn{d^2} in the function 
#' \code{\link[stats]{hclust}} when the ordering of the (unconstrained) 
#' clustering (in \code{\link[stats]{hclust}}) is compatible with the natural 
#' ordering of objects used as a constraint. It is also equivalent (under the
#' same assumption or orderings) to the option \code{"ward.D2"} performed on the
#' distance matrix \eqn{d} itself, except for the final heights of the merges 
#' that are equal to the square of the heights obtained with \code{"ward.D2"} in
#' \code{\link[stats]{hclust}}. See the 
#' \href{https://pneuvial.github.io/adjclust/articles/notesCHAC.html#relations-with-hclust-and-rioja}{vignette on implementation} 
#' and (Murtagh and Legendre, 2014) for further details.
#'
#' @examples
#' sim <- matrix(
#'   c(1.0, 0.1, 0.2, 0.3,
#'     0.1, 1.0 ,0.4 ,0.5,
#'     0.2, 0.4, 1.0, 0.6,
#'     0.3, 0.5, 0.6, 1.0), 
#'  nrow = 4)
#'
#' ## similarity, full width
#' fit1 <- adjClust(sim, "similarity")
#' plot(fit1)
#'
#' ## similarity, h < p-1
#' fit2 <- adjClust(sim, "similarity", h = 2)
#' plot(fit2)
#'
#' ## dissimilarity
#' dist <- as.dist(sqrt(2-(2*sim)))
#'
#' ## dissimilarity, full width
#' fit3 <- adjClust(dist, "dissimilarity")
#' plot(fit3)
#'
#' ## dissimilarity, h < p-1
#' fit4 <- adjClust(dist, "dissimilarity", h = 2)
#' plot(fit4)
#'
#' @export
#'
#' @importFrom sparseMatrixStats rowCumsums colCumsums
#' @importFrom Matrix diag
#' @importFrom Matrix t
adjClust <- function(mat, type = c("similarity", "dissimilarity"), 
                     h = ncol(mat) - 1, strictCheck = TRUE, nthreads = 1L) {
  UseMethod("adjClust")
}

#' @importFrom Matrix isSymmetric forceSymmetric
#' @export
adjClust.matrix <- function(mat, type = c("similarity", "dissimilarity"), 
                            h = ncol(mat) - 1, strictCheck = TRUE,
                            nthreads = 1L) {
  if (!is.numeric(mat))
    stop("Input matrix is not numeric")
  if (!(isSymmetric(mat)))
    stop("Input matrix is not symmetric")
  res <- run.adjclust(mat, type = type, h = h, strictCheck = strictCheck,
                      nthreads = nthreads)
  x <- sys.call()
  res$call <- update_call(x, "adjClust")
  return(res)
}

#' @export
adjClust.dsyMatrix <- function(mat, type = c("similarity", "dissimilarity"), 
                            h = ncol(mat) - 1, strictCheck = TRUE, 
                            nthreads = 1L) {
  # RcppArmadillo functions don't support dsyMatrix, so convert to matrix
  res <- run.adjclust(as.matrix(mat), type = type, h = h, 
                      strictCheck = strictCheck, nthreads = nthreads)
  x <- sys.call()
  res$call <- update_call(x, "adjClust")
  return(res)
}

#' @export
adjClust.dgeMatrix <- function(mat, type = c("similarity", "dissimilarity"), 
                               h = ncol(mat) - 1, strictCheck = TRUE, 
                               nthreads = 1L) {
  type <- match.arg(type)
  if (!(isSymmetric(mat))) {
    stop("Input matrix is not symmetric")
  } else mat <- forceSymmetric(mat)
  res <- adjClust(mat, type = type, h = h, strictCheck = strictCheck, 
                  nthreads = nthreads)
  x <- sys.call()
  res$call <- update_call(x, "adjClust")
  return(res)
}

#' @export
adjClust.dsCMatrix <- function(mat, type = c("similarity", "dissimilarity"), 
                               h = ncol(mat) - 1, strictCheck = TRUE, 
                               nthreads = 1L) {
  type <- match.arg(type)
  if (type == "dissimilarity")
    stop("'type' can only be 'similarity' with sparse Matrix inputs")
  res <- run.adjclust(mat, type = type, h = h, strictCheck = strictCheck,
                      nthreads = nthreads)
  x <- sys.call()
  res$call <- update_call(x, "adjClust")
  return(res)
}

#' @export
adjClust.dgCMatrix <- function(mat, type = c("similarity", "dissimilarity"), 
                               h = ncol(mat) - 1, strictCheck = TRUE,
                               nthreads = 1L) {
  if (!(isSymmetric(mat))) {
    stop("Input matrix is not symmetric")
  } else {
    mat <- forceSymmetric(mat)
  }
  res <- adjClust(mat, type = type, h = h, strictCheck = strictCheck, 
                  nthreads = nthreads)
  x <- sys.call()
  res$call <- update_call(x, "adjClust")
  return(res)
}

#' @export
adjClust.dsTMatrix <- function(mat, type = c("similarity", "dissimilarity"), 
                               h = ncol(mat) - 1, strictCheck = TRUE, 
                               nthreads = 1L) {
  type <- match.arg(type)
  if (type == "dissimilarity")
    stop("'type' can only be 'similarity' with sparse Matrix inputs")
  res <- run.adjclust(mat, type = type, h = h, strictCheck = strictCheck, 
                      nthreads = nthreads)
  x <- sys.call()
  res$call <- update_call(x, "adjClust")
  return(res)
}

#' @export
adjClust.dgTMatrix <- function(mat, type = c("similarity", "dissimilarity"), 
                               h = ncol(mat) - 1, strictCheck = TRUE,
                               nthreads = 1L) {
  type <- match.arg(type)
  if (!(isSymmetric(mat))) {
    stop("Input matrix is not symmetric")
  } else {
    mat <- forceSymmetric(mat)
  }
  res <- adjClust(mat, type = type, h = h, strictCheck = strictCheck,
                  nthreads = nthreads)
  x <- sys.call()
  res$call <- update_call(x, "adjClust")
  return(res)
}

#' @export
adjClust.dist <- function(mat, type = c("similarity", "dissimilarity"), 
                          h = ncol(mat) - 1, strictCheck = TRUE, 
                          nthreads = 1L) {
  type <- match.arg(type)
  if (type != "dissimilarity")
    message("Note: input class is 'dist' so 'type' is supposed to be 'dissimilarity'")
  mat <- as.matrix(mat)
  res <- adjClust.matrix(mat, type = "dissimilarity", h = h, 
                         strictCheck = strictCheck, nthreads = nthreads)
  x <- sys.call()
  res$call <- update_call(x, "adjClust")
  return(res)
}

#' @importFrom methods is
#' @import Rcpp
run.adjclust <- function(mat, type = c("similarity", "dissimilarity"), h, 
                         strictCheck = TRUE, nthreads = 1L) {
  # sanity checks
  type <- match.arg(type)
  if (any(is.na(mat)))
    stop("Missing values in the input are not allowed")
  
  p <- nrow(mat)
  
  # 'h'
  if (!is.numeric(h))
    stop("Input band width 'h' must be numeric")
  if (h != as.integer(h))
    stop("Input band width 'h' must be an integer")
  if (h < 0)
    stop("Input band width 'h' must be non negative")
  if (h >= p) 
    stop("Input band width 'h' must be strictly less than dimensions of matrix")
  
  # data preprocessing
  if (type == "dissimilarity") {
    mat <- 1 - 0.5*(mat^2)
  }
  
  if (strictCheck) {
    # this check is *very* expensive for a large number of features
    res_cc <- checkCondition(mat)
    if (is.numeric(res_cc)) {
      message(paste("Note: modifying similarity to ensure positive heights...
        added", res_cc, "to diagonal (merges will not be affected)"))
      mat <- mat + diag(rep(res_cc, ncol(mat)))
    } else res_cc <- 0
  } else res_cc <- 0
  
  if (is(mat, "sparseMatrix")) { 
    # left  
    rCumL <- matL_sparse_rowCumsums(mat, h, nthreads = nthreads)
    rcCumL <- colCumsums(rCumL, useNames = FALSE) # p x (h+1) matrix
    rm(rCumL)

    # right
    rCumR <- matR_sparse_rowCumsums(mat, h, nthreads = nthreads)
    rcCumR <- colCumsums(rCumR, useNames = FALSE) # p x (h+1) matrix
    rm(rCumR)

    out_matL <- matL_sparse(mat, 2)
  } else {
    # left
    rCumL <- matL_full_rowCumsums(mat, h, nthreads = nthreads)
    rcCumL <- colCumsums(rCumL, useNames = FALSE) # p x (h+1) matrix
    rm(rCumL)

    # right
    rCumR <- matR_full_rowCumsums(mat, h, nthreads = nthreads)
    rcCumR <- colCumsums(rCumR, useNames = FALSE) # p x (h+1) matrix
    rm(rCumR)

    out_matL <- matL_full(mat, 2)
  }

  ## Initialization:
  maxSize <- 3*p
  ## NB: The size of heap, D and chainedL are set to the maximum size required,
  ## ie 3*p. This comes from: size p at initialization; at each of the p
  ## iterations, at most 2 new merges are added (only one if the last merge
  ## involves the first or last cluster).

  gains <- rep(0, p-1)
  merge <- matrix(0, nrow = p-1, ncol = 2)  # matrix of the merges
  traceW <- matrix(0, nrow = p-1, ncol = 2) # matrix of traceW
  sd1 <- out_matL[1:(p-1), 2] / 2 # similarity of objects with their right neighbors
  sii <- out_matL[1:p, 1] # auto-similarity of objects
  rm(out_matL)
    
  ## initialization of the heap
  heap <- as.integer(rep(-1, maxSize))
  lHeap <- length(heap)
  v <- 1:(p-1)
  heap[v] <- v
  D <- rep(-1, maxSize)
  
  ## linkage value of objects with their right neighbors
  D[v] <- (sii[v] + sii[v+1]) / 2 - sd1 
  ## initialization of the length of the Heap
  lHeap <- p-1
  ## each element contains a vector: c(cl1, cl2, label1, label2, posL, posR, valid)
  chainedL <- matrix(-1, nrow = 12, ncol = maxSize)
  rownames(chainedL) <- c("minCl1", "maxCl1", "minCl2", "maxCl2", "lab1", 
                          "lab2", "posL", "posR", "Sii", "Sjj", "Sij", "valid")
  w <- as.integer(v + 1)
  chainedL[1,v] <- v
  chainedL[2,v] <- v
  chainedL[3,v] <- w
  chainedL[4,v] <- w
  chainedL[5,v] <- -v
  chainedL[6,v] <- -w
  chainedL[7,v] <- v - 1
  chainedL[8,v] <- w
  chainedL[9,v] <- sii[v]
  chainedL[10,v] <- sii[v+1]
  chainedL[11,v] <- sd1
  chainedL[12,v] <- 1
  chainedL[7,1] <- -1
  chainedL[8,p-1] <- -1

  heap <- buildHeap(heap, D, lHeap)
  # performing clustering  
  res <- .Call("cWardHeaps", rcCumR, rcCumL, as.integer(h), as.integer(p), 
               chainedL, heap, D, as.integer(lHeap), merge, gains, traceW, 
               PACKAGE = "adjclust")
    
  # free memory as soon as it is not needed
  rm(rcCumL)
  rm(rcCumR)

  # formatting outputs (as 'hclust') and checking if decreasing heights are present
  height <- gains
  if (is.null(rownames(mat))) {
    labels <- as.character(1:p)
  } else {
    labels <- rownames(mat)
  }
  tree <- list(merge = res,
               height = height,
               order = 1:p,
               labels = labels,
               call = match.call(),
               method = "adjClust",
               dist.method = attr(D, "method"),
               data = mat, 
               correction = res_cc)
  class(tree) <- c("chac")
    
  if (any(diff(tree$height) < 0)) 
    message(paste("Note:", sum(diff(tree$height) < 0),
                  "merges with non increasing heights."))
    
  return(tree)
}

