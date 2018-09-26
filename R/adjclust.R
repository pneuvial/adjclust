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
#' implementation are available in the package vignette entitled "Notes on CHAC 
#' implementation in adjclust".
#' 
#' @param mat A similarity matrix or a dist object
#' @param type Type of matrix : similarity or dissimilarity. Defaults to 
#'   \code{"similarity"}
#' @param h band width. It is assumed that the similarity between two items is 0
#'   when these items are at a distance of more than band width h. Default value
#'   is \code{ncol(mat)-1}
#'   
#' @return An object of class \code{\link{chac}} which describes the tree 
#' produced by the clustering process. The object a list with the same elements 
#' as an object of class \code{\link{chac}} (\code{merge}, \code{height}, 
#' \code{order}, \code{labels}, \code{call}, \code{method}, \code{dist.method}),
#' and an extra element \code{mat}: the data on which the clustering is 
#' performed, possibly after pre-transformations described in the vignette 
#' entitled "Notes on CHAC implementation in adjclust".
#'   
#' @seealso \code{\link{snpClust}} to cluster SNPs based on linkage disequilibrium
#' @seealso \code{\link{hicClust}} to cluster Hi-C data
#'   
#' @references Dehman A. (2015) \emph{Spatial Clustering of Linkage 
#'   Disequilibrium Blocks for Genome-Wide Association Studies}, PhD thesis, 
#'   Universite Paris Saclay.
#'   
#' @examples
#' sim <- matrix(
#' c(1.0, 0.1, 0.2, 0.3,
#'   0.1, 1.0 ,0.4 ,0.5,
#'   0.2, 0.4, 1.0, 0.6, 
#'   0.3, 0.5, 0.6, 1.0), nrow = 4)
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
#' @exportClass chac
#' 
#' @importFrom matrixStats rowCumsums
#' @importFrom matrixStats colCumsums
#' @importFrom Matrix diag
#' @importFrom Matrix t
#' @importFrom Matrix isSymmetric

adjClust <- function(mat, type = c("similarity", "dissimilarity"), 
                     h = ncol(mat) - 1) {
  UseMethod("adjClust")
}

#' @export
adjClust.matrix <- function(mat, type = c("similarity", "dissimilarity"), 
                            h = ncol(mat) - 1) {
  if (!is.numeric(mat))
    stop("Input matrix is not numeric")
  if (!(isSymmetric(mat)))
    stop("Input matrix is not symmetric")
  res <- run.adjclust(mat, type = type, h = h)
  return(res)
}

#' @export
adjClust.Matrix <- function(mat, type = c("similarity", "dissimilarity"), 
                            h = ncol(mat) - 1) {
  if (!(isSymmetric(mat)))
    stop("Input matrix is not symmetric")
  res <- run.adjclust(mat, type = type, h = h)
  return(res)
}

#' @export
adjClust.dgCMatrix <- function(mat, type = c("similarity", "dissimilarity"), 
                            h = ncol(mat) - 1) {
  type <- match.arg(type)
  if (type == "dissimilarity")
    stop("'type' can only be 'similarity' with sparse Matrix inputs")
  res <- run.adjclust(mat, type = type, h = h)
  return(res)
}

#' @export
adjClust.dsCMatrix <- function(mat, type = c("similarity", "dissimilarity"), 
                               h = ncol(mat) - 1) {
  type <- match.arg(type)
  if (!(isSymmetric(mat)))
    stop("Input matrix is not symmetric")
  if (type == "dissimilarity")
    stop("'type' can only be 'similarity' with sparse Matrix inputs")
  res <- run.adjclust(mat, type = type, h = h)
  return(res)
}

#' @export
adjClust.dgeMatrix <- function(mat, type = c("similarity", "dissimilarity"), 
                               h = ncol(mat) - 1) {
  type <- match.arg(type)
  if (!(isSymmetric(mat)))
    stop("Input matrix is not symmetric")
  if (type == "dissimilarity")
    stop("'type' can only be 'similarity' with sparse Matrix inputs")
  res <- run.adjclust(mat, type = type, h = h)
  return(res)
}

#' @export
adjClust.dgTMatrix <- function(mat, type = c("similarity", "dissimilarity"), 
                               h = ncol(mat) - 1) {
  type <- match.arg(type)
  if (!(isSymmetric(mat)))
    stop("Input matrix is not symmetric")
  if (type == "dissimilarity")
    stop("'type' can only be 'similarity' with sparse Matrix inputs")
  res <- run.adjclust(mat, type = type, h = h)
  return(res)
}

#' @export
adjClust.dist <- function(mat, type = c("similarity", "dissimilarity"), 
                          h = ncol(mat) - 1) {
  type <- match.arg(type)
  if (type != "dissimilarity")
    message("Note: input class is 'dist' so 'type' is supposed to be 'dissimilarity'")
  mat <- as.matrix(mat)
  res <- adjClust.matrix(mat, type = "dissimilarity", h = h)
  return(res)
}

run.adjclust <- function(mat, type = c("similarity", "dissimilarity"), h) {
  # sanity checks
  type <- match.arg(type)
  if (!(nrow(mat) == ncol(mat)))
    stop("Input matrix is not a square matrix")
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
  
  res_cc <- checkCondition(mat)
  if (is.numeric(res_cc)) {
    message(paste("Note: modifying similarity to ensure positive heights...
      added", res_cc, "to diagonal (merges will not be affected)"))
    mat <- mat + diag(rep(res_cc, ncol(mat)))
  }
  
  out_matL <- matL(mat, h)
  out_matR <- matR(mat, h)
  
  ## computing pencils
  rCumL <- rowCumsums(out_matL) # p x (h+1) matrix
  rcCumL <- colCumsums(rCumL) # p x (h+1) matrix
    
  rCumR <- rowCumsums(out_matR) # p x (h+1) matrix
  rcCumR <- colCumsums(rCumR) # p x (h+1) matrix
  
  ## initialization (heap, D and chainedL are too large in order to avoid memory pb in C)
  gains <- rep(0, p-1)
  merge <- matrix(0, nrow = p-1, ncol = 2) # matrix of the merges
  traceW <- matrix(0, nrow = p-1, ncol = 2) # matrix of traceW
  sd1 <- out_matL[1:(p-1),2]/2 # similarity of objects with their right neighbors
  sii <- out_matL[1:p,1] # auto-similarity of objects
    
  ## initialization of the heap
  heap <- as.integer(rep(-1, 3*p))
  lHeap <- length(heap)
  v <- 1:(p-1)
  heap[v] <- v
  D <- rep(-1, 3*p)
  
  ## linkage value of objects with their right neighbors
  D[v] <- (sii[1:(p-1)] + sii[2:p]) / 2 - sd1 
  ## initialization of the length of the Heap
  lHeap <- p-1
  ## each element contains a vector: c(cl1, cl2, label1, label2, posL, posR, valid)
  chainedL <- matrix(-1, nrow = 12, ncol = 3*p)
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
  chainedL[9,v] <- sii[1:(p-1)]
  chainedL[10,v] <- sii[2:p]
  chainedL[11,v] <- sd1
  chainedL[12,v] <- 1
  chainedL[7,1] <- -1
  chainedL[8,p-1] <- -1

  heap <- buildHeap(heap, D, lHeap)
  
  # performing clustering  
  res <- .Call("cWardHeaps", rcCumR, rcCumL, as.integer(h), as.integer(p), 
               chainedL, heap, D, as.integer(lHeap), merge, gains, traceW, 
               PACKAGE = "adjclust")
    
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
               data = mat)
  class(tree) <- c("chac")
    
  if (any(diff(tree$height) < 0)) 
    message(paste("Note:", sum(diff(tree$height) < 0),
                  "merges with non increasing heights."))
    
  return(tree)
}

