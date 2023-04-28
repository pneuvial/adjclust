#' Adjacency-constrained Clustering of Single Nucleotide Polymorphisms
#'
#' Adjacency-constrained hierarchical agglomerative clustering of Single
#' Nucleotide Polymorphisms based on Linkage Disequilibrium
#'
#' Adjacency-constrained hierarchical agglomerative clustering (HAC) is HAC in
#' which each observation is associated to a position, and the clustering is
#' constrained so as only adjacent clusters are merged. SNPs are clustered based
#' on their similarity as measured by the linkage disequilibrium.
#'
#' In the special case where genotypes are given as input and the corresponding
#' LD matrix has missing entries, the clustering cannot be performed. This can
#' typically happen when there is insufficient variability in the sample
#' genotypes. In this special case, the indices of the SNP pairs which yield
#' missing values are returned.
#'
#' @param x either a genotype matrix of class
#'   \code{\link[snpStats:SnpMatrix-class]{SnpMatrix}}/\code{\link{matrix}} or a
#'   linkage disequilibrium matrix of class
#'   \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}. In the latter case the LD
#'   values are expected to be in [0,1]
#'
#' @param h band width. If not provided, \code{h} is set to default value `p-1`
#'   where `p` is the number of columns of \code{x}
#'
#' @param stats a character vector specifying the linkage disequilibrium
#'   measures to be calculated (using the \code{\link[snpStats:ld]{ld}}
#'   function) when \code{x} is a genotype matrix. Only "R.squared" and
#'   "D.prime" are allowed, see Details.
#'
#' @return An object of class \code{\link{chac}} (when no LD value is missing)
#'
#' @seealso \code{\link{adjClust}} \code{\link[snpStats:ld]{ld}}
#'
#' @references Dehman A. (2015) \emph{Spatial Clustering of Linkage
#'   Disequilibrium Blocks for Genome-Wide Association Studies}, PhD thesis,
#'   Universite Paris Saclay.
#'
#' @references Dehman, A. Ambroise, C. and Neuvial, P. (2015). Performance of a
#'   blockwise approach in variable selection using linkage disequilibrium
#'   information. *BMC Bioinformatics* 16:148.
#'   
#' @references Ambroise C., Dehman A., Neuvial P., Rigaill G., and Vialaneix N
#'   (2019). \emph{Adjacency-constrained hierarchical clustering of a band
#'   similarity matrix with application to genomics}, Algorithms for Molecular
#'   Biology 14(22)"
#'
#'
#' @details If \code{x} is of class
#'   \code{\link[snpStats:SnpMatrix-class]{SnpMatrix}} or \code{\link{matrix}},
#'   it is assumed to be a \eqn{n \times p} matrix of \eqn{p} genotypes for
#'   \eqn{n} individuals. This input is converted to a LD similarity matrix
#'   using the \code{snpStats::ld}. If \code{x} is of class
#'   \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}, it is assumed to be a
#'   (squared) LD matrix.
#'
#'   Clustering on a LD similarity other than "R.squared" or "D.prime" can be
#'   performed by providing the LD values directly as argument \code{x}. These
#'   values are expected to be in [0,1], otherwise they are truncated to [0,1].
#'
#'
#' @examples
#' ## a very small example
#' if (requireNamespace("snpStats", quietly = TRUE)) {
#'   data(testdata, package = "snpStats")
#'
#'   # input as snpStats::SnpMatrix
#'   fit1 <- snpClust(Autosomes[1:200, 1:5], h = 3, stats = "R.squared")
#'
#'   # input as base::matrix
#'   fit2 <- snpClust(as.matrix(Autosomes[1:200, 1:5]), h = 3, stats = "R.squared")
#'
#'   # input as Matrix::dgCMatrix
#'   ldres <- snpStats::ld(Autosomes[1:200, 1:5], depth = 3, stats = "R.squared", symmetric = TRUE)
#'   fit3 <- snpClust(ldres, 3)
#' }
#'
#' @export
#'
#' @importFrom methods as
#'   

snpClust <- function(x, h = ncol(x) - 1, stats = c("R.squared", "D.prime")) {
  UseMethod("snpClust")
}

#' @export
snpClust.matrix <- function(x, h = ncol(x) - 1, 
                            stats = c("R.squared", "D.prime")) {
  if (!requireNamespace("snpStats"))
    stop("Package 'snpStats' not available. This function cannot be used with 'matrix' data.")
  if (is.null(rownames(x)))
    rownames(x) <- 1:nrow(x)
  if (is.null(colnames(x)))
    colnames(x) <- 1:ncol(x)
  x <- as(x, "SnpMatrix")
  res <- snpClust.snpStats(x, h = h, stats = stats)
  x <- sys.call()
  res$call <- update_call(x, "snpClust")
  return(res)
}

#' @export
snpClust.dgCMatrix <- function(x, h = ncol(x) - 1, 
                               stats = c("R.squared", "D.prime")) {
  res <- run.snpClust(x, h = h, stats = stats)
  x <- sys.call()
  res$call <- update_call(x, "snpClust")
  return(res)
}

#' @export
snpClust.dsCMatrix <- function(x, h = ncol(x) - 1, 
                               stats = c("R.squared", "D.prime")) {
  res <- run.snpClust(x, h = h, stats = stats)
  x <- sys.call()
  res$call <- update_call(x, "snpClust")
  return(res)
}


#' @export
snpClust.snpStats <- function(x, h = ncol(x) - 1, 
                              stats = c("R.squared", "D.prime")) {
  if (!requireNamespace("snpStats"))
    stop("Package 'snpStats' not available. This function cannot be used with 'matrix' data.")
  p <- ncol(x)
  if (h >= p) {
    stop("h should be strictly less than p")
  }
  stats <- match.arg(stats)
  x <- snpStats::ld(x, stats = stats, depth = h, symmetric = TRUE)
  diag(x) <- rep(1, p)
  if (any(is.na(x))) {
    ww <- which(is.na(as.matrix(x)), arr.ind = TRUE)
    warning("Clustering could not be performed due to missing value(s) or NaN(s) in LD estimates. Returning these indices.")
    return(ww)
  }
  res <- run.snpClust(x, h = h, stats = stats)
  x <- sys.call()
  res$call <- update_call(x, "snpClust")
  return(res)
}

run.snpClust <- function(x, h, stats) {
  if (any(is.na(x))) {
    stop("Missing value(s) or NaN(s) not allowed in LDs (and found some).")    
  }
  if (!all(diag(x) == 1)) {
    message("Note: forcing the diagonal of the LD similarity matrix to be 1")
    diag(x) <- rep(1, length(diag(x)))
  }
  out <- (x > 1)
  if (any(out)) {
    warning("Forcing the LD similarity to be smaller than or equal to 1")
    x[out] <- 1
  }
  out <- x < 0
  if (any(out)) {
    warning("Forcing the LD similarity to be larger than or equal to 0")
    x[out] <- 0
  }

  res <- adjClust(x, type = "similarity", h = h)
  res$method <- "snpClust"
  x <- sys.call()
  res$call <- update_call(x, "snpClust")
  return(res)
}