#' Constrained Hierarchical Agglomerative Clustering of Single Nucleotide Polymorphisms
#' 
#' Function to perform adjacency-constrained hierarchical agglomerative clustering of Single Nucleotide Polymorphisms
#' 
#' \code{snpClust} performs constrained hierarchichal agglomerative clustering of single nucleotide 
#' polymorphisms. Constrained Hierarchical Agglomerative Clustering is hierarchical agglomerative
#' clustering in which each observation is associated to a position, and the clustering is 
#' constrained so as to merge only adjacent clusters.
#' 
#' @param x either a genotype matrix of class snpStats::SnpMatrix/base::matrix or a linkage disequlibrium matrix of 
#' class Matrix::dgCMatrix
#' @param h band width. It is assumed that the similarity between two items is 0 when these items are at a 
#' distance of more than band width h
#' @param stats a character vector specifying the linkage disequilibrium measure to be calculated. This should contain one or 
#' more of the strings: "LLR", "OR", "Q", "Covar", "D.prime", "R.squared" or "R"
#' 
#' @return Function \code{snpClust} returns an object of class \code{\link[stats]{hclust}}.  
#' 
#' @examples
#' library(snpStats)
#' data(testdata)
#' 
#' #Input as snpStats::SnpMatrix
#' fit1 <- snpClust(Autosomes[1:200, 1:5], 3, "R.squared")
#' 
#' #Input as base::matrix
#' fit2 <- snpClust(as.matrix(Autosomes[1:200, 1:5]), 3, "R.squared")
#' 
#' #Input as Matrix::dgCMatrix
#' ld <- ld(Autosomes[1:200, 1:5], depth=3, stats="R.squared")
#' fit3 <- snpClust(ld, 3)
#' 
#' @export
#' 
#' @importFrom methods as
#' @importFrom snpStats ld
snpClust <- function(x ,h ,stats) {
  
  if (!is.numeric(h))
    stop("h should be numeric")
  
  CLASS <- c("dgCMatrix", "matrix", "SnpMatrix")
  classcheck <- pmatch(class(x), CLASS)
  
  if(is.na(classcheck))
    stop("Input matrix class not supported")
  if(classcheck == -1)
    stop("Ambiguous matrix class")
  class <- CLASS[classcheck]  

  if(class != "dgCMatrix" )
  {
    p <- ncol(x)
    if (h >= p)
      stop("h should be strictly less than p")
    
    if (class == "matrix") {
      rownames(x) <- 1:nrow(x)
      colnames(x) <- 1:ncol(x)
      x <- as(x, "SnpMatrix")
    }
    
    if (missing(stats)) 
      stop("LD stats must be specified")
    
    STATS <- c("LLR", "OR", "Q", "Covar", "D.prime", "R.squared", "R")
    statscheck <- pmatch(stats, STATS)
    
    if(is.na(statscheck))
      stop("Invalid stats")
    if(statscheck == -1)
      stop("Ambiguous stats")
    stats <- STATS[statscheck] 
    
    x <- ld(x, stats = stats, depth = h)
    x <- round(x, digits=10)
    
    if (any(is.na(x)))
      stop("Missing value(s) or NaN(s) generated in ld calculation. Please check the input for missing observations.")    
  }
  
  diag(x) <- 1
  res <- adjClust(x, "similarity", h, 1, FALSE)
  
  return(res)
}
