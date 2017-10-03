#' Constrained Hierarchical Agglomerative Clustering of Single Nucleotide 
#' Polymorphisms
#' 
#' Function to perform adjacency-constrained hierarchical agglomerative 
#' clustering of Single Nucleotide Polymorphisms
#' 
#' \code{snpClust} performs constrained hierarchichal agglomerative clustering 
#' of single nucleotide polymorphisms. Constrained Hierarchical Agglomerative 
#' Clustering is hierarchical agglomerative clustering in which each observation
#' is associated to a position, and the clustering is constrained so as to merge
#' only adjacent clusters.
#' 
#' @param x either a genotype matrix of class snpStats::SnpMatrix/base::matrix 
#'   or a linkage disequlibrium matrix of class Matrix::dgCMatrix
#'   
#' @param h band width. If not provided, \code{h} is set to default value `p-1`
#'   where `p` is the number of columns of `x`
#'   
#' @param \dots Further arguments to be passed to the \code{snpStats::ld}
#'   function
#'   
#' @return Function \code{snpClust} returns an object of class 
#'   \code{\link[stats]{hclust}}.
#'   
#' @seealso \code{\link{adjClust}} \code{\link[snpStats:ld]{ld}}
#'   
#' @examples
#' ## a very small example
#' data(testdata, package = "snpStats")
#' 
#' # input as snpStats::SnpMatrix
#' fit1 <- snpClust(Autosomes[1:200, 1:5], h = 3, stats = "R.squared")
#' 
#' # input as base::matrix
#' fit2 <- snpClust(as.matrix(Autosomes[1:200, 1:5]), h = 3, stats = "R.squared")
#' 
#' # input as Matrix::dgCMatrix
#' ld <- snpStats::ld(Autosomes[1:200, 1:5], depth = 3, stats = "R.squared")
#' fit3 <- snpClust(ld, 3)
#' 
#' @export
#' 
#' @importFrom methods as
#' @importFrom snpStats ld
#'   
snpClust <- function(x, h = ncol(x) - 1, ...) {
    
    if (!is.numeric(h)) {
        stop("h should be numeric")
    }
    
    CLASS <- c("dgCMatrix", "matrix", "SnpMatrix")
    classcheck <- pmatch(class(x), CLASS)
    
    if (is.na(classcheck)) {
        stop("Input matrix class not supported")
    }
    if (classcheck == -1) {
        stop("Ambiguous matrix class")
    }
    class <- CLASS[classcheck]  
    
    if (class != "dgCMatrix" ) {
        p <- ncol(x)
        if (h >= p) {
            stop("h should be strictly less than p")
        }
        if (class == "matrix") {
            rownames(x) <- 1:nrow(x)
            colnames(x) <- 1:ncol(x)
            x <- as(x, "SnpMatrix")
        }
        
        x <- ld(x, ..., depth = h)
        x <- round(x, digits = 10)
        if (any(is.na(x))) {
            ww <- which(is.na(as.matrix(x)), arr.ind = TRUE)
            print(utils::str(ww))
            stop("Missing value(s) or NaN(s) in LD estimates (indices are printed above)")
        }
    }
    if (any(is.na(x))) {
        stop("Missing value(s) or NaN(s) not allowed in similarity matrix.")    
    }
    diag(x) <- 1
    res <- adjClust(x, type = "similarity", h = h)
    
    return(res)
}
