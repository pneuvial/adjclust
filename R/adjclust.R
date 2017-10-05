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
#' @param blMin depth of clustering. It is the number of clusters below which 
#'   the algorithm stops. Default value is 1
#' @param verbose Currently not used
#'   
#' @return An object of class \code{\link{chac}}.
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

#' @export
#' @exportClass chac
#' 
#' @importFrom matrixStats rowCumsums
#' @importFrom matrixStats colCumsums

adjClust <- function(mat, type = c("similarity", "dissimilarity"), 
                     h = ncol(mat) - 1, blMin = 1, verbose = FALSE) {
    
    CLASS <- c("matrix", "dgCMatrix", "dsCMatrix", "dist")
    classcheck <- pmatch(class(mat), CLASS)
    
    if (is.na(classcheck)) {
        stop("Input matrix class not supported")
    }
    if (classcheck == -1) {
        stop("Ambiguous matrix class")
    }
    class <- CLASS[classcheck]
    
    type <- match.arg(type)
    
    if (class == "matrix") {
        if (!(nrow(mat) == ncol(mat)))
            stop("Input matrix is not a square matrix")
        if (!is.numeric(mat))
            stop("Input matrix is not numeric")
        if (!(isSymmetric.matrix(mat)))
            stop("Input matrix is not symmetric")
        if (any(is.na(mat)))
            stop("Missing values in the input")
        
        p <- nrow(mat)
        if (h >= p) 
            stop("Input band width should be strictly less than dimensions of matrix")
        
        
        if (type == "dissimilarity") {
            mat <- 1 - 0.5*(mat^2)
        }
        
        ## modifiy (if required) input similarity matrix
        ##  and return a similiarity matrix with diagonal 1
        mat <- modify(mat, as.integer(p), as.integer(h)) 
        
        matL <- findMatL(mat, as.integer(p), as.integer(h))
        rotatedMatR <- findRMatR(mat, as.integer(p), as.integer(h))
        
    } else if (class == "dgCMatrix" || class == "dsCMatrix") {   
        ## dgC/dsC sparse matrices
        
        if (mat@Dim[1] != mat@Dim[2])
            stop("Input matrix is not a square matrix")
        if (any(!(is.numeric(mat@x))))
            stop("Input matrix is not numeric")
        
        p <- mat@Dim[1]
        
        mat <- sparseBand(mat@x, mat@p, mat@i, as.integer(p), as.integer(h))
        mat <- modifySparse(mat, as.integer(p), as.integer(h))
        
        matL <- findSparseMatL(mat, as.integer(p), as.integer(h))
        rotatedMatR <- findSparseRMatR(mat, as.integer(p), as.integer(h))
        
    } else if (class=="dist") { 
        ## for dist objects 
        mat <- as.matrix(mat)
        p <- nrow(mat)
        if (length(h)==0) { 
            ## h defaults to 'ncol(mat)-1' which is 'numeric(0)'
            ## if the *input* 'mat' is of type 'dist'
            h <- ncol(mat)-1
        }
        
        if (h >= p)
            stop("Input band width should be strictly less than dimensions of matrix")
        
        mat <- 1 - 0.5*(mat^2)
        
        matL <- findMatL(mat, p, h)
        rotatedMatR <- findRMatR(mat, p, h)        
    } else {
        stop("Input matrix class not supported")
    }
    
    rCumL <- rowCumsums(matL)          ## p x h matrix
    rcCumL <- colCumsums(rCumL)        ## p x h matrix
    
    rCumR <- rowCumsums(rotatedMatR)   ## p x h matrix
    rcCumR <- colCumsums(rCumR)        ## p x h matrix
    
    ## initialization
    gains <- rep(0, p-blMin)
    merge <- matrix(0, nrow = p-blMin, ncol = 2)   ## matrix of the merges
    traceW <- matrix(0, nrow = p-blMin, ncol = 2)  ## matrix of traceW
    sd1 <- matL[1:(p-1),1]
    
    ## initialization of the heap
    heap <- as.integer(rep(-1, 3*p))
    lHeap <- length(heap)
    v <- 1:(p - 1)
    heap[v] <- v
    D <- rep(-1, 3*p)
    D[v] <- 1 - sd1
    ## initialization of the length of the Heap
    lHeap <- p - 1
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
    chainedL[9,v] <- 1
    chainedL[10,v] <- 1
    chainedL[11,v] <- sd1
    chainedL[12,v] <- 1
    chainedL[7,1] <- -1
    chainedL[8,p-1] <- -1
    heap <- buildHeap(heap, D, lHeap)
    
    res <- .Call("cWardHeaps", rcCumR, rcCumL, as.integer(h), as.integer(p), 
                 chainedL, heap, D, as.integer(lHeap), merge, gains, traceW, 
                 as.integer(blMin), PACKAGE = "adjclust")
    
    height <- cumsum(gains)
    tree <- list(traceW = traceW,
                 gains = gains,
                 merge = res,
                 height = height,
                 seqdist = height,
                 order = 1:p,
                 labels = paste("",1:p),
                 method = "adjClust",
                 call = match.call(),
                 dist.method = attr(D, "method"),
                 data = mat)
    class(tree) <- c("chac")
    return(tree)
}

