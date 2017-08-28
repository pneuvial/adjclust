#' Constrained Hierarchical Agglomerative Clustering
#' 
#' Function to perform adjacency-constrained hierarchical agglomerative clustering
#' 
#' \code{adjClust} performes constrained hierarchichal agglomerative clustering which is   
#' hierarchical agglomerative clustering in which each observation is associated to a position, 
#' and the clustering is constrained so as only adjacent clusters are merged. These methods are  
#' useful in various application fields, including ecology (Quaternary data) and bioinformatics 
#' (e.g. in Genome-Wide Association Studies (GWAS))
#' 
#' @param mat A similarity matrix or a dist object
#' @param type Type of matrix : similarity or dissimilarity
#' @param h band width. It is assumed that the similarity between two items is 
#' 0 when these items are at a distance of more than band width h
#' @param blMin depth of clustering. It is number of clusters at which the algorithm stops. Default value is 1.
#' @param verbose Currently not used
#' 
#' @return Function \code{adjClust} returns an object of
#' class \code{\link[stats]{hclust}}.  
#'
#' @examples
#' sim <- matrix(c(1,0.1,0.2,0.3,0.1,1,0.4,0.5,0.2,0.4,1,0.6,0.3,0.5,0.6,1), nrow=4)
#' h <- 3
#' fit1 <- adjClust(sim, "similarity", h, 1, FALSE)
#' plot(fit1)
#' 
#' dist <- as.dist(sqrt(2-(2*sim)))
#' 
#' #Compatibility with dist objects
#' fit2 <- adjClust(dist, "dissimilarity", h, 1, FALSE)
#' plot(fit2)
#' 
#' @export
#' 
#' @importFrom matrixStats rowCumsums
#' @importFrom matrixStats colCumsums

adjClust <- function(mat, type = "similarity", h, blMin=1, verbose=FALSE){
    
    if (!is.numeric(h))
      stop("Input band width is not numeric")
    
  
    CLASS <- c("matrix","dgCMatrix","dsCMatrix","dist")
    classcheck <- pmatch(class(mat), CLASS)
    
    if(is.na(classcheck))
      stop("Input matrix class not supported")
    if(classcheck == -1)
      stop("Ambiguous matrix class")
    class <- CLASS[classcheck]
    
    
    TYPES <- c("similarity", "dissimilarity")
    typecheck <- pmatch(type, TYPES)
    
    if(is.na(typecheck))
      stop("Invalid matrix type")
    if(typecheck == -1)
      stop("Ambiguous matrix type")
    type <- TYPES[typecheck]
    
    
    if (class == "matrix")  ## for standard matrices
    {
      if (!(nrow(mat) == ncol(mat)))
        stop("Input matrix is not a square matrix")
      if (!is.numeric(mat))
        stop("Input matrix is not numeric")
      if (!(isSymmetric.matrix(mat)))
        stop("Input matrix is not symmetric")
      if (any(is.na(mat)))
        stop("Missing values in the input")
      
      p <- nrow(mat)
      if (h >= p) {
        stop("Input band width should be strictly less than dimensions of matrix")
      }
      
      
      if (type == "dissimilarity") {
        mat <- 1 - 0.5*(mat^2)
      }
      
      mat <- modify(mat, as.integer(p), as.integer(h)) ## modifies, if required, input similarity matrix and returns a similiarity matrix with diagonal 1
      
      matL <- findMatL(mat, as.integer(p), as.integer(h))
      rotatedMatR <- findRMatR(mat, as.integer(p), as.integer(h))
        
    } else if(class == "dgCMatrix" || class == "dsCMatrix") {   ## for dgC/dsC sparse matrices
      
      if (mat@Dim[1] != mat@Dim[2])
        stop("Input matrix is not a square matrix")
      if (any(!(is.numeric(mat@x))))
        stop("Input matrix is not numeric")
      
      p <- mat@Dim[1]
      
      mat <- sparseBand(mat@x, mat@p, mat@i, as.integer(p), as.integer(h))
      mat <- modifySparse(mat, as.integer(p), as.integer(h))
      
      matL <- findSparseMatL(mat, as.integer(p), as.integer(h))
      rotatedMatR <- findSparseRMatR(mat, as.integer(p), as.integer(h))
        
    } else { ## for dist objects 
      
      mat <- as.matrix(mat)
      p <- nrow(mat)
      if (h >= p)
        stop("Input band width should be strictly less than dimensions of matrix")
      
      mat <- 1 - 0.5*(mat^2)
      
      matL <- findMatL(mat, p, h)
      rotatedMatR <- findRMatR(mat, p, h)        
    }
    
    rCumL <- rowCumsums(matL)         ## p x h matrix
    rcCumL <- colCumsums(rCumL)        ## p x h matrix

    rCumR <- rowCumsums(rotatedMatR)  ## p x h matrix
    rcCumR <- colCumsums(rCumR)  ## p x h matrix

    ## initialization
    gains <- rep(0, p-blMin)
    merge <- matrix(0, nrow=p-blMin, ncol=2)  ## matrix of the merges
    traceW <- matrix(0, nrow=p-blMin, ncol=2)  ## matrix of traceW
    sd1 <- matL[1:(p-1),1]

    ## initialization of the heap
    heap <- as.integer(rep(-1, 3*p))
    lHeap <- length(heap)
    heap[1:(p-1)] <- 1:(p-1)
    D <- rep(-1, 3*p)
    D[1:(p-1)] <- 1-sd1
    ## initialization of the length of the Heap
    lHeap <- p-1
    ## each element contains a vector: c(cl1, cl2, label1, label2, posL, posR, valid)
    chainedL <- matrix(-1, nrow=12, ncol=3*p)
    rownames(chainedL) <- c("minCl1", "maxCl1", "minCl2", "maxCl2", "lab1", "lab2", "posL", "posR", "Sii", "Sjj", "Sij", "valid")
    v <- 1:(p-1)
    w <- as.integer(v+1)
    chainedL[1,v] <- v
    chainedL[2,v] <- v
    chainedL[3,v] <- w
    chainedL[4,v] <- w
    chainedL[5,v] <- -v
    chainedL[6,v] <- -w
    chainedL[7,v] <- v-1
    chainedL[8,v] <- w
    chainedL[9,v] <- 1
    chainedL[10,v] <- 1
    chainedL[11,v] <- sd1
    chainedL[12,v] <- 1
    chainedL[7,1] <- -1
    chainedL[8,p-1] <- -1
    heap <- buildHeap(heap, D, lHeap)

    res <- .Call("cWardHeaps", rcCumR, rcCumL, as.integer(h), as.integer(p), chainedL, heap, D, as.integer(lHeap), merge, gains, traceW, as.integer(blMin), PACKAGE="adjclust")

    height <- cumsum(gains)
    tree <- list(traceW=traceW,
                 gains=gains,
                 merge = res,
                 height = height,
                 seqdist = height,
                 order = 1:p,
                 labels = paste("",1:p),
                 method = "adjclust-heaps",
                 call = match.call(),
                 dist.method = attr(D, "method"))
    class(tree) <- "hclust"
    return(tree)
}
