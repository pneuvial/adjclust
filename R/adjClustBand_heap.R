#' Constrained Hierarchical Agglomerative Clustering
#' 
#' Function to perform adjacency-constrained hierarchical agglomerative clustering
#' 
#' \code{adjClustBand_heap} performes constrained hierarchichal agglomerative clustering which is   
#' hierarchical agglomerative clustering in which each observation is associated to a position, 
#' and the clustering is constrained so as only adjacent clusters are merged. These methods are  
#' useful in various application fields, including ecology (Quaternary data) and bioinformatics 
#' (e.g. in Genome-Wide Association Studies (GWAS))
#' 
#' @param mat A similarity matrix or a dist object
#' @param h band width. It is assumed that the similarity between two items is 
#' 0 when these items are at a distance of more than band width h
#' @param blMin to be added
#' @param verbose to be added
#' 
#' @return Function \code{adjClustBand_heap} returns an object of
#' class \code{\link[stats]{hclust}}.  
#'
#' @examples
#' 
#' h <- 100
#' fit1 <- adjClustBand_heap(mat, h,1,FALSE)
#' plot(fit1, \dots)
#' 
#' #Compatibility with dist objects
#' fit2 <- adjClustBand_heap(dist, h,1,FALSE)
#' plot(fit2, \dots)
#' 
#' @export
#' 
#' @importFrom matrixStats rowCumsums
#' @importFrom matrixStats colCumsums
adjClustBand_heap <- function(mat, h, blMin=1, verbose=FALSE){

    if (!is.numeric(h))
      stop("Input band width is not numeric")
    
    mat <- modify(mat) # modifies input, if required, to convert into similarity matrix with diagonal 1
    
    p <- nrow(mat)
    
    if (h >= p)
      stop("Input band width should be strictly less than dimensions of matrix")    
    
    x <- band(mat, h)
    
    len <- length(x)
    stopifnot(len==(p-1)*h-h*(h-1)/2)
    xt <- transpose(x, p, h)

    ## sum of the "rectangles" beginning from the left
    matL <- .toMatLeft(xt, p, h)     ## a matrix p x h (with zeros at the bottom) of the LD values
    rCumL <- rowCumsums(matL)         ## p x h matrix
    rcCumL <- colCumsums(rCumL)        ## p x h matrix

    ## sum of the "rectangles" beginning from the right
    matR <- .toMatRight(x, p, h)
    rotatedMatR <- .rotate(.rotate(matR))  ## WTF ??
    rCumR <- rowCumsums(rotatedMatR)  ## p x h matrix
    rcCumR <- colCumsums(rCumR)  ## p x h matrix

    rm(x, xt)

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
