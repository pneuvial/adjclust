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
#' @param x either: 1. A pxp contact map of class Matrix::dsCMatrix in which the
#'   entries are the number of counts of physical interactions observed between 
#'   all pairs of loci 2. An object of class HiTC::HTCexp. The corresponding 
#'   Hi-C data is stored as a Matrix::dsCMatrix object in the intdata slot 3. A 
#'   text file with one line per pair of loci for which an interaction has been 
#'   observed (in the format: locus1<tab>locus2<tab>signal).
#'   
#' @param h band width. If not provided, \code{h} is set to default value `p-1`.
#'   
#' @param \dots further arguments to be passed to \code{\link{read.table}} 
#'   function when \code{x} is a text file name. If not provided, the text file 
#'   is supposed to be separated by tabulations, with no header.
#'   
#' @return An object of class \code{\link{chac}}.
#'   
#' @seealso \code{\link{adjClust}} \code{\link[HiTC:HTCexp]{HTCexp}}
#'   
#' @references Dehman A. (2015) \emph{Spatial Clustering of Linkage 
#'   Disequilibrium Blocks for Genome-Wide Association Studies}, PhD thesis, 
#'   Universite Paris Saclay.
#'   
#' @references Servant N. \emph{et al} (2012). \emph{HiTC : Exploration of 
#'   High-Throughput 'C' experiments. Bioinformatics}.
#'   
#' @examples
#' # input as HiTC::HTCexp object
#' data("hic_imr90_40_XX", package="adjclust")
#' 
#' # input as HiTC::HTCexp object
#' res1 <- hicClust(hic_imr90_40_XX)
#' 
#' \dontrun{
#' # input as Matrix::dsCMatrix contact map
#' mat <- HiTC::intdata(hic_imr90_40_XX) 
#' res2 <- hicClust(mat)
#' } 
#' 
#' # input as text file
#' res3 <- hicClust(system.file("extdata", "sample.txt", package = "adjclust"))
#' 
#' @export
#' 
#' @importFrom utils read.table
#' @importFrom HiTC intdata

hicClust <- function(x, h = NULL, ...) {
    
    if (!is.null(h)) {
        if (!is.numeric(h))
            stop("h should be numeric")
    }
    
    class <- class(x)
    if( (class != "dsCMatrix") && (class !=  "HTCexp")&&(!file.exists(x)) )
        stop("Invalid Input:x should be a text file or an object of class Matrix::dsCMatrix/HiTC::HTCexp")
    
    if(class == "dsCMatrix" || class  == "HTCexp") {
        
        if (class == "HTCexp") {
            x <- intdata(x)
        }
        p <- x@Dim[1]
        if(is.null(h)) h <- p-1  
        res <- adjClust(x, type = "similarity", h)
        return(res)
        
    } else {
        
        inoptions <- list(...)
        inoptions$file <- x
        if (is.null(inoptions$sep)) {
            inoptions$sep <- "\t"
        }
        if (is.null(inoptions$header)) {
            inoptions$header <- FALSE
        }
        df <- do.call("read.table", inoptions) 
        
        lis <- sort(unique(c(df[,1], df[,2])))
        p <- length(lis)
        rowindx <- match(df[,1], lis)
        colindx <- match(df[,2], lis)
        
        m <- matrix(0, nrow = p, ncol = p)
        m[cbind(rowindx,colindx)] <- m[cbind(colindx,rowindx)] <- df[,3]
        
        if (is.null(h)) h <- p-1  
        res <- adjClust(m, type = "similarity", h = h)
        res$method <- "hicClust"
        
        return(res)
    }
}
