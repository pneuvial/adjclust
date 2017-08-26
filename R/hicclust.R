#' Constrained Hierarchical Agglomerative Clustering of Genomic regions(loci)
#' 
#' Function to perform adjacency-constrained hierarchical agglomerative clustering 
#' of genomic regions(loci)
#' 
#' \code{hicclust} performs constrained hierarchichal agglomerative clustering of 
#' genomic regions (loci) according to information provided by high-throughput conformation capture
#'  data (Hi-C). Constrained Hierarchical Agglomerative Clustering is hierarchical agglomerative 
#'  clustering in which each observation is associated to a position, and the clustering is 
#'  constrained so as only adjacent clusters are merged.
#' 
#' @param x x can be:
#' 1. A pxp contact map of class Matrix::dsCMatrix in which the entries are the number of counts of physical interactions observed between all pairs of loci
#' 2. An object of class HiTC::HTCexp. The corresponding Hi-C data is stored as a Matrix::dsCMatrix object in the intdata slot
#' 3. A text file with one line per pair of loci for which an interaction has been observed (in the format: locus1<tab>locus2<tab>signal). 
#' 
#' @param h band width. If NULL `h` is set to default value `p-1`.It is assumed that the similarity between two items is 0 when these 
#' items are at a distance of more than band width h
#' 
#' @param \dots further arguments to be passed to \code{\link{read.table}} function. If NULL the text file is supposed to be separated by tab with no header.
#'  
#' @return Function \code{hicclust} returns an object of class \code{\link[stats]{hclust}}.  
#' 
#' @examples
#' #Input as HiTC::HTCexp object
#' library("HiTC")
#' data("hic_imr90_40$chrXchrX", package="adjclust")
#' 
#' #Removing rows and columns containing only zeros
#' selected <- apply(intdata(obj), 1, sum) > 0
#' intd <- intdata(obj)[selected,selected]
#' x_int <- x_intervals(obj)[selected,]
#' y_int <- y_intervals(obj)[selected,]
#' obj <- new("HTCexp", intd, x_int, y_int)
#' 
#' \dontrun{
#' #Input as HiTC::HTCexp object
#' res1 <- hicclust(obj)
#' 
#' #Input as Matrix::dsCMatrix contact map
#' mat <- intdata(obj) 
#' res2 <- hicclust(mat)
#' } 
#'
#' #Input as text file
#' res3 <- hicclust(system.file("extdata", "sample.txt", package = "adjclust"))
#' 
#' @export 
#' 
#' @importFrom utils read.table
#' @importFrom HiTC intdata

hicclust <- function(x, h = NULL, ...) {

  if (!is.null(h)) {
    if (!is.numeric(h))
    stop("h should be numeric")
  }
  
  class <- class(x)
  if( (class != "dsCMatrix")&&(class !=  "HTCexp")&&(!file.exists(x)) )
    stop("Invalid Input:x should be a text file or an object of class Matrix::dsCMatrix/HiTC::HTCexp")
  
  if(class == "dsCMatrix" || class  == "HTCexp") {
  
    if(class == "HTCexp")
      x <- intdata(x)
    
    p <- x@Dim[1]
    if(is.null(h)) h <- p-1  
    res <- adjClustBand_heap(x, type = "similarity", h)
    return(res)
  
  } else {
  
  inoptions <- list(...)
  inoptions$file <- x
  if (is.null(inoptions$sep)) inoptions$sep <- "\t"
  if (is.null(inoptions$header)) inoptions$header <- FALSE  
  df <- do.call("read.table", inoptions) 
    
  lis <- sort(unique(c(df[,1], df[,2])))
  p <- length(lis)
  rowindx <- match(df[,1], lis)
  colindx <- match(df[,2], lis)
  
  m <- matrix(0, nrow = p, ncol = p)
  m[cbind(rowindx,colindx)] <- m[cbind(colindx,rowindx)] <- df[,3]

  if(is.null(h)) h <- p-1  
  res <- adjClustBand_heap(m, type = "similarity", h)
  return(res)
  
  }
}
