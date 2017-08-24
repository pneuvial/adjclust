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
#' @param h band width. It is assumed that the similarity between two items is 0 when these 
#' items are at a distance of more than band width h
#' 
#' @param \dots further arguments to be passed to read.table function.It can contain any argument(s) of 
#' \code{\link{read.table}} except the file argument.
#'  
#' @return Function \code{hicclust} returns an object of class \code{\link[stats]{hclust}}.  
#' 
#' @examples
#' #Input as HiTC::HTCexp object
#' library("HiTC")
#' data("hic_imr90_40$chrXchrX", package="adjclust")
#' h <- 3881
#' res1 <- hicclust(obj, h)
#' 
#' #Input as Matrix::dsCMatrix contact map
#' mat <- intdata(obj) 
#' res2 <- hicclust(mat, h)
#' 
#' #Input as text file
#' h <- 5
#' res3 <- hicclust(system.file("extdata", "sample.txt", package = "adjclust"), h)
#' 
#' @export 
#' 
#' @importFrom utils read.table
#' @importFrom HiTC intdata

hicclust <- function(x, h = NULL, ...) {

  if ((h!=NULL)&&(!is.numeric(h)))
    stop("h should be numeric")
  
  class <- class(x)
  if( (class != "dsCMatrix")&&(class !=  "HTCexp")&&(!file.exists(x)) )
    stop("Invalid Input:x should be a text file or an object of class Matrix::dsCMatrix/HiTC::HTCexp")
  
  if(class == "dsCMatrix" || class  == "HTCexp") {
  
    if(class == "HTCexp")
      x <- intdata(x)
    
    p <- mat@Dim[1]
    if(h == NULL) h <- p-1  
    res <- adjClustBand_heap(x, type = "similarity", h)
    return(res)
  
  } else {
  
  inoptions <- list(...)
  inoptions$file <- x
  if (is.null(inoptions$sep)) inoptions$sep <- "\t"
  if (is.null(inoptions$header)) inoptions$header <- FALSE  
#  rtargs <- c(x, inoptions)
#  df <- do.call("read.table", rtargs)
  df <- do.call("read.table", inoptions) 
    
  lis <- sort(unique(c(df[,1], df[,2])))
  p <- length(lis)
  rowindx <- match(df[,1], lis)
  colindx <- match(df[,2], lis)
  
  m <- matrix(0, nrow = p, ncol = p)
  m[cbind(rowindx,colindx)] <- m[cbind(colindx,rowindx)] <- df[,3]

  if(h == NULL) h <- p-1  
  res <- adjClustBand_heap(m, type = "similarity", h)
  return(res)
  
  }
}
