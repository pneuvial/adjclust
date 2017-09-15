#' Constrained Hierarchical Agglomerative Clustering of Genomic regions(loci)
#' 
#' Function to perform adjacency-constrained hierarchical agglomerative clustering 
#' of genomic regions(loci)
#' 
#' \code{hicClust} performs constrained hierarchichal agglomerative clustering of 
#' genomic regions (loci) according to information provided by high-throughput conformation capture
#'  data (Hi-C). Constrained Hierarchical Agglomerative Clustering is hierarchical agglomerative 
#'  clustering in which each observation is associated to a position, and the clustering is 
#'  constrained so as only adjacent clusters are merged.
#' 
#' @param x either:
#' 1. A pxp contact map of class Matrix::dsCMatrix in which the entries are the number of counts of physical interactions observed between all pairs of loci
#' 2. An object of class HiTC::HTCexp. The corresponding Hi-C data is stored as a Matrix::dsCMatrix object in the intdata slot
#' 3. A text file with one line per pair of loci for which an interaction has been observed (in the format: locus1<tab>locus2<tab>signal). 
#' 
#' @param h band width. If not provided, `h` is set to default value `p-1`.It is assumed that the similarity between two items is 0 when these 
#' items are at a distance of more than band width h
#' 
#' @param \dots further arguments to be passed to \code{\link{read.table}} function. If not provided, the text file is supposed to be separated by tabulations, with no header.
#'  
#' @return Function \code{hicClust} returns an object of class \code{\link[stats]{hclust}}.  
#' 
#' @examples
#' #Input as HiTC::HTCexp object
#' data("hic_imr90_40_XX", package="adjclust")
#' 
#' #Input as HiTC::HTCexp object
#' res1 <- hicClust(hic_imr90_40_XX)
#' 
#' \dontrun{
#' #Input as Matrix::dsCMatrix contact map
#' mat <- HiTC::intdata(hic_imr90_40_XX) 
#' res2 <- hicClust(mat)
#' } 
#'
#' #Input as text file
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
  if( (class != "dsCMatrix")&&(class !=  "HTCexp")&&(!file.exists(x)) )
    stop("Invalid Input:x should be a text file or an object of class Matrix::dsCMatrix/HiTC::HTCexp")
  
  if(class == "dsCMatrix" || class  == "HTCexp") {
  
    if(class == "HTCexp")
      x <- intdata(x)
    
    p <- x@Dim[1]
    if(is.null(h)) h <- p-1  
    res <- adjClust(x, type = "similarity", h)
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
  res <- adjClust(m, type = "similarity", h)
  return(res)
  
  }
}
