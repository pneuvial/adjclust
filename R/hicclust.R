#' Constrained Hierarchical Agglomerative Clustering of Genomic regions(loci)
#' 
#' Function to perform adjacency-constrained hierarchical agglomerative clustering 
#' of genomic regions(loci)
#' 
#' \code{hicclust} performes constrained hierarchichal agglomerative clustering of 
#' genomic regions(loci).Constrained Hierarchical Agglomerative Clustering is hierarchical
#'  agglomerative clustering in which each observation is associated to a position, and 
#'  the clustering is constrained so as only adjacent clusters are merged.
#' 
#' @param x x can be:
#' 1. A contact map matrix of class Matrix::dsCMatrix.A contact map is a p x p similarity
#' matrix, where the similarity between any pair of loci quantifies the intensity of the 
#' physical interaction between these loci
#' 2. An object of class HiTC::HTCexp.The corresponding Hi-C data should be stored as a 
#' Matrix::dsCMatrix in the intdata slot
#' 3. A text file with one line per pair of loci for which an interaction has been observed 
#' (in the format: locus1<space>locus2<space>signal). 
#' 
#' @param h band width. It is assumed that the similarity between two items is 0 when these 
#' items are at a distance of more than band width h
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
#' V3 <- c(201,300,203,103, 34, 21, 41, 22, 35, 66)
#' V1 <- c( 31, 32, 35, 36, 32, 33, 35, 34, 35, 35)          #loci1names
#' V2 <- c( 31, 32, 35, 36, 31, 31, 31, 32, 32, 34)          #loci2names
#' h <- 5
#' content <- cbind(V1, V2, V3)
#' tf <- tempfile(fileext = ".txt")
#' write.table(content, tf, sep = " ", col.names = FALSE, row.names = FALSE)
#' res3 <- hicclust(tf, h)
#' 
#' @export 
#' 
#' @importFrom utils read.table
#' @importFrom HiTC intdata

hicclust <- function(x, h) {

  if (!is.numeric(h))
    stop("h should be numeric")
  
  class <- class(x)
  if( (class != "dsCMatrix")&&(class !=  "HTCexp")&&(!file.exists(x)) )
    stop("Invalid Input:x should be a text file or an object of class Matrix::dsCMatrix/HiTC::HTCexp")
  
  if(class == "dsCMatrix" || class  == "HTCexp") {
  
    if(class == "HTCexp")
      x <- intdata(x)
    
    res <- adjClustBand_heap(x, type = "similarity", h)
    return(res)
  
  } else {
  
  x <- file(x, "r")
  df <- read.table(x, header = FALSE)
  close(x)
  lis <- sort(unique(c(df$V1, df$V2)))
  p <- length(lis)
  rowindx <- match(df$V1, lis)
  colindx <- match(df$V2, lis)
  
  m <- matrix(0, nrow = p, ncol = p)
  m[cbind(rowindx,colindx)] <- m[cbind(colindx,rowindx)] <- df$V3
  res <- adjClustBand_heap(m, type = "similarity", h)
  return(res)
  
  }
}
