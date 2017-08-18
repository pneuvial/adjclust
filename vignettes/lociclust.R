## ---- message = FALSE----------------------------------------------------
library("adjclust")
library("HiTC")

## ------------------------------------------------------------------------
data("hic_imr90_40$chrXchrX", package="adjclust")
mat <- intdata(obj)
image(mat, lwd=0)

## ------------------------------------------------------------------------
  h <- 3881
  class(obj)
  fit1 <- lociclust(obj, h)

## ------------------------------------------------------------------------
  class(mat)
  fit2 <- lociclust(mat, h)

## ------------------------------------------------------------------------
  V3 <- mat@x
  V1 <- mat@Dimnames[[1]][mat@i+1]          #loci1names
  V2 <- rep(mat@Dimnames[[2]], diff(mat@p)) #loci2names
  
  l1 <- unique(c(as.numeric(V1), as.numeric(V2)))  
  
  l2 <- setdiff(c(72037:75918), l1) 
  V1 <- append(V1, l2)
  V2 <- append(V2, l2)
  V3 <- append(V3, rep(0, length(l2)))
  
  content <- cbind(as.numeric(V1), as.numeric(V2), as.numeric(V3))
  
  tf <- tempfile(fileext = ".txt")
  write.table(content, tf, sep = " ", col.names = FALSE, row.names = FALSE)

## ------------------------------------------------------------------------
  fit3 <- lociclust(tf, h)  

## ------------------------------------------------------------------------
head(cbind(fit1$merge, fit1$gains))

## ------------------------------------------------------------------------
all.equal(fit1$merge, fit2$merge)
all.equal(fit1$gains, fit2$gains)

all.equal(fit2$merge, fit3$merge)
all.equal(fit2$gains, fit3$gains)

## ------------------------------------------------------------------------
sessionInfo()

