## ---- message = FALSE----------------------------------------------------
library("adjclust")
library("HiTC")

## ------------------------------------------------------------------------
data("hic_imr90_40$chrXchrX", package="adjclust")

## ------------------------------------------------------------------------
selected <- apply(intdata(obj), 1, sum) > 0
obj <- new("HTCexp", intdata(obj)[selected,selected], x_intervals(obj)[selected,], y_intervals(obj)[selected, ])

## ------------------------------------------------------------------------
mat <- intdata(obj)
image(mat, lwd=0)

## ------------------------------------------------------------------------
  class(obj)
  fit1 <- hicclust(obj)

## ------------------------------------------------------------------------
  class(mat)
  fit2 <- hicclust(mat)

## ------------------------------------------------------------------------
  V3 <- mat@x
  V1 <- mat@Dimnames[[1]][mat@i+1]          #loci1names
  V2 <- rep(mat@Dimnames[[2]], diff(mat@p)) #loci2names
  
  content <- cbind(as.numeric(V1), as.numeric(V2), as.numeric(V3))
  
  tf <- tempfile(fileext = ".txt")
  write.table(content, tf, sep = " ", col.names = FALSE, row.names = FALSE)

## ------------------------------------------------------------------------
  fit3 <- hicclust(tf, sep = " ")  

## ------------------------------------------------------------------------
head(cbind(fit1$merge, fit1$gains))

## ------------------------------------------------------------------------
all.equal(fit1$merge, fit2$merge)
all.equal(fit1$gains, fit2$gains)

all.equal(fit2$merge, fit3$merge)
all.equal(fit2$gains, fit3$gains)

## ------------------------------------------------------------------------
sessionInfo()

