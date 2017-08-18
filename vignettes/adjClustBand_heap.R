## ------------------------------------------------------------------------
library("adjclust")

## ---- results="hide", message=FALSE--------------------------------------
library("matrixStats")
library("snpStats")
data("ld.example", package="snpStats")

## ---- echo=FALSE---------------------------------------------------------
p <- ncol(ceph.1mb)
nSamples <- nrow(ceph.1mb)

## ------------------------------------------------------------------------
ceph.1mb

## ------------------------------------------------------------------------
ld.ceph <- ld(ceph.1mb, stats="R.squared", depth=p-1)
image(ld.ceph, lwd=0)

## ------------------------------------------------------------------------
h <- 100
ld.ceph <- ld(ceph.1mb, stats="R.squared", depth=h)
image(ld.ceph, lwd=0)

## ------------------------------------------------------------------------
any(is.na(ld.ceph))

#which(is.na(ld.ceph), arr.ind = TRUE)

ceph.1mb[4,286]@.Data[1,1] <- as.raw(3)  ## to avoid NaNs
ld.ceph <- ld(ceph.1mb, stats="R.squared", depth=h)  ##calculate LD values again after NaN correction
ld.ceph <- round(ld.ceph, digits=10)  ## precision above 1e-10 is not required for this example
any(is.na(ld.ceph))
diag(ld.ceph) <- 1

## ------------------------------------------------------------------------
fit <- adjClustBand_heap(ld.ceph, "similarity", h, blMin=1)

## ------------------------------------------------------------------------
plot(fit)

## ------------------------------------------------------------------------
head(cbind(fit$merge, fit$gains))

## ------------------------------------------------------------------------
sessionInfo()

