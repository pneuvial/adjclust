library("adjclust") #primitive version 0.3.0
library("snpStats")

data("ld.example", package="snpStats")
p <- ncol(ceph.1mb)
nSamples <- nrow(ceph.1mb)
h <- 100
ceph.1mb[4,286]@.Data[1,1] <- as.raw(3) #to avoid NaNs

ld.ceph <- ld(ceph.1mb, stats="R.squared", depth=h)
ld.ceph <- round(ld.ceph, digits=10)

prevfit <- adjclust:::adjClustBand_heap(ld.ceph@x, p, h, blMin=1)

save(prevfit, file="data/prevfit.rda", compress="xz")
