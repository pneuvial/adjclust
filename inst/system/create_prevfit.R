if (packageVersion("adjclust") > "0.3.0") {
    stop("Please install a version of 'adjclust' <= 0.3.0 to run this script")
}
library("adjclust")
library("snpStats")

data("ld.example", package="snpStats")
p <- ncol(ceph.1mb)
nSamples <- nrow(ceph.1mb)
h <- 100
ceph.1mb[4,286]@.Data[1,1] <- as.raw(3) #to avoid NaNs

ld.ceph <- ld(ceph.1mb, stats="R.squared", depth=h)
ld.ceph <- round(ld.ceph, digits=10)

res_adjclust_0.3.0 <- adjclust:::adjClustBand_heap(ld.ceph@x, p, h, blMin=1)
save(res_adjclust_0.3.0, file="data/res_adjclust_0.3.0.rda", compress="xz")
