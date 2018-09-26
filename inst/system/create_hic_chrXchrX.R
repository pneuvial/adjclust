library("HiTC")
## source("https://bioconductor.org/biocLite.R")
## biocLite("HiCDataHumanIMR90")
library("HiCDataHumanIMR90")
data("Dixon2012_IMR90")

## `hic_imr90_40` is a list of objects of class `HTCexp` which has been obtained from the `HiTC` package. 
## Each of these objects corresponds to a Hi-C contact map between one chromosome and another. 
## 
## 
obj <- hic_imr90_40$chrXchrX ## contact map corresponding to chromosome X vs chromosome X. 
## The corresponding Hi-C data is stored as a Matrix::dsCMatrix in the slot named intdata.

## Remove all the rows and columns containing only zeros from the dataset.
selected <- apply(intdata(obj), 1, sum) > 0

## keep only the first 500 indices (10x smaller object)
idxs <- which(selected)[1:500]
intd <- intdata(obj)[idxs, idxs]
x_int <- x_intervals(obj)[idxs, ]
y_int <- y_intervals(obj)[idxs, ]
hic_imr90_40_XX <- new("HTCexp", intd, x_int, y_int)

save(hic_imr90_40_XX, file = "inst/extdata/hic_imr90_40_XX.rda", compress = "xz")
