source("https://bioconductor.org/biocLite.R")
biocLite("HiCDataHumanIMR90")
require(HiCDataHumanIMR90)
data(Dixon2012_IMR90)
obj <- hic_imr90_40$chrXchrX
save(obj, file="data/hic_imr90_40$chrXchrX.rda", compress="xz")
