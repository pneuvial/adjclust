# adjclust

Adjacency-constrained clustering of a block-diagonal similarity matrix

## Low-level clustering function

library("adjclust")
data("R2.100", package="adjclust")
str(R2.100)
h <- 100
p <- 603

## default flavor: "crayons"
res <- adjClustBand(R2.100, p, h)

## alternative flavor: "PseudoMatrix"
resP <- adjClustBand(R2.100, p, h, flavor="PseudoMatrix")
