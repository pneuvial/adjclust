# adjclust

Adjacency-constrained clustering of a block-diagonal similarity matrix

## Low-level clustering function

data("ld_ceph", package="adjclust")

x <- R2.100
x <- Dprime.100

h <- 100
p <- 603

## default flavor: "crayons"
res <- adjClustBand(x, p, h)

## alternative flavor: "pseudoMatrix"
resP <- adjClustBand(x, p, h, flavor="pseudoMatrix")
