# adjclust

adjclust is a package that provides methods to perform adjacency-constrained hierarchical agglomerative clustering. Adjacency-constrained hierarchical agglomerative clustering is hierarchical agglomerative clustering in which each observation is associated to a position, and the clustering is constrained so as only adjacent clusters are merged. It is a common method used widely in various application fields including ecology (Quaternary data) and bioinformatics (for instance in Genome Wide Association Studies).

<em> Version 0.4.0 of this package was completed as a part of the [Google Summer of Code 2017](https://summerofcode.withgoogle.com/projects/#4961904920363008) program.</em>

Present version adjclust package provides three user level functions: `adjClust`, `snpClust` and `hicClust`.

`adjClust` function performs adjacency-constrained hierarchichal agglomerative clustering for standard and sparse, similarity and dissimilarity matrices and dist objects.Matrix::dgCMatrix and Matrix::dsCMatrix are the supported sparse matrix classes. Lets look at an example

```r
library("adjclust")

sim <- matrix(c(1,0.1,0.2,0.3,0.1,1,0.4,0.5,0.2,0.4,1,0.6,0.3,0.5,0.6,1), nrow=4)
h <- 3
fit1 <- adjClust(sim, "similarity", h, 1, FALSE)
plot(fit1)

dist <- as.dist(sqrt(2-(2*sim)))

#Compatibility with dist objects
fit2 <- adjClust(dist, "dissimilarity", h, 1, FALSE)
plot(fit2)

```


`snpClust` function performs adjacency-constrained hierarchichal agglomerative clustering for specific application of Genome Wide Association Studies (GWAS). See [GWAS Vignette](vignettes/snpClust.Rmd) for details.

`hicClust`function performs adjacency-constrained hierarchichal agglomerative clustering for specific application of Hi-C data analysis. See [Hi-C Vignette](vignettes/hicClust.Rmd) for details.
