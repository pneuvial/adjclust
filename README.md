[![Travis Build Status](https://travis-ci.org/pneuvial/adjclust.svg?branch=develop)](https://travis-ci.org/pneuvial/adjclust)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/pneuvial/adjclust?branch=develop&svg=true)](https://ci.appveyor.com/project/pneuvial/adjclust)
[![Coverage Status](https://img.shields.io/codecov/c/github/pneuvial/adjclust/develop.svg)](https://codecov.io/github/pneuvial/adjclust?branch=develop)

# adjclust

adjclust is a package that provides methods to perform adjacency-constrained hierarchical agglomerative clustering. Adjacency-constrained hierarchical agglomerative clustering is hierarchical agglomerative clustering in which each observation is associated to a position, and the clustering is constrained so as only adjacent clusters are merged. It is a common method used widely in various application fields including ecology (Quaternary data) and bioinformatics (for instance in Genome Wide Association Studies).

<em> Version 0.4.0 of this package was completed as a part of the [Google Summer of Code 2017](https://summerofcode.withgoogle.com/projects/#4961904920363008) program.</em>

Present version adjclust package provides three user level functions: `adjClust`, `snpClust` and `hicClust`.

`adjClust` function performs adjacency-constrained hierarchichal agglomerative clustering for standard and sparse, similarity and dissimilarity matrices and dist objects.Matrix::dgCMatrix and Matrix::dsCMatrix are the supported sparse matrix classes. See the [documentation](man/adjClust.Rd) for details.

`snpClust` function performs adjacency-constrained hierarchichal agglomerative clustering for specific application of Genome Wide Association Studies (GWAS). See [GWAS Vignette](vignettes/snpClust.Rmd) for details.

`hicClust`function performs adjacency-constrained hierarchichal agglomerative clustering for specific application of Hi-C data analysis. See [Hi-C Vignette](vignettes/hicClust.Rmd) for details.
