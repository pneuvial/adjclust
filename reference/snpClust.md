# Adjacency-constrained Clustering of Single Nucleotide Polymorphisms

Adjacency-constrained hierarchical agglomerative clustering of Single
Nucleotide Polymorphisms based on Linkage Disequilibrium

## Usage

``` r
snpClust(x, h = ncol(x) - 1, stats = c("R.squared", "D.prime"))
```

## Arguments

- x:

  either a genotype matrix of class
  [`SnpMatrix`](https://rdrr.io/pkg/snpStats/man/SnpMatrix-class.html)/[`matrix`](https://rdrr.io/r/base/matrix.html)
  or a linkage disequilibrium matrix of class
  [`dgCMatrix`](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html). In
  the latter case the LD values are expected to be in \[0,1\]

- h:

  band width. If not provided, `h` is set to default value \`p-1\` where
  \`p\` is the number of columns of `x`

- stats:

  a character vector specifying the linkage disequilibrium measures to
  be calculated (using the
  [`ld`](https://rdrr.io/pkg/snpStats/man/ld.html) function) when `x` is
  a genotype matrix. Only "R.squared" and "D.prime" are allowed, see
  Details.

## Value

An object of class
[`chac`](https://pneuvial.github.io/adjclust/reference/chac.md) (when no
LD value is missing)

## Details

Adjacency-constrained hierarchical agglomerative clustering (HAC) is HAC
in which each observation is associated to a position, and the
clustering is constrained so as only adjacent clusters are merged. SNPs
are clustered based on their similarity as measured by the linkage
disequilibrium.

In the special case where genotypes are given as input and the
corresponding LD matrix has missing entries, the clustering cannot be
performed. This can typically happen when there is insufficient
variability in the sample genotypes. In this special case, the indices
of the SNP pairs which yield missing values are returned.

If `x` is of class
[`SnpMatrix`](https://rdrr.io/pkg/snpStats/man/SnpMatrix-class.html) or
[`matrix`](https://rdrr.io/r/base/matrix.html), it is assumed to be a
\\n \times p\\ matrix of \\p\\ genotypes for \\n\\ individuals. This
input is converted to a LD similarity matrix using the
[`snpStats::ld`](https://rdrr.io/pkg/snpStats/man/ld.html). If `x` is of
class
[`dgCMatrix`](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html), it
is assumed to be a (squared) LD matrix.

Clustering on a LD similarity other than "R.squared" or "D.prime" can be
performed by providing the LD values directly as argument `x`. These
values are expected to be in \[0,1\], otherwise they are truncated to
\[0,1\].

## References

Dehman A. (2015) *Spatial Clustering of Linkage Disequilibrium Blocks
for Genome-Wide Association Studies*, PhD thesis, Universite Paris
Saclay.

Dehman, A. Ambroise, C. and Neuvial, P. (2015). Performance of a
blockwise approach in variable selection using linkage disequilibrium
information. \*BMC Bioinformatics\* 16:148.

Ambroise C., Dehman A., Neuvial P., Rigaill G., and Vialaneix N (2019).
*Adjacency-constrained hierarchical clustering of a band similarity
matrix with application to genomics*, Algorithms for Molecular Biology
14(22)"

## See also

[`adjClust`](https://pneuvial.github.io/adjclust/reference/adjClust.md)
[`ld`](https://rdrr.io/pkg/snpStats/man/ld.html)

## Examples

``` r
## a very small example
if (requireNamespace("snpStats", quietly = TRUE)) {
  data(testdata, package = "snpStats")

  # input as snpStats::SnpMatrix
  fit1 <- snpClust(Autosomes[1:200, 1:5], h = 3, stats = "R.squared")

  # input as base::matrix
  fit2 <- snpClust(as.matrix(Autosomes[1:200, 1:5]), h = 3, stats = "R.squared")

  # input as Matrix::dgCMatrix
  ldres <- snpStats::ld(Autosomes[1:200, 1:5], depth = 3, stats = "R.squared", symmetric = TRUE)
  fit3 <- snpClust(ldres, 3)
}
#> Note: 1 merges with non increasing heights.
#> Note: 1 merges with non increasing heights.
#> Note: forcing the diagonal of the LD similarity matrix to be 1
#> Note: 1 merges with non increasing heights.
```
