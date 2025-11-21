# Adjacency-constrained Clustering of Hi-C contact maps

Adjacency-constrained hierarchical agglomerative clustering of Hi-C
contact maps

## Usage

``` r
hicClust(x, h = NULL, log = FALSE, ...)
```

## Arguments

- x:

  either: 1. A pxp contact sparse or dense matrix (classes matrix,
  Matrix, dscMatrix, dgTMatrix, dgCMatrix, dgeMatrix). Its entries are
  the number of counts of physical interactions observed between all
  pairs of loci. 2. An object of class HiTC::HTCexp. The corresponding
  Hi-C data is stored as a Matrix::dsCMatrix object in the intdata
  slot. 3. A text file path with one line per pair of loci for which an
  interaction has been observed (in the format:
  locus1\<tab\>locus2\<tab\>signal) or a matrix or data frame with
  similar data (3 columns).

- h:

  band width. If not provided, `h` is set to default value \`p-1\`.

- log:

  logical. Whether to log-transform the count data. Default to `FALSE`.

- ...:

  further arguments to be passed to
  [`read.table`](https://rdrr.io/r/utils/read.table.html) function when
  `x` is a text file name. If not provided, the text file is supposed to
  be separated by tabulations, with no header.

## Value

An object of class
[`chac`](https://pneuvial.github.io/adjclust/reference/chac.md).

## Details

Adjacency-constrained hierarchical agglomerative clustering (HAC) is HAC
in which each observation is associated to a position, and the
clustering is constrained so as only adjacent clusters are merged.
Genomic regions (loci) are clustered according to information provided
by high-throughput conformation capture data (Hi-C).

## References

Ambroise C., Dehman A., Neuvial P., Rigaill G., and Vialaneix N (2019).
*Adjacency-constrained hierarchical clustering of a band similarity
matrix with application to genomics*, Algorithms for Molecular Biology
14(22)"

Servant N. *et al* (2012). *HiTC : Exploration of High-Throughput 'C'
experiments. Bioinformatics*.

## See also

[`adjClust`](https://pneuvial.github.io/adjclust/reference/adjClust.md)

## Examples

``` r
# input as HiTC::HTCexp object
if (FALSE) { # \dontrun{
if (require("HiTC", quietly = TRUE)) {
  load(system.file("extdata", "hic_imr90_40_XX.rda", package = "adjclust"))
  res1 <- hicClust(hic_imr90_40_XX)
}
} # }

# input as Matrix::dsCMatrix contact map
if (FALSE) { # \dontrun{
mat <- HiTC::intdata(hic_imr90_40_XX) 
res2 <- hicClust(mat)
} # }

# input as text file
res3 <- hicClust(system.file("extdata", "sample.txt", package = "adjclust"))
#> Note: 1 merges with non increasing heights.
```
