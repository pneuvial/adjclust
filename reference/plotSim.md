# Plot (dis)similarity matrix

Heatmap of the (dis)similarity matrix

## Usage

``` r
plotSim(
  mat,
  type = c("similarity", "dissimilarity"),
  clustering = NULL,
  dendro = NULL,
  k = NULL,
  log = TRUE,
  legendName = "intensity",
  main = NULL,
  priorCount = 0.5,
  stats = c("R.squared", "D.prime"),
  h = NULL,
  axis = FALSE,
  naxis = min(10, nrow(mat)),
  axistext = NULL,
  xlab = "objects",
  cluster_col = "darkred",
  mode = c("standard", "corrected", "total-disp", "within-disp", "average-disp")
)
```

## Arguments

- mat:

  matrix to plot. It can be of class `'matrix'`, `'dgCMatrix'`,
  `'dsCMatrix'`, `'dist'`, `'HTCexp'`, `'snpMatrix'`.

- type:

  input matrix type. Can be either `"similarity"` or `"dissimilarity"`
  (kernels are supposed to be of type `"similarity"`).

- clustering:

  vector of clusters to display on the matrix (if not `NULL`). If
  `clustering` is provided, it must be a numeric vector of length
  identical to the matrix size with clusters identified as consecutive
  integers 1, 2, 3, ...

- dendro:

  dendrogram provided as an `hclust` object, or as another type of
  object that can be converted to `hclust` with `as.hclust`.

- k:

  number of clusters to display. Used only when `dendro` is not `NULL`
  and `clustering` is `NULL`. The clustering is then deduced from the
  dendrogram by a standard cut.

- log:

  logical. Should the breaks be based on log-scaled values of the matrix
  entries. Default to `TRUE`.

- legendName:

  character. Title of the legend. Default to `"intensity"`.

- main:

  character. Title of the plot. Default to `NULL` (no title).

- priorCount:

  numeric. Average count to be added to each entry of the matrix to
  avoid taking log of zero. Used only if `log == TRUE` and default to
  0.5.

- stats:

  input SNP correlation type. Used when `mat` is of type `SnpMatrix`.

- h:

  positive integer. Threshold distance for SNP correlation computation.
  Used when `mat` is of type `SnpMatrix`.

- axis:

  logical. Should x-axis be displayed on the plot? Default to `FALSE`.

- naxis:

  integer. If `axis == TRUE`, number of ticks to display on the x-axis.

- axistext:

  character vector. If `axis == TRUE`, labels to display of the x-axis
  (its length has to be equal to `naxis`).

- xlab:

  character. If `axis == TRUE`, x-axis title.

- cluster_col:

  color for the cluster line if `clustering` is not `NULL`.

- mode:

  type of dendrogram to plot (see
  [`plot.chac`](https://pneuvial.github.io/adjclust/reference/chac.md)).
  Default to `"standard"`.

## Details

This function produces a heatmap for the used (dis)similarity matrix
that can be used as a diagnostic plot to check the consistency between
the obtained clustering and the original (dis)similarity

## See also

[`select`](https://pneuvial.github.io/adjclust/reference/select.md),
[`adjClust`](https://pneuvial.github.io/adjclust/reference/adjClust.md)

## Examples

``` r
if (FALSE) { # \dontrun{
clustering <- rep(1:3, each = 50)
dist_data <- as.matrix(dist(iris[, 1:4]))
dendro_iris <- adjClust(dist_data, type = "dissimilarity")
plotSim(dist_data, type = "dissimilarity", dendro = dendro_iris, axis = TRUE)
plotSim(dist_data, type = "dissimilarity", dendro = dendro_iris,
        clustering = clustering)
plotSim(dist_data, type = "dissimilarity", dendro = dendro_iris, axis = TRUE,
        k = 3)
plotSim(dist_data, type = "dissimilarity", legendName = "IF", axis = TRUE, 
        clustering = clustering)
p <- plotSim(dist(iris[, 1:4]), type = "dissimilarity", log = FALSE, 
             clustering = clustering, cluster_col = "blue")
# custom palette
p + scale_fill_gradient(low = "yellow", high = "red")
# dsCMatrix
m <- Matrix(c(0, 0, 2, 0, 3, 0, 2, 0, 0), ncol = 3)
res <- adjClust(m)
plotSim(m, axis = TRUE)
plotSim(m, dendro = res)
# dgCMatrix
m <- as(m, "generalMatrix")
plotSim(m)
m <- as.dist(m)
if (require("HiTC", quietly = TRUE)) {
  load(system.file("extdata", "hic_imr90_40_XX.rda", package = "adjclust"))
  res <- hicClust(hic_imr90_40_XX, log = TRUE)
  plotSim(hic_imr90_40_XX, axis = TRUE)
}
if (requireNamespace("snpStats", quietly = TRUE)) {
  data(testdata, package = "snpStats")
  plotSim(Autosomes[1:200, 1:5], h = 3, stats = "R.squared", axis = TRUE,
          axistext = c("A", "B", "C", "D", "E"))
}
} # }
```
