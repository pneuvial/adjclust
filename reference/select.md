# Clustering selection

Clustering selection from a chac object with the slope heuristic or the
broken stick heuristic

## Usage

``` r
select(
  x,
  type = c("capushe", "bstick"),
  k.max = NULL,
  graph = FALSE,
  pct = 0.15
)
```

## Arguments

- x:

  an object of class 'chac'

- type:

  model selection approach between slope heuristic (`"capushe"`) and
  broken stick approach (`"bstick"`)

- k.max:

  maximum number of clusters that can be selected. Default to `NULL`, in
  which case it is set to \\\min(\max(100, \frac{n}{\log(n)}),
  \frac{n}{2})\\ where \\n\\ is the number of objects to be clustered
  for capushe and to \\n\\ for the broken stick model

- graph:

  logical. Whether the diagnostic plot for the capushe selection is
  displayed or not. Default to `FALSE`

- pct:

  minimum percentage of points for the plateau selection in capushe
  selection. See [`DDSE`](https://rdrr.io/pkg/capushe/man/DDSE.html) for
  further details

## Value

The function returns the clustering selected by the slope heuristic, as
implemented in the R package `capushe`.

## References

Baudry J.P., Maugis C., and Michel B. (2012). Slope heuristics: overview
and implementation. *Statistics and Computing*, **22**(2), 355-470.
MacArthur, R.H. (1957) On the relative abundance of bird species.
*Proceedings of the National Academy of Sciences*, **43**, 293-295.

## Examples

``` r
if (FALSE) if (require("HiTC", quietly = TRUE)) {
  load(system.file("extdata", "hic_imr90_40_XX.rda", package = "adjclust"))
  res <- hicClust(hic_imr90_40_XX, log = TRUE)
  selected.capushe <- select(res)
  table(selected.capushe)
  selected.bs <- select(res, type = "bstick")
  table(selected.bs)
} # \dontrun{}

res <- adjClust(dist(iris[, 1:4]))
#> Note: input class is 'dist' so 'type' is supposed to be 'dissimilarity'
#> Note: 30 merges with non increasing heights.
select.clust <- select(res, "bs")
table(select.clust)
#> select.clust
#>  1  2  3 
#> 50 50 50 
```
