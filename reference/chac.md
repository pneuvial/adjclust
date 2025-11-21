# Class chac

S3 class for Constrained Hierarchical Agglomerative Clustering results

## Usage

``` r
# S3 method for class 'chac'
as.hclust(x, ...)

# S3 method for class 'chac'
print(x, ...)

# S3 method for class 'chac'
head(x, ...)

# S3 method for class 'chac'
summary(object, ...)

# S3 method for class 'chac'
plot(
  x,
  y,
  ...,
  mode = c("standard", "corrected", "total-disp", "within-disp", "average-disp"),
  nodeLabel = FALSE
)

diagnose(x, graph = TRUE, verbose = TRUE)

correct(x)

cutree_chac(tree, k = NULL, h = NULL)
```

## Arguments

- x, object, tree:

  an object of class 'chac'

- ...:

  for [`plot`](https://rdrr.io/r/graphics/plot.default.html), arguments
  passed to the function
  [`plot.dendrogram`](https://rdrr.io/r/stats/dendrogram.html). Default
  values for `type` and `leaflab` are respectively set to `"triangle"`
  and `"none"`

- y:

  not used

- mode:

  type of dendrogram to plot (see Details). Default to `"standard"`

- nodeLabel:

  (logical) whether the order of merging has to be displayed or not.
  `nodeLabel=TRUE` prints orders of fusion at corresponding nodes.
  Default to `FALSE`

- graph:

  (logical) whether the diagnostic plot has to be displayed or not.
  Default to `TRUE`

- verbose:

  (logical) whether to print a summary of the result or not. Default to
  `TRUE`

- k:

  an integer scalar or vector with the desired number of groups

- h:

  numeric scalar or vector with heights where the tree should be cut.
  Only available when the heights are increasing

## Value

The function `plot.chac` displays the dendrogram and additionally
invisibly returns an object of class
[`dendrogram`](https://rdrr.io/r/stats/dendrogram.html) with heights as
specified by the user through the option `mode`.

`diagnose` invisibly exports a data frame with the numbers of decreasing
merges described by the labels of the clusters being merged at this step
and at the previous one, as well as the corresponding merge heights.

The function `correct` returns a `chac` objects with modified heights so
as they are increasing. The new heights are calculated in an way
identical to the option `mode = "corrected"` of the function `plot.chac`
(see Details). In addition, the `chac` object has its field `method`
modified from `adjClust` to `adjClust-modified`.

The function `cutree_chac` returns the clustering with `k` groups or
with the groups obtained by cutting the tree at height `h`. If the
heights are not increasing, the cutting of the tree is based on the
corrected heights as provided by the function `correct`.

## Details

Methods for class 'chac'

When `plot.chac` is called with `mode = "standard"`, the standard
dendrogram is plotted, even though, due to contingency constrains, some
branches are reversed (decreasing merges). When `plot.chac` is called
with `mode = "corrected"`, a correction is applied to original heights
so as to have only non decreasing merges). It does not change the result
of the clustering, only the look of the dendrogram for easier
interpretation.  
  
Other modes are provided that correspond to different alternatives
described in Grimm (1987):

- in `mode = "within-disp"`, heights correspond to within-cluster
  dispersion, *i.e.*, for a corresponding cluster, its height is
  \$\$I(C) = \sum\_{i \in C} d(i,g_C)\$\$ where \\d\\ is the
  dissimilarity used to cluster objects and \\g_C\\ is the center of
  gravity of cluster \\C\\. In this case, heights are always non
  decreasing;

- in `mode = "total-disp"`, heights correspond to the total
  within-cluster dispersion. It is obtained from `mode = "standard"` by
  the cumulative sum of its heights. In this case, heights are always
  non decreasing;

- in `mode = "average-disp"`, heights correspond to the within-cluster
  dispersion divided by the cluster size. In this case, there is no
  guaranty that the heights are non decreasing. When reversals are
  detected, a warning is printed to advice the user to change the mode
  of the representation.

Grimm (1987) indicates that heights as provided by
`mode = "within-disp"` are highly dependent on cluster sizes and that
the most advisable representation is the one provided by
`mode = "total-disp"`. Further details are provided in the vignette
"Notes on CHAC implementation in adjclust".

## References

Grimm, E.C. (1987) CONISS: a fortran 77 program for stratigraphically
constrained analysis by the method of incremental sum of squares.
*Computer & Geosciences*, **13**(1), 13-35.
