# Version 0.6.9  [2024-10-07]

* Bug fix on OMP option (#81)

# Version 0.6.9  [2024-02-07]

* Properly handled OMP threads in C++ code (now default to 1 but with an option
to increase this value)
* Reintroduced tests and examples
* Removed WCSS function that was not exported or documented
* Fixed a problem in S3class for a non exported function

# Version 0.6.8  [2024-01-10]

* Fix CRAN error on useNames (deprecated NA)
* Fix CRAN note on itemize (unecessary use of itemize)
* Limited OMP threads to 2 in examples, vignettes and tests
* Updated citation of the package

# Version 0.6.7  [2023-04-24]

* Fix #60 (increase test coverage)
* Fix #61 (NAMESPACE error in 'SnpClust.matrix')
* Fix #45 (update 'calls' for 'adjclust.*' methods)
* Fix #49 (calls to 'library' in tests)
* Fix #55 (pkgdown action)
* Fix #40 (moved plotSim to a ggplot2 version)
* Fix #30 (improved plotSim, k and correct have been added as new arguments)
* Fix #62 (citation file corrected)
* Fix #68 (CRAN note)

# Version 0.6.6 [2022-09-13]

* Minor updates in tests to comply with updates in Matrix 1.5.0

# Version 0.6.5 [2022-08-18]

* Code improvement: replaced `class( ) !=` by equivalent expression using `inherits( )`

# Version 0.6.4 [2022-03-30]

* Update URL in DESCRIPTION
* Fixed bug in tests due to recent upgrade of Matrix
* Slightly improved vignettes and references here and there

# Version 0.6.3 [2021-07-26]

* Fix issues for CRAN submission: RcppArmadillo moved to LinkingTo, version number
* Cleanup: move res_adjclust_0.3.0.rda from data/ to inst/extdata/
* Cleanup: remove useless data entries from pkgdown index, and add missing ones
* Update docs (using pkgdown)
* Move site to dedicated branch (release with `pkgdown::deploy_to_branch()`)
* Add inst/CITATION and add Randriamihamison et al 2020

# Version 0.6.2 [2021-07-22]

* Increase in speed by code optimization and by using RcppArmadillo (update by Gabriel Hoffman), at least for linux machines (uses OpenMP)
* reduces memory usage (update by Gabriel Hoffman)
* option to disable expensive checking code (strictCheck = TRUE by default, update by Gabriel Hoffman)

# Version 0.5.99 [2020-06-08]

* Updated citation in DESCRIPTION and man files (almob paper) and added a CITATION file.
* removed exportation of S3 classes.
* changed cutree into a simple function 'cutree_chac' rather than a method (because stats::cutree is not a method).

# Version 0.5.9 [2019-12-10]

* Clarified types of inputs handled by adjclust ('S3 methods').
* Shortened some examples.
* Rewrote plotSim to avoid CRAN error on devel (chose to use 'S3 methods' now).
* Fixed a similar problem in helpers.
* Clarified code and comments.
* Fixed Issue #35 (probably due to a change in the upstream snpStats package).
* Added a test on plotSim.
* Changed test on NA in snpStats package (probably due to a change in the upstream snpStats package) and fixed Issue #43.
* Minor updates to tests.
* The package passes R CMD check with [0;0;0].

# Version 0.5.7 [2018-09-26]

* Example Hi-C data now 10x smaller (subset of the original one). The package 
is smaller and tests are faster.
* implemented a model selection approach based on slope heuristic or on the 
broken stick heuristic to select a relevant number of clusters
* fixed minor problems in some method definition for class 'chac'
* proposed a log-transformation of data in the wrapper 'hicClust'
* implemented a heatmap with possible highlighting of the constrained 
clustering
* implemented an option to display number of the merge on the dendrogram

# Version 0.5.6 [2018-02-08]

* changed dependencies to Bioconductor packages 'HiTC' and 'snpStats' into 
Suggest and conditionally used them

# Version 0.5.5 [2018-01-30]

* simplified code (replaced many C functions by a unique R function using 
Matrix)
* adjClust now properly handles similarities with diagonal entries different 
from 1
* removed arguments that were not used (blMin and verbose)
* simplified Hi-C example

# Version 0.5.4 [2018-01-12]

* More tests for modify and modifySparse 
* BUG FIX in condnCheck

# Version 0.5.3 [2017-12-04]

* 'height' is now defined as the value of the linkage criterion (as is done in
'hclust'), rather than the total inertia of the clustering (as is done in
'rioja').
* Added several representations for the dendrogram corresponding to different
choices for the height.
* Improved documentation and vignettes.
* Removed non-standard fields in the output of 'adjclust' (#13).
* Added tests for: equivalence with 'hclust',  comparing sum of heights and 
   pseudo inertia, plots, non-increasing heights, cutree (#14).
* Fixed #13 (man).
* Fixed #15 (cutree with decreasing merges).
* Fixed #3 (Non-positive 'gains').
* Using BiocStyle::html_document2 as a temporary fix for vignette 
  compilation errors.

# Version 0.5.2 [2017-10-17]

* Added citation to Alia Dehman's PhD thesis to DESCRIPTION.

# Version 0.5.1 [2017-10-16]

* More informative 'Description' of the method in DESCRIPTION
* Updates to test scripts to pass R CMD check on all windows platforms
* Moved README-*.png files to man/figures

# Version 0.5.0 [2017-10-13]

* Bump version number for CRAN submission

# Version 0.4.2 [2017-10-05]

* Added 'chac' S3 class and corresponding 'plot' and 'summary' methods
* Documentation cleanups
* Removed objects "R2.100" and "Dprime.100" (can be obtained from the 
  imported 'snpStats' package)
* In 'snpClust': argument 'stat' is now passed to the 'snpStats::ld' function 
  through '...'
* Some code cleanups
* Improved handling of default value for 'h' in 'adjclust' for 'dist' objects
* Renamed 'prevfit' into the more explicit 'res_adjclust_0.3.0'
* Dropped 'simmatrix' toy data set (now generated on the fly in tests)

# Version 0.4.1 [2017-09-15]

* Cleanups in Hi-C and LD vignettes and corresponding tests
* Dropped outdated BALD test script
* Added test script for NA values in LD
* Renamed Hi-C data sets and updated corresponding documentation
* Added package website generated by pkgdown

# Version 0.4.0 [2017-08-29]

* Implemented interface to handle standard and sparse matrices in adjClust
* Implemented interface to handle either kernel or dissimilarities
* Implemented wrapper for SNP and Hi-C data
* Documented the package and created vignettes for the different use cases
* Added scripts to increase package coverage and test the equivalence with 
  rioja for the small dimensional case
* Cleaned up code to improve efficiency and removed unnecessary scripts and functions

# Version 0.3.0 [2017-02-13]

* Removed 'adjClustBand': main entry points are now 'HeapHop' and 'adjClustBand_heap'.
* Updated test scripts and LD vignette accordingly.
* Added Travis CI and AppVeyor support.

# Version 0.2.*

## Version 0.2.3 [2017-02-02]

* Updated LD vignette
* In adjClustBand, renamed flavor "Koskas" to "PseudoMatrix"

## Version 0.2.2 [2016-12-01]

* Added dummy R/adjclust.R so that document() adds 'importFrom Rcpp evalCpp' to NAMESPACE
* "Fixed" warning at check due to .hpp file in src (this warning should not exist IMHO)

## Version 0.2.1 [2016-11-09]

* Added minimal documentation
* Replaced "std::cout" by "Rcpp::Rcout", and so on for "exit()" and "cerr".

## Version 0.2.0 [2016-06-24]

* Incorporated Michel's implementation (R function 'HeapHop')
* 'adjClustBand' is now a wrapper to call either Alia's or Michel's
  implementation

# Version 0.1.0 [2016-06-24]

* Created from BALD
* Added a test to check that we are reproducing the results of BALD::cWard