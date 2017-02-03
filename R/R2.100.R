#' Linkage disequilibrium on chromosome 22 in the European population
#'
#' A dataset containing R-squared linkage disequilibrium (LD) statistics for
#' \eqn{p=603} SNPs spanning a one megabase regions on chromosome 22, in a
#' sample of 90 Europeans (CEPH). LD values were only calculated for entries in
#' a diagonal band of size \code{h=100}; the other values are assumed to be 0.
#'
#' @format A sparse matrix of class \code{\link[Matrix]{dgCMatrix-class}}.
#'
#' @source Data extracted from \code{\link[snpStats]{ld.example}}. The original data is
#'   from the HapMap project.
#'
"R2.100"
