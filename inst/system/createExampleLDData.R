library("snpStats")

data(ld.example, package="snpStats")
str(ceph.1mb)

ld.ceph <- ld(ceph.1mb, stats=c("D.prime", "R.squared"), depth=100)

R2.100 <- ld.ceph$R.squared
Dprime.100 <- ld.ceph$D.prime

save(R2.100, file="data/R2.100.rda", compress="xz")
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


save(Dprime.100, file="data/Dprime.100.rda", compress="xz")
#' Linkage disequilibrium on chromosome 22 in the European population
#'
#' A dataset containing D' linkage disequilibrium (LD) statistics for
#' \eqn{p=603} SNPs spanning a one megabase regions on chromosome 22, in a
#' sample of 90 Europeans (CEPH). LD values were only calculated for entries in
#' a diagonal band of size \code{h=100}; the other values are assumed to be 0.
#'
#' @format A sparse matrix of class \code{\link[Matrix]{dgCMatrix-class}}.
#'
#' @source Data extracted from \code{\link[snpStats]{ld.example}}. The original data is
#'   from the HapMap project.
#'
