check_snp <- function() {
  if (!require("snpStats")) {
    skip("'snpStats' package not available")
  }
}

context("Ascending compatibility of the adjclust algorithm")

test_that("snpClust gives results identical to those of adjclust 0.3.0", {
  check_snp()
  #Sys.setenv("OMP_THREAD_LIMIT" = 2)
  
  ## Note: this test depends on external data (genotypes) and functions 
  ## (snpStats::ld) which may change over time
  ## skip_on_cran()
  
  pathname <- system.file("extdata", "res_adjclust_0.3.0.rds", package = "adjclust")
  prevfit <- readRDS(pathname)

  data("ld.example", package = "snpStats")
  h <- 100
  
  ld.ceph <- snpStats::ld(ceph.1mb, depth = h, stats = "R.squared")
  
  p <- ncol(ceph.1mb)
  nSamples <- nrow(ceph.1mb)
  h <- 100
  ceph.1mb[4,286]@.Data[1,1] <- as.raw(3) ## to avoid NaNs
  ld.ceph <- snpStats::ld(ceph.1mb, depth = h, stats = "R.squared", symmetric = TRUE)
  ld.ceph <- round(ld.ceph, digits = 10)
  
  ## diagonal elements are 0 
  expect_identical(unname(diag(ld.ceph)), rep(0, p))
  expect_message(snpClust(ld.ceph, h = 100), 
                 "Note: forcing the diagonal of the LD similarity matrix to be 1",
                 all = FALSE)
  
  diag(ld.ceph) <- rep(1, p)
  fit <- snpClust(ld.ceph, h = 100)
  # expect_equal(fit$merge, prevfit$merge) 
  expect_equal(cumsum(fit$height), prevfit$height, tolerance = 1e-5)
})
