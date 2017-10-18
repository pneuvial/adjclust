library("adjclust")
library("snpStats")

context("Ascending compatibility of the adjclust algorithm")

test_that("snpClust gives results identical to those of adjclust 0.3.0", {
    
    ## Note: this test depends on external data (genotypes) and functions 
    ## (snpStats::ld) which may change over time
    skip_on_cran() 
    
    data("res_adjclust_0.3.0", package = "adjclust")
    prevfit <- res_adjclust_0.3.0
    
    data("ld.example", package = "snpStats")
    h <- 100
    
    ld.ceph <- ld(ceph.1mb, depth = h, stats = "R.squared")
    
    p <- ncol(ceph.1mb)
    nSamples <- nrow(ceph.1mb)
    h <- 100
    ceph.1mb[4,286]@.Data[1,1] <- as.raw(3) ## to avoid NaNs
    ld.ceph <- ld(ceph.1mb, depth = h, stats = "R.squared")
    ld.ceph <- round(ld.ceph, digits = 10)

    fit <- snpClust(ld.ceph, h = 100)
    expect_equal(fit$merge, prevfit$merge)
    expect_equal(fit$height, prevfit$height, tolerance=1e-5)
})
