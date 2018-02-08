library("adjclust")

context("Case of NA values in LD estimates")

check_snp <- function() {
  if (!require("snpStats")) {
    skip("'snpStats' package not available")
  }
}

check_missing_ld <- function() {
  data("ld.example", package = "snpStats")
  p <- ncol(ceph.1mb)
  h <- p-1
  ld.ceph <- snpStats::ld(ceph.1mb, depth = h, stats = "R.squared")
  if (!any(is.na(ld.ceph)))  {
      skip("No NA value: nothing to test here!")
  }
}

test_that("NA values in LD estimates gives a warning/error in 'snpClust'", {
  check_snp()
  
  check_missing_ld()
  
  data("ld.example", package = "snpStats")
  p <- ncol(ceph.1mb)
  h <- p-1
    
  ld.ceph <- snpStats::ld(ceph.1mb, depth = h, stats = "R.squared")
  ## In this unrealistically small example with only 90 subjects, it turns out
  ## that one of the LD estimates is NA due to the lack of genetic diversity in
  ## the sample:
  expect_warning(snpClust(ceph.1mb, h = h, stats = "R.squared"))
  expect_error(snpClust(ld.ceph, h = h))
})

## One way to correct this is to drop one of the incriminated SNPs
test_that("Dropping a SNP yielding NA values in LD fixes the NA problem", {
  check_snp()
  
  check_missing_ld()
  
  geno <- ceph.1mb[, -316]  ## drop one SNP leading to one missing LD value
  p <- ncol(geno)
  h <- p-1
  
  ld.ceph <- snpStats::ld(geno, depth = h, stats = "R.squared")
  expect_true(all(!is.na(ld.ceph)))
  
  ## avoid LD values slightly greater than 1
  # ld.ceph10 <- round(ld.ceph, digits = 10)  
  
  fit <- snpClust(geno, h = h, stats = "R.squared")
  fit10 <- snpClust(ld.ceph, h = h, stats = "R.squared")
  expect_identical(fit, fit10)
})

## Another way to correct this is to modify the genotype of a single sample
test_that("Modifying one genotype also fixes the NA problem", {
  check_snp()
  
  check_missing_ld()
  
  data("ld.example", package = "snpStats")
  p <- ncol(ceph.1mb)
  h <- p-1
  nSamples <- nrow(ceph.1mb)

  ceph.1mb[4,286]@.Data[1,1] <- as.raw(3) ## to avoid NaNs
  ld.ceph <- ld(ceph.1mb, depth = h, stats = "R.squared")
  expect_true(all(!is.na(ld.ceph)))
  
  ## avoid LD values slightly greater than 1
  # ld.ceph10 <- round(ld.ceph, digits = 10)  
  
  fit <- snpClust(ceph.1mb, h = h, stats = "R.squared")
  fit10 <- snpClust(ld.ceph, h = h, stats = "R.squared")
  expect_identical(fit, fit10)
})
