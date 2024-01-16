context("Case of NA values in LD estimates")

check_missing_ld <- function() {
  skip_if_not_installed("snpStats")
  
  #Sys.setenv("OMP_THREAD_LIMIT" = 2)
  data("ld.example", package = "snpStats")
  p <- ncol(ceph.1mb)
  h <- p - 1
  ld.ceph <- snpStats::ld(ceph.1mb, depth = h, stats = "R.squared")
  if (!any(is.na(ld.ceph)))  {
    skip("No NA value: nothing to test here!")
  }
  sf <- system.file("data/ld.example.RData", package="snpStats")
  expected <- "497fcd532b5c2bcb082a0dad7ca0d44d"
  if (!(tools::md5sum(sf) == expected)) {
    skip("Different version of data('ld.example', package = 'snpStats')")
  }
}

test_that("NA values in LD estimates gives a warning/error in 'snpClust'", {
  skip_if_not_installed("snpStats")
  check_missing_ld()
  
  #Sys.setenv("OMP_THREAD_LIMIT" = 2)
  data("ld.example", package = "snpStats")
  p <- ncol(ceph.1mb)
  h <- p - 1
    
  ld.ceph <- snpStats::ld(ceph.1mb, depth = h, stats = "R.squared")
  ## In this unrealistically small example with only 90 subjects, it turns out
  ## that one of the LD estimates is NA due to the lack of genetic diversity in
  ## the sample:
  expect_warning(snpClust(ceph.1mb, h = h, stats = "R.squared"))
  expect_error(snpClust(ld.ceph, h = h))
})

test_that("NA values in LD estimates gives a warning/error in 'snpClust' (second version)", {
  # when check_missing_ld() skips the previous test: it means that snpClust does not produce NA
  skip_if_not_installed("snpStats")

  #Sys.setenv("OMP_THREAD_LIMIT" = 2)
  data("ld.example", package = "snpStats")
  p <- ncol(ceph.1mb)
  h <- p - 1
  
  ld.ceph <- snpStats::ld(ceph.1mb, depth = h, stats = "R.squared")
  ld.ceph[9,8] <- ld.ceph[8,9] <- NA # manually add NA
  
  expect_error(snpClust(ld.ceph, h = h))

  mat <- matrix(1, nrow = 10, ncol = 2)
  mat <- as(mat, "SnpMatrix")
  expect_warning(snpClust(mat))
})

## One way to correct this is to drop one of the incriminated SNPs
test_that("Dropping a SNP yielding NA values in LD fixes the NA problem", {
  skip_if_not_installed("snpStats")
  check_missing_ld()
  
  #Sys.setenv("OMP_THREAD_LIMIT" = 2)
  geno <- ceph.1mb[, -316]  ## drop one SNP leading to one missing LD value
  p <- ncol(geno)
  h <- p - 1
  
  ld.ceph <- snpStats::ld(geno, depth = h, stats = "R.squared", 
                          symmetric = TRUE)
  expect_true(all(!is.na(ld.ceph)))
  
  ## avoid LD values slightly greater than 1
  expect_warning(fit <- snpClust(geno, h = h, stats = "R.squared"),
                 "Forcing the LD similarity to be smaller than or equal to 1")
  expect_warning(fit10 <- snpClust(ld.ceph, h = h, stats = "R.squared"),
                 "Forcing the LD similarity to be smaller than or equal to 1")
  expect_identical(fit, fit10)
})

## Another way to correct this is to modify the genotype of a single sample
test_that("Modifying one genotype also fixes the NA problem", {
  skip_if_not_installed("snpStats")
  check_missing_ld()
  
  #Sys.setenv("OMP_THREAD_LIMIT" = 2)
  data("ld.example", package = "snpStats")
  p <- ncol(ceph.1mb)
  h <- p - 1
  nSamples <- nrow(ceph.1mb)

  ceph.1mb[4,286]@.Data[1,1] <- as.raw(3) ## to avoid NaNs
  ld.ceph <- snpStats::ld(ceph.1mb, depth = h, stats = "R.squared",
                          symmetric = TRUE)
  expect_true(all(!is.na(ld.ceph)))
  
  ## avoid LD values slightly greater than 1
  # ld.ceph10 <- round(ld.ceph, digits = 10)  
  expect_warning(fit <- snpClust(ceph.1mb, h = h, stats = "R.squared"),
                 "Forcing the LD similarity to be smaller than or equal to 1")
  expect_warning(fit10 <- snpClust(ld.ceph, h = h, stats = "R.squared"),
                 "Forcing the LD similarity to be smaller than or equal to 1")
  expect_identical(fit, fit10)
})
