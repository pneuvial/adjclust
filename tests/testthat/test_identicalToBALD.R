library("adjclust")

context("Comparison between the results of the 'BALD' and 'adjclust' packages")

test_that("BALD::cWard and adjClustBand give idenctical results on CEPH data (R.squared and D.prime)", {
              library("adjclust")
              data("ld_ceph", package="adjclust")
              h <- 100
              p <- 603
              dat <- list("R2"=R2.100, "Dprime"=Dprime.100)

              elts <- c("traceW", "gains", "merge", "height", "seqdist", "order", "labels")

              for (stat in c("R2", "Dprime")) {
                  x <- dat[[stat]]
                  res <- adjClustBand(x, p, h)

                  filename <- sprintf("cWard_ceph_%s,h=%s.rds", stat, h)
                  pathname <- system.file("extdata", filename, package="adjclust")
                  resB <- readRDS(pathname)

                  for (ee in elts) {
                      expect_equal(resB[[ee]],
                                   res[[ee]])
                  }
              }
          })
