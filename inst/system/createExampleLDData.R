library("snpStats")

data(ld.example, package="snpStats")
str(ceph.1mb)

ld.ceph <- ld(ceph.1mb, stats=c("D.prime", "R.squared"), depth=100)

R2.100 <- ld.ceph$R.squared
Dprime.100 <- ld.ceph$D.prime

save(R2.100, file="data/R2.100.rda")
save(Dprime.100, file="data/Dprime.100.rda")
