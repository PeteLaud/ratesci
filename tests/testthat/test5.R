library(ratesci)
context("Consistency")

# Exploring the extent of a problem with the non-iterative SCAS method when x=1 or 2
# (or n-1 or n-2)
# and level>0.99. Checking against iterative version shows it doesn't match under these
# conditions.
# Seems to be OK for level≤0.99

#n <- round(runif(1000, 15, 9000), 0)
n <- seq(15, 9000, length.out = 1000)
xs <- rep_len(0:10, 1000)
rounded <- 10
for (level in c(0.9, 0.95, 0.99, 0.999)) {
  test_that("noniterative scas matches iterative version", {
    expect_equal(
      round(scoreci(x1 = xs, n1 = n, contrast = "p", level = level, precis = rounded + 1)$estimates[, c(1:3)], rounded),
      round(scaspci(x = xs, n = n, level = level)[, c(1:3)], rounded) #Env test bug 9Nov2021 relates to this line
    )
    expect_equal(
      round(scoreci(x1 = xs, n1 = n, contrast = "p", cc = T, level = level, precis = rounded + 1)$estimates[, c(1:3)], rounded),
      round(scaspci(x = xs, n = n, cc = T, level = level)[, c(1:3)], rounded) #Or is it this one? Fixed by removing Rmpfr code
    )
    expect_equal(
      round(scoreci(x1 = xs, n1 = n, contrast = "p", level = level, distrib = "poi", precis = rounded + 1)$estimates[, c(1:3)], rounded),
      round(scaspci(x = xs, n = n, level = level, distrib = "poi")[, c(1:3)], rounded)
    )
    expect_equal(
      round(scoreci(x1 = xs, n1 = n, contrast = "p", cc = T, level = level, distrib = "poi", precis = rounded + 1)$estimates[, c(1:3)], rounded),
      round(scaspci(x = xs, n = n, cc = T, level = level, distrib = "poi")[, c(1:3)], rounded)
    )
  })
}
