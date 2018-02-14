library(ratesci)
context("consistency")

#Exploring the extent of a problem with the non-iterative SCAS method when x=1 or 2
#(or n-1 or n-2)
#and level>0.99. Checking against iterative version shows it doesn't match under these
#conditions.
#Seems to be OK for level≤0.99

n <- 1000
xs <- 1
level <- 0.99
test_that("noniterative scas matches iterative version", {
  expect_equal(
    round(scoreci(x1=xs, n1=n, contrast="p", level=level, precis=11, plot=T)$estimates[ ,c(1:3)], 8),
    round(scasci.nonit(x=xs, n=n, level=level)[ ,c(1:3)], 8)
  )
})
