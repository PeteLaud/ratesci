library(ratesCI)
context("Examples")

test_that("legacy methods match published examples", {
  expect_equal(
    unname(round(scoreCI(x1=5,x2=0,n1=56,n2=29,skew=F)$estimates[,c(1,3)],4)),
    c(-0.0326,0.1933)
  )
  expect_equal(
    unname(round(scoreCI(x1=5,x2=0,n1=56,n2=29,skew=F,bcf=F)$estimates[,c(1,3)],4)),
    c(-0.0313,0.1926)
  )
  expect_equal(
    unname(round(scoreCI(x1=7,x2=2,n1=25,n2=25,bcf=F,skew=F)$estimates[,c(1,3)],3)),
    c(-0.014,0.412)
  )
  expect_equal(
    unname(round(scoreCI(x1=7,x2=2,n1=25,n2=25,bcf=F,skew=T)$estimates[,c(1,3)],3)),
    c(-0.014,0.414)
  )
})
