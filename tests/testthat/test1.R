library(ratesCI)
context("Symmetry")

n1 <- 150
n2 <- 50
xs<-expand.grid(0:n1,0:n2)
x1 <- xs[,1]
x2 <- xs[,2]

lcl <- scoreCI(x1,x2,n1,n2,contrast="RR",precis=10)$estimates[,1]
matrix(lcl,nrow=51,byrow=T)[1:10,1:10]
ucl <- scoreCI(x1,x2,n1,n2,contrast="RR",precis=10)$estimates[,3]
1/matrix(ucl,nrow=51,byrow=T)[1:10,1:10]

test_that("Transposed inputs produce inverted intervals", {
  expect_equal(
    unname(scoreCI(x1,x2,n1,n2,contrast="RD",precis=10)$estimates[,c(1,3)]), 
    -unname(scoreCI(x2,x1,n2,n1,contrast="RD",precis=10)$estimates[,c(3,1)])
  )
  expect_equal(
    unname(scoreCI(x1,x2,n1,n2,contrast="RD",distrib="poi",precis=10)$estimates[,c(1,3)]),
    -unname(scoreCI(x2,x1,n2,n1,contrast="RD",distrib="poi",precis=10)$estimates[,c(3,1)])
  )
  expect_equal(
    unname(scoreCI(x1,x2,n1,n2,contrast="RR",precis=10)$estimates)[,c(1,3)],
    1/unname(scoreCI(x2,x1,n2,n1,contrast="RR",precis=10)$estimates)[,c(3,1)]
  )
  expect_equal(
    unname(scoreCI(x1,x2,n1,n2,contrast="RR",distrib="poi",precis=10)$estimates)[,c(1,3)],
    1/unname(scoreCI(x2,x1,n2,n1,contrast="RR",distrib="poi",precis=10)$estimates)[,c(3,1)]
  )
  expect_equal(
    signif(unname(scoreCI(x1,x2,n1,n2,contrast="OR",precis=1000)$estimates)[,c(1,3)],digits=4), 
    signif(1/unname(scoreCI(x2,x1,n2,n1,contrast="OR",precis=1000)$estimates)[,c(3,1)],digits=4)
  )
})

test_that("p-values consistent with confidence interval", {
  expect_equal(
    scoreCI(x1,x2,n1,n2,contrast="RD",precis=10)$estimates[,"Lower"] > 0 | 
      scoreCI(x1,x2,n1,n2,contrast="RD",precis=10)$estimates[,"Upper"] < 0, 
    scoreCI(x1,x2,n1,n2,contrast="RD",precis=10)$pval[,"pval2sided"]<0.05 
  )
  expect_equal(
    scoreCI(x1,x2,n1,n2,contrast="RD",precis=10)$estimates[,"Lower"] > 0, 
    scoreCI(x1,x2,n1,n2,contrast="RD",precis=10)$pval[,"pval.right"] < 0.025 
  )
  expect_equal(
    scoreCI(x1,x2,n1,n2,contrast="RD",precis=10)$estimates[,"Upper"] < 0, 
    scoreCI(x1,x2,n1,n2,contrast="RD",precis=10)$pval[,"pval.left"] < 0.025 
  )
  expect_equal(
    scoreCI(x1,x2,n1,n2,contrast="RR",distrib="poi",precis=10)$estimates[,"Lower"] > 1 | 
      scoreCI(x1,x2,n1,n2,contrast="RR",distrib="poi",precis=10)$estimates[,"Upper"] < 1, 
    scoreCI(x1,x2,n1,n2,contrast="RR",distrib="poi",precis=10)$pval[,"pval2sided"]<0.05 
  )
  expect_equal(
    scoreCI(x1,x2,n1,n2,contrast="RR",distrib="poi",precis=10,delta=2)$estimates[,"Lower"] > 2, 
    scoreCI(x1,x2,n1,n2,contrast="RR",distrib="poi",precis=10,delta=2)$pval[,"pval.right"] < 0.025 
  )
  expect_equal(
    scoreCI(x1,x2,n1,n2,contrast="RR",distrib="poi",precis=10,delta=2)$estimates[,"Upper"] < 2, 
    scoreCI(x1,x2,n1,n2,contrast="RR",distrib="poi",precis=10,delta=2)$pval[,"pval.left"] < 0.025 
  )
  expect_equal(
    scoreCI(x1,x2,n1,n2,contrast="OR",precis=10)$estimates[,"Lower"] > 1 | 
      scoreCI(x1,x2,n1,n2,contrast="OR",precis=10)$estimates[,"Upper"] < 1, 
    scoreCI(x1,x2,n1,n2,contrast="OR",precis=10)$pval[,"pval2sided"]<0.05 
  )
  expect_equal(
    scoreCI(x1,x2,n1,n2,contrast="OR",precis=10)$estimates[,"Lower"] > 1, 
    scoreCI(x1,x2,n1,n2,contrast="OR",precis=10)$pval[,"pval.right"] < 0.025 
  )
  expect_equal(
    scoreCI(x1,x2,n1,n2,contrast="OR",precis=10)$estimates[,"Upper"] < 1, 
    scoreCI(x1,x2,n1,n2,contrast="OR",precis=10)$pval[,"pval.left"] < 0.025 
  )
})

test_that("legacy methods match published examples", {
  expect_equal(
    unname(round(scoreCI(5,0,56,29,skew=F)$estimates[,c(1,3)],4)),
    c(-0.0326,0.1933)
  )
  expect_equal(
    unname(round(scoreCI(5,0,56,29,skew=F,bcf=F)$estimates[,c(1,3)],4)),
    c(-0.0313,0.1926)
  )
  expect_equal(
    unname(round(scoreCI(7,2,25,25,bcf=F,skew=F)$estimates[,c(1,3)],3)),
    c(-0.014,0.412)
  )
  expect_equal(
    unname(round(scoreCI(7,2,25,25,bcf=F,skew=T)$estimates[,c(1,3)],3)),
    c(-0.014,0.414)
  )
})
  

scoreCI(0,0,20,10,contrast="RR",plot=T)
scoreCI(100,0,200,5,plot=T)
scoreCI(0,100,5,200,plot=T)
