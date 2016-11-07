library(ratesCI)
context("Symmetry")

#n1 <- 150
#n2 <- 50
n1 <- 80
n2 <- 40
xs<-expand.grid(0:n1,0:n2)
x1 <- xs[,1]
x2 <- xs[,2]

#xs[2317,]
#scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="RD",distrib="poi",precis=10)$estimates
#scoreCI(x1=x2,n1=n2,x2=x1,n2=n1,contrast="RD",distrib="poi",precis=10,plot=T)
#scoretheta(-1.000,x1=28,x2=48,n1=40,n2=80,contrast="RD",distrib="poi",skew=T)
#scoretheta(1.000,x1=48,x2=28,n1=80,n2=40,contrast="RD",distrib="poi",skew=T)

#mle <- scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="RD",distrib="poi,precis=10)$estimates[,2]
#sum(mle == 1)
lcl <- scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="RD",distrib="poi",precis=10)$estimates[,1]
sum(lcl == 1)
#matrix(lcl,nrow=51,byrow=T)[1:10,1:10]
ucl <- scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="RD",distrib="poi",precis=10)$estimates[,3]
sum(ucl == 1)
#1/matrix(ucl,nrow=51,byrow=T)[1:10,1:10]

test_that("Transposed inputs produce inverted intervals", {
  expect_equal(
    unname(scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="RD",precis=10)$estimates[,c(1,3)]), 
    -unname(scoreCI(x1=x2,n1=n2,x2=x1,n2=n1,contrast="RD",precis=10)$estimates[,c(3,1)])
  )
  expect_equal(
    unname(scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="RD",distrib="poi",precis=10)$estimates[,c(1,3)]),
    -unname(scoreCI(x1=x2,n1=n2,x2=x1,n2=n1,contrast="RD",distrib="poi",precis=10)$estimates[,c(3,1)])
  )
  expect_equal(
    unname(scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="RR",precis=10)$estimates)[,c(1,3)],
    1/unname(scoreCI(x1=x2,n1=n2,x2=x1,n2=n1,contrast="RR",precis=10)$estimates)[,c(3,1)]
  )
  expect_equal(
    unname(scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="RR",distrib="poi",precis=10)$estimates)[,c(1,3)],
    1/unname(scoreCI(x1=x2,n1=n2,x2=x1,n2=n1,contrast="RR",distrib="poi",precis=10)$estimates)[,c(3,1)]
  )
  expect_equal(
    signif(unname(scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="OR",precis=1000)$estimates)[,c(1,3)],digits=4), 
    signif(1/unname(scoreCI(x1=x2,n1=n2,x2=x1,n2=n1,contrast="OR",precis=1000)$estimates)[,c(3,1)],digits=4)
  )
})

test_that("p-values consistent with confidence interval", {
  expect_equal(
    scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="RD",precis=10)$estimates[,"Lower"] > 0 | 
      scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="RD",precis=10)$estimates[,"Upper"] < 0,
    scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="RD",precis=10)$pval[,"pval2sided"]<0.05 
  )
  expect_equal(
    scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="RD",precis=10)$estimates[,"Lower"] > 0, 
    scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="RD",precis=10)$pval[,"pval.right"] < 0.025 
  )
  expect_equal(
    scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="RD",precis=10)$estimates[,"Upper"] < 0, 
    scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="RD",precis=10)$pval[,"pval.left"] < 0.025 
  )
  expect_equal(
    scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="RR",distrib="poi",precis=10)$estimates[,"Lower"] > 1 | 
      scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="RR",distrib="poi",precis=10)$estimates[,"Upper"] < 1, 
    scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="RR",distrib="poi",precis=10)$pval[,"pval2sided"]<0.05 
  )
  expect_equal(
    scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="RR",distrib="poi",precis=10,delta=2)$estimates[,"Lower"] > 2, 
    scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="RR",distrib="poi",precis=10,delta=2)$pval[,"pval.right"] < 0.025 
  )
  expect_equal(
    scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="RR",distrib="poi",precis=10,delta=2)$estimates[,"Upper"] < 2, 
    scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="RR",distrib="poi",precis=10,delta=2)$pval[,"pval.left"] < 0.025 
  )
  expect_equal(
    scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="OR",precis=10)$estimates[,"Lower"] > 1 | 
      scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="OR",precis=10)$estimates[,"Upper"] < 1, 
    scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="OR",precis=10)$pval[,"pval2sided"]<0.05 
  )
  expect_equal(
    scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="OR",precis=10)$estimates[,"Lower"] > 1, 
    scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="OR",precis=10)$pval[,"pval.right"] < 0.025 
  )
  expect_equal(
    scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="OR",precis=10)$estimates[,"Upper"] < 1, 
    scoreCI(x1=x1,n1=n1,x2=x2,n2=n2,contrast="OR",precis=10)$pval[,"pval.left"] < 0.025 
  )
})

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
  
if (FALSE) {
  
scoreCI(x1=0,x2=0,n1=20,n2=10,contrast="RR",plot=T)
scoreCI(x1=100,x2=0,n1=200,n2=5,plot=T)
scoreCI(x1=0,x2=100,n1=5,n2=200,plot=T)


scoreCI(x1=0,n1=80000,contrast="p",plot=T)

x1 <- 0
n1 <- 800
x2 <- NULL
n2 <- NULL
contrast="p"
distrib="bin"
bcf=TRUE
skew=TRUE
myfun <- function(theta, ...) {
  #  scoretheta(theta = theta, x1 = x1, x2 = x2, n1 = n1, n2 = n2, bcf = bcf,
  #             contrast = contrast, distrib = distrib, skew = skew)$score 
  scoretheta(theta = theta, x1 = x1, n1 = n1, contrast = "p", skew = T)$score 
}
myfun(theta=0.00064)
myfun(theta=0.00065)

myseq <- seq(0, 0.001, length.out = 400)
dim1 <- 1
sc <- array(sapply(myseq, function(x) myfun(x)), dim = c(dim1, length(myseq)))
myfun(theta=0.00064)
myfun(theta=0.00065)
scoretheta(theta=0.00064,x1=0,n1=800,distrib="bin",contrast="p",skew=T)$score
scoretheta(theta=0.00065,x1=0,n1=800,distrib="bin",contrast="p",skew=T)$score

}
