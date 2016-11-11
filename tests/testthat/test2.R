library(ratesci)
context("p-values")

#n1 <- 150
#n2 <- 50
n1 <- 40
n2 <- 20
xs<-expand.grid(0:n1,0:n2)
x1 <- xs[,1]
x2 <- xs[,2]


for (skew in c(TRUE, FALSE)) {
  
  test_that("p-values consistent with confidence interval", {
    expect_equal(
      scoreci(x1=x1,n1=n1,x2=x2,n2=n2,skew=skew,contrast="RD",precis=10)$estimates[,"Lower"] > 0 | 
        scoreci(x1=x1,n1=n1,x2=x2,n2=n2,skew=skew,contrast="RD",precis=10)$estimates[,"Upper"] < 0,
      scoreci(x1=x1,n1=n1,x2=x2,n2=n2,skew=skew,contrast="RD",precis=10)$pval[,"pval2sided"]<0.05 
    )
    expect_equal(
      scoreci(x1=x1,n1=n1,x2=x2,n2=n2,skew=skew,contrast="RD",precis=10)$estimates[,"Lower"] > 0, 
      scoreci(x1=x1,n1=n1,x2=x2,n2=n2,skew=skew,contrast="RD",precis=10)$pval[,"pval.right"] < 0.025 
    )
    expect_equal(
      scoreci(x1=x1,n1=n1,x2=x2,n2=n2,skew=skew,contrast="RD",precis=10)$estimates[,"Upper"] < 0, 
      scoreci(x1=x1,n1=n1,x2=x2,n2=n2,skew=skew,contrast="RD",precis=10)$pval[,"pval.left"] < 0.025 
    )
    expect_equal(
      scoreci(x1=x1,n1=n1,x2=x2,n2=n2,skew=skew,contrast="RR",distrib="poi",precis=10)$estimates[,"Lower"] > 1 | 
        scoreci(x1=x1,n1=n1,x2=x2,n2=n2,skew=skew,contrast="RR",distrib="poi",precis=10)$estimates[,"Upper"] < 1, 
      scoreci(x1=x1,n1=n1,x2=x2,n2=n2,skew=skew,contrast="RR",distrib="poi",precis=10)$pval[,"pval2sided"]<0.05 
    )
    expect_equal(
      scoreci(x1=x1,n1=n1,x2=x2,n2=n2,skew=skew,contrast="RR",distrib="poi",precis=10,delta=2)$estimates[,"Lower"] > 2, 
      scoreci(x1=x1,n1=n1,x2=x2,n2=n2,skew=skew,contrast="RR",distrib="poi",precis=10,delta=2)$pval[,"pval.right"] < 0.025 
    )
    expect_equal(
      scoreci(x1=x1,n1=n1,x2=x2,n2=n2,skew=skew,contrast="RR",distrib="poi",precis=10,delta=2)$estimates[,"Upper"] < 2, 
      scoreci(x1=x1,n1=n1,x2=x2,n2=n2,skew=skew,contrast="RR",distrib="poi",precis=10,delta=2)$pval[,"pval.left"] < 0.025 
    )
    expect_equal(
      scoreci(x1=x1,n1=n1,x2=x2,n2=n2,skew=skew,contrast="OR",precis=10)$estimates[,"Lower"] > 1 | 
        scoreci(x1=x1,n1=n1,x2=x2,n2=n2,skew=skew,contrast="OR",precis=10)$estimates[,"Upper"] < 1, 
      scoreci(x1=x1,n1=n1,x2=x2,n2=n2,skew=skew,contrast="OR",precis=10)$pval[,"pval2sided"]<0.05 
    )
    expect_equal(
      scoreci(x1=x1,n1=n1,x2=x2,n2=n2,skew=skew,contrast="OR",precis=10)$estimates[,"Lower"] > 1, 
      scoreci(x1=x1,n1=n1,x2=x2,n2=n2,skew=skew,contrast="OR",precis=10)$pval[,"pval.right"] < 0.025 
    )
    expect_equal(
      scoreci(x1=x1,n1=n1,x2=x2,n2=n2,skew=skew,contrast="OR",precis=10)$estimates[,"Upper"] < 1, 
      scoreci(x1=x1,n1=n1,x2=x2,n2=n2,skew=skew,contrast="OR",precis=10)$pval[,"pval.left"] < 0.025 
    )
  })
  
}