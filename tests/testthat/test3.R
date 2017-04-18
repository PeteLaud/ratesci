library(ratesci)
context("Examples")

test_that("legacy & new methods match published examples", {
  #Newcombe example (d)
  #Miettinen-Nurminen
  expect_equal(
    unname(round(scoreci(x1=5,x2=0,n1=56,n2=29,skew=F)$estimates[,c(1,3)],4)),
    c(-0.0326,0.1933)
  )
  #Mee
  expect_equal(
    unname(round(scoreci(x1=5,x2=0,n1=56,n2=29,skew=F,bcf=F)$estimates[,c(1,3)],4)),
    c(-0.0313,0.1926)
  )
  #Newcombe/'Score'/Square&add
  expect_equal(
    unname(round(moverci(x1=5,x2=0,n1=56,n2=29,type="wilson")[,c(1,3)],4)),
    c(-0.0381,0.1926)
  )
  #MOVER-J
  expect_equal(
    unname(round(moverci(x1=5,x2=0,n1=56,n2=29,type="jeff")[,c(1,3)],3)),
    c(-0.010,0.177)
  )
  #Gart Nam for skewness corrected RD (unable to replicate other examples)
  expect_equal(
    unname(round(scoreci(x1=7,x2=2,n1=25,n2=25,bcf=F,skew=F)$estimates[,c(1,3)],3)),
    c(-0.014,0.412)
  )
  expect_equal(
    unname(round(scoreci(x1=7,x2=2,n1=25,n2=25,bcf=F,skew=T)$estimates[,c(1,3)],3)),
    c(-0.014,0.414)
  )
  #Gart Nam for skewness corrected RR
  expect_equal(
    unname(round(scoreci(x1=c(8,6),n1=c(15,10),x2=c(4,6),n2=c(15,20),contrast="RR",bcf=F,skew=F)$estimates[,c(1,3)],3)),
    matrix(c(0.815,5.336,0.844,4.594),byrow=T,nrow=2)
  )
  expect_equal(
    unname(round(scoreci(x1=c(8,6),n1=c(15,10),x2=c(4,6),n2=c(15,20),contrast="RR",bcf=F,skew=T)$estimates[,c(1,3)],3)),
    matrix(c(0.806,6.148,0.822,4.954),byrow=T,nrow=2)
  )
  #Li Tang Wong for Poisson RR (don't quite match!)
  expect_equal(
    unname(round(moverci(x1=15,x2=41,n1=19017,n2=28010,distrib="poi",contrast="RR",type="jeff")[,c(1,3)],4)),
    c(0.2933, 0.9571)
  )
  expect_equal(
    unname(round(moverci(x1=60,x2=30,n2=54308.8,n1=51477.5,level=0.9,distrib="poi",contrast="RR",type="jeff")[,c(1,3)],4)),
    c(1.4668, 3.0586)
  )
  #Zou Donner for RR - they pre-round the Wilson interval to 3dps
  expect_equal(
    unname(round(moverci(x1=6,n1=22,x2=112,n2=842,distrib="bin",contrast="RR",type="wilson")[,c(1,3)],4)),
    c(0.9735, 3.7295)
  )
  #Fagerland Newcombe for OR
  expect_equal(
    unname(round(moverci(x1=c(24,29,7),n1=c(73,55,18),x2=c(53,11,1),n2=c(65,11,18),contrast="OR",type="wilson")[,c(1,3)],3)),
    matrix(c(0.050,0.242,0,0.427,1.395,75.890),byrow=T,nrow=3)
  )
  expect_equal(
    unname(round(scoreci(x1=c(24,29,7),n1=c(73,55,18),x2=c(53,11,1),n2=c(65,11,18),contrast="OR",skew=F)$estimates[,c(1,3)],3)),
    matrix(c(0.050,0.245,0,0.417,1.416,76.428),byrow=T,nrow=3)
  )
  #Fagerland Newcombe for RR
  expect_equal(
    unname(round(moverci(x1=c(24,29,7),n1=c(73,55,18),x2=c(53,11,1),n2=c(65,11,18),contrast="RR",type="wilson")[,c(1,3)],3)),
    matrix(c(0.282,0.563,0.398,0.761,1.277,40.751),byrow=T,nrow=3)
  )
  expect_equal(
    unname(round(scoreci(x1=c(24,29,7),n1=c(73,55,18),x2=c(53,11,1),n2=c(65,11,18),contrast="RR",skew=F)$estimates[,c(1,3)],3)),
    matrix(c(0.280,0.560,0.397,0.742,1.296,42.544),byrow=T,nrow=3)
  )
  #Hartung & Knapp stratified example
  #SCAS
  expect_equal(
    unname(round(scoreci(x1 = c(15,12,29,42,14,44,14,29,10,17,38,19,21),
                  x2 = c(9,1,18,31,6,17,7,23,3,6,12,22,19),
                  n1 = c(16,16,34,56,22,54,17,58,14,26,44,29,38),
                  n2 = c(16,16,34,56,22,55,15,58,15,27,45,30,38),
                  stratified = TRUE)$estimates[,c(1,3)],4)),
    c(0.2441,0.3698)
  )
  #TDAS
  expect_equal(
    unname(round(scoreci(x1 = c(15,12,29,42,14,44,14,29,10,17,38,19,21),
                  x2 = c(9,1,18,31,6,17,7,23,3,6,12,22,19),
                  n1 = c(16,16,34,56,22,54,17,58,14,26,44,29,38),
                  n2 = c(16,16,34,56,22,55,15,58,15,27,45,30,38),
                  stratified = TRUE, tdas = TRUE)$estimates[,c(1,3)],4)),
    c(0.1928,0.4647)
  )
  #paired methods examples from Agresti & Min 2005
  expect_equal(
    unname(round(pairbinci(x=c(53,16,8,9),contrast="RD")$estimates[,c(1,3)],4)),
    c(-0.0191,0.2052)
  )
  expect_equal(
    unname(round(pairbinci(x=c(53,16,8,9),contrast="OR")$estimates[c(1,3)],4)),
    c(0.8718,4.8816)
  )
})
