#rm(list=ls())

#validation code for ratesCI methods


#Example from Newcombe 1998, used in Table 2 and Table S1 of Laud 2016
x1<-5; x2 <- 0; n1 <- 56; n2 <- 29


fround <- function(x,digits=6) paste("(",paste(format(round(x,digits=digits),nsmall=digits),collapse=", "),")",sep="")
mytab <-rbind(SCAS=
                c(
                  RDbin=fround(scoreCI(x1=x1,x2=x2,n1=n1,n2=n2,skew=T,contrast="RD")$estimates[,c(1,3)],3),
                  RDpoi=fround(scoreCI(x1=x1,x2=x2,n1=n1,n2=n2,skew=T,contrast="RD",dist="poi")$estimates[,c(1,3)],3),
                  RRbin=fround(scoreCI(x1=x1,x2=x2,n1=n1,n2=n2,skew=T,contrast="RR")$estimates[,c(1,3)],3),
                  RRpoi=fround(scoreCI(x1=x1,x2=x2,n1=n1,n2=n2,skew=T,contrast="RR",dist="poi")$estimates[,c(1,3)],3),
                  OR=fround(scoreCI(x1=x1,x2=x2,n1=n1,n2=n2,skew=T,contrast="OR")$estimates[,c(1,3)],3)
                ), MN=
                c(
                  RDbin=fround(scoreCI(x1=x1,x2=x2,n1=n1,n2=n2,skew=F,ontrast="RD")$estimates[,c(1,3)],3),
                  RDpoi=fround(scoreCI(x1=x1,x2=x2,n1=n1,n2=n2,skew=F,contrast="RD",dist="poi")$estimates[,c(1,3)],3),
                  RRbin=fround(scoreCI(x1=x1,x2=x2,n1=n1,n2=n2,skew=F,contrast="RR")$estimates[,c(1,3)],3),
                  RRpoi=fround(scoreCI(x1=x1,x2=x2,n1=n1,n2=n2,skew=F,contrast="RR",dist="poi")$estimates[,c(1,3)],3),
                  OR=fround(scoreCI(x1=x1,x2=x2,n1=n1,n2=n2,skew=F,contrast="OR")$estimates[,c(1,3)],3)
                ), MOVERJ=
                c(
                  RDbin=fround(MOVERCI(x1,x2,n1,n2,dist="bin",contrast="RD"),3),
                  RDpoi=fround(MOVERCI(x1,x2,n1,n2,dist="poi",contrast="RD"),3),
                  RRbin=fround(MOVERCI(x1,x2,n1,n2,dist="bin",contrast="RR"),3),
                  RRpoi=fround(MOVERCI(x1,x2,n1,n2,dist="poi",contrast="RR"),3),
                  OR=fround(MOVERCI(x1,x2,n1,n2,dist="bin",contrast="OR"),3)
                ), SCAScomp=
                c(
                  RDbin=fround(scoreCI(x1=x1,x2=x2,n1=n1,n2=n2,skew=T,contrast="RD",cc=0.25)$estimates[,c(1,3)],3),
                  RDpoi=fround(scoreCI(x1=x1,x2=x2,n1=n1,n2=n2,skew=T,contrast="RD",dist="poi",cc=0.25)$estimates[,c(1,3)],3),
                  RRbin=fround(scoreCI(x1=x1,x2=x2,n1=n1,n2=n2,skew=T,contrast="RR",cc=0.25)$estimates[,c(1,3)],3),
                  RRpoi=fround(scoreCI(x1=x1,x2=x2,n1=n1,n2=n2,skew=T,contrast="RR",dist="poi",cc=0.25)$estimates[,c(1,3)],3),
                  OR=fround(scoreCI(x1=x1,x2=x2,n1=n1,n2=n2,skew=T,contrast="OR",cc=0.25)$estimates[,c(1,3)],3)
                ), MOVERcomp=
                c(
                  RDbin=fround(MOVERCI(x1,x2,n1,n2,dist="bin",contrast="RD",cc=0.25),3),
                  RDpoi=fround(MOVERCI(x1,x2,n1,n2,dist="poi",contrast="RD",cc=0.25),3),
                  RRbin=fround(MOVERCI(x1,x2,n1,n2,dist="bin",contrast="RR",cc=0.25),3),
                  RRpoi=fround(MOVERCI(x1,x2,n1,n2,dist="poi",contrast="RR",cc=0.25),3),
                  OR=fround(MOVERCI(x1,x2,n1,n2,dist="bin",contrast="OR",cc=0.25),3)
                )
)
mytab

MOVERCI(5,0,56,29,type="Wilson",contrast="RD") #Newcombe 1998
MOVERCI(5,0,56,29,type="exact",contrast="RD")
MOVERCI(5,0,56,29,type="Jeff",contrast="RD",cc=0.25)

scoreCI(x1,x2,n1,n2,delta=0.5)

#Hartung & Knapp stratified example: (Table 3 of Laud 2016)
x1 <- c(15,12,29,42,14,44,14,29,10,17,38,19,21)
x2 <- c(9,1,18,31,6,17,7,23,3,6,12,22,19)
n1 <- c(16,16,34,56,22,54,17,58,14,26,44,29,38)
n2 <- c(16,16,34,56,22,55,15,58,15,27,45,30,38)

scoreCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=F,skew=TRUE,cc=F) $estimates[c(2,12),]

fround <- function(x,digits=6) paste(format(round(x[1],digits=digits),nsmall=digits)," (",paste(format(round(x[2:3],digits=digits),nsmall=digits),collapse=", "),")",sep="")
mytab <-rbind(
  SCASmh=c(
    RDbin=fround(scoreCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="MH",skew=T,tdas=F)$estimates[,c(2,1,3)],3),
    RRbin=fround(scoreCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="MH",skew=T,tdas=F)$estimates[,c(2,1,3)],2),
    OR=fround(scoreCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="MH",skew=T,tdas=F)$estimates[,c(2,1,3)],2)
  ),
  SCASiv=c(
    RDbin=fround(scoreCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="IV",skew=T,tdas=F)$estimates[,c(2,1,3)],3),
    RRbin=fround(scoreCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="IV",skew=T,tdas=F)$estimates[,c(2,1,3)],2),
    OR=fround(scoreCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="IV",skew=T,tdas=F)$estimates[,c(2,1,3)],2)
  ),
  TDAS=c(
    RDbin=fround(scoreCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="IV",tdas=T)$estimates[,c(2,1,3)],3),
    RRbin=fround(scoreCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="IV",tdas=T)$estimates[,c(2,1,3)],2),
    OR=fround(scoreCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="IV",tdas=T)$estimates[,c(2,1,3)],2)
  )
)
mytab





if(FALSE) {
  
  MOVERCI(5,0,56,29,type="Wilson",contrast="RD") #Newcombe 1998
  MOVERCI(5,0,56,29,type="exact",contrast="RD")
  MOVERCI(5,0,56,29,type="Jeff",contrast="OR",cc=0)
  
  
  
  MOVER(0:10,0,10,10,contrast="RR",dist="poi")
  
  x1 <- 3; x2<-2; n1=n2=10
  
  MOVER(x1,x2,n1,n2,contrast="OR",type="Wilson")
  1/MOVER(n1-x1,n2-x2,n1,n2,contrast="OR",type="Wilson")
  1/MOVER(49,49,50,50,contrast="OR")
  
  1/MOVER(41,15,28010,19017,contrast="RR",dist="poi")
  MOVER(15,41,19017,28010,contrast="RR",dist="poi")
  MOVER(c(15,20),41,19017,28010,contrast="RD",dist="poi")
  
  MOVER(0,5,15,20,contrast="OR",dist="bin")
  1/MOVER(5,0,20,15,contrast="OR",dist="bin")
  1/MOVER(15,15,15,20,contrast="OR",dist="bin")
  MOVER(15,15,20,15,contrast="OR",dist="bin")
  
  c=26
  MOVER(c,1,50,50,contrast="OR",dist="bin")
  1/MOVER(50-c,50-1,50,50,contrast="OR",dist="bin")
  1/MOVER(15,15,15,20,contrast="OR",dist="bin")
  MOVER(15,15,20,15,contrast="OR",dist="bin")
  
  x1=c(0,1,15,15);x2=c(1,0,15,15);n1=c(15,20,15,20);n2=c(20,15,20,15)
  
  MOVER(c(1,15,14,1),c(15,1,15,0),15,15,contrast="OR",dist="bin")
  
  xy <- expand.grid(0:15,0:15)
  MOVER(xy[,1],xy[,2],15,15,contrast="RR",dist="poi",type="Wil")
  MOVER(0,5,15,15,contrast="RR",dist="poi")
  
  1/(41/28010 / (15/19017))
  
  #diffBinconf.NJ(c(0,1),c(1,0),c(10,10),c(10,10),adj=F)
  
  
  #diffBinconf.NJ(56,48,70,80)
  diffBinconf.NJ(5,0,56,29)	
  diffBinconf.NJ(5,0,56,29,dist="poi")	
  #diffBinconf.NJ(0,5,10,10,k=0)	
}



root="/Users/petelaud/Documents/"
newpath=paste(root,"Main/Courses_papers/skewscore/R/",sep="")
source(paste(newpath,"ratesCI iterative weights skewness on point estimate.R",sep=""))
#source(paste(newpath,"ratesCI iterative weights.R",sep=""))
source(paste(newpath,"MOVERCI.R",sep=""))
source(paste(newpath,'diffBinconf.all.ssc.R',sep=""))
source(paste(newpath,'diffBinconf.all.OR.R',sep=""))
source(paste(newpath,'diffBinconf.all.RR.R',sep=""))
source(paste(newpath,"diffBinconf.metabin.R",sep=""))
source(paste(newpath,'diffBinconf.GNID.R',sep=""))
source(paste(newpath,'diffBinconf.GNIDRR.R',sep=""))

ratesCI(0,0,15,15,skew=F,cc=0,plot=T)
ratesCI(5,1,50,10,cc=0,plot=T,skew=T,contrast="RR")
ratesCI(0,5,10,100,cc=0,plot=T,skew=T,contrast="RD")
ratesCI(0,5,50,500,cc=0,plot=T,skew=T,contrast="RD",dist="poi")
ratesCI(0,0,15,15,cc=0.5,plot=T,skew=F)
ratesCI(15,0,15,15,plot=T,skew=F,cc=0,contrast="RD")#,dist="poi")
ratesCI(15,0,15,15,plot=T,skew=T,cc=0,contrast="RD")#,dist="poi")
score(1,5,1,50,10,skew=T,contrast="RR")
score(0.9894,15,0,15,15,skew=T,contrast="RD")
score(0,15,0,15,15,skew=T,contrast="RD")


#install.packages("xtable")
library(xtable)
library(meta)
fround <- function(x,digits=6) paste("(",paste(format(round(x,digits=digits),nsmall=digits),collapse=", "),")",sep="")
alpha<-0.05

#Hartung & Knapp stratified example:
x1 <- c(15,12,29,42,14,44,14,29,10,17,38,19,21)
x2 <- c(9,1,18,31,6,17,7,23,3,6,12,22,19)
n1 <- c(16,16,34,56,22,54,17,58,14,26,44,29,38)
n2 <- c(16,16,34,56,22,55,15,58,15,27,45,30,38)

#Hartung & Knapp stratified example, modified for an unequal allocation example:
x1 <- round(1.5*c(15,12,29,42,14,44,14,29,10,17,38,19,21))
x2 <- round(0.5*c(9,1,18,31,6,17,7,23,3,6,12,22,19))
n1 <- round(1.5*c(16,16,34,56,22,54,17,58,14,26,44,29,38))
n2 <- round(0.5*c(16,16,34,56,22,55,15,58,15,27,45,30,38))


#Hartung & Knapp stratified example:
x1hk <- c(15,12,29,42,14,44,14,29,10,17,38,19,21)
x2hk <- c(9,1,18,31,6,17,7,23,3,6,12,22,19)
n1hk <- c(16,16,34,56,22,54,17,58,14,26,44,29,38)
n2hk <- c(16,16,34,56,22,55,15,58,15,27,45,30,38)

#Example from Newcombe
x1<-5; x2 <- 0; n1 <- 56; n2 <- 29

#Example from Gart1990
x1<-4; x2 <- 5; n1 <- 10; n2 <- 40

i <- 12
myxtab <- list()
for (i in c(2,12)) {
x1 <- x1hk[i]
x2 <- x2hk[i]
n1 <- n1hk[i]
n2 <- n2hk[i]

mytab <-rbind(SCAS=
c(
RDbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0.5,skew=T,leftside=T,random=T,contrast="RD")$estimates[,c(1,3)],3),
RDpoi=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0.5,skew=T,leftside=T,random=T,contrast="RD",dist="poi")$estimates[,c(1,3)],3),
RRbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0.5,skew=T,leftside=T,random=T,contrast="RR")$estimates[,c(1,3)],3),
RRpoi=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0.5,skew=T,leftside=T,random=T,contrast="RR",dist="poi")$estimates[,c(1,3)],3),
OR=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0.5,skew=T,leftside=T,random=T,contrast="OR")$estimates[,c(1,3)],3)
), MN=
c(
RDbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0.5,skew=F,leftside=T,random=T,contrast="RD")$estimates[,c(1,3)],3),
RDpoi=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0.5,skew=F,leftside=T,random=T,contrast="RD",dist="poi")$estimates[,c(1,3)],3),
RRbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0.5,skew=F,leftside=T,random=T,contrast="RR")$estimates[,c(1,3)],3),
RRpoi=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0.5,skew=F,leftside=T,random=T,contrast="RR",dist="poi")$estimates[,c(1,3)],3),
OR=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0.5,skew=F,leftside=T,random=T,contrast="OR")$estimates[,c(1,3)],3)
), MOVERJ=
c(
RDbin=fround(MOVER(x1,x2,n1,n2,dist="bin",contrast="RD"),3),
RDpoi=fround(MOVER(x1,x2,n1,n2,dist="poi",contrast="RD"),3),
RRbin=fround(MOVER(x1,x2,n1,n2,dist="bin",contrast="RR"),3),
RRpoi=fround(MOVER(x1,x2,n1,n2,dist="poi",contrast="RR"),3),
OR=fround(MOVER(x1,x2,n1,n2,dist="bin",contrast="OR"),3)
), AN=
#RRbin=fround(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="Inverse",level=1-alpha,hakn=T,contrast="RR")$estimates[,c(1,3)],3), 
#RRpoi=fround(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="Inverse",level=1-alpha,hakn=T,contrast="RR",dist="poi")$estimates[,c(1,3)],3),
c(
RDbin=fround(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="Inverse",level=1-alpha,hakn=T,contrast="RD",incr=0)$estimates[,c(1,3)],3), 
RDpoi=fround(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="Inverse",level=1-alpha,hakn=T,contrast="RD",dist="poi")$estimates[,c(1,3)],3) ,
RRbin=fround(diffBinconf.all.RR(x1=x1,x2=x2,n1=n1,n2=n2)["Katz",1:2,],3), 
RRpoi=fround(crude.RRID(x1,x2,n1,n2,alpha=alpha),3),
OR=fround(exp(diffBinconf.all.OR(x1=x1,x2=x2,n1=n1,n2=n2)["Woolf",1:2,]),3) 
), SCAScomp=
c(
RDbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0.5,skew=T,leftside=T,random=T,contrast="RD",cc=0.25)$estimates[,c(1,3)],3),
RDpoi=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0.5,skew=T,leftside=T,random=T,contrast="RD",dist="poi",cc=0.25)$estimates[,c(1,3)],3),
RRbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0.5,skew=T,leftside=T,random=T,contrast="RR",cc=0.25)$estimates[,c(1,3)],3),
RRpoi=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0.5,skew=T,leftside=T,random=T,contrast="RR",dist="poi",cc=0.25)$estimates[,c(1,3)],3),
OR=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0.5,skew=T,leftside=T,random=T,contrast="OR",cc=0.25)$estimates[,c(1,3)],3)
), MOVERcomp=
c(
RDbin=fround(MOVER(x1,x2,n1,n2,dist="bin",contrast="RD",type="exact",k=0.25),3),
RDpoi=fround(MOVER(x1,x2,n1,n2,dist="poi",contrast="RD",type="exact",k=0.25),3),
RRbin=fround(MOVER(x1,x2,n1,n2,dist="bin",contrast="RR",type="exact",k=0.25),3),
RRpoi=fround(MOVER(x1,x2,n1,n2,dist="poi",contrast="RR",type="exact",k=0.25),3),
OR=fround(MOVER(x1,x2,n1,n2,dist="bin",contrast="OR",type="exact",k=0.25),3)
)
)
mytab
myxtab<-list(myxtab,xtable(mytab))
}

diffBinconf.OR(x1,x2,n1,n2)

x1 <- x1hk
x2 <- x2hk
n1 <- n1hk
n2 <- n2hk

mytab <-rbind(
I2=c(
#RDbin=round(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD")$homog[,c(3)],0),
#RRbin=round(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR")$homog[,c(3)],0),
#OR=round(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR")$homog[,c(3)],0)
RDbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="IV",skew=T,random=F,hk=F,fixtau=T)$Qtest["I2"],0),
RRbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="IV",skew=T,random=F,hk=F,fixtau=T)$Qtest["I2"],0),
OR=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="IV",skew=T,random=F,hk=F,fixtau=T)$Qtest["I2"],0)
),
#pval=c(
##RDbin=round(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD")$homog[,c(4)],3),
##RRbin=round(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR")$homog[,c(4)],3),
##OR=round(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR")$homog[,c(4)],3)
#RDbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="IV",skew=T,random=F,hk=F,fixtau=T)$Qtest["pval.het"],3),
#RRbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="IV",skew=T,random=F,hk=F,fixtau=T)$Qtest["pval.het"],3),
#OR=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="IV",skew=T,random=F,hk=F,fixtau=T)$Qtest["pval.het"],3)
#),
ANiv=c(
RDbin=fround(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",weight="Inverse")$estimates[,c(2,1,3)],3),
RRbin=fround((diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",weight="Inverse")$estimates[,c(2,1,3)]),2),
OR=fround((diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",weight="Inverse")$estimates[,c(2,1,3)]),2)
),
ANmh=c(
RDbin=fround(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",weight="MH")$estimates[,c(2,1,3)],3),
RRbin=fround(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",weight="MH")$estimates[,c(2,1,3)],2),
OR=fround(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",weight="MH")$estimates[,c(2,1,3)],2)
),
#MNmh=c(
#RDbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="MH",skew=F,random=F,hk=T,fixtau=T)$estimates[,c(2,1,3)],3),
#RRbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="MH",skew=F,random=F,hk=T,fixtau=T)$estimates[,c(2,1,3)],2),
#OR=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="MH",skew=F,random=F,hk=T,fixtau=T)$estimates[,c(2,1,3)],2)
#),
SCASmh=c(
RDbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="MH",skew=T,random=F,hk=T,fixtau=T)$estimates[,c(2,1,3)],3),
RRbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="MH",skew=T,random=F,hk=T,fixtau=T)$estimates[,c(2,1,3)],2),
OR=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="MH",skew=T,random=F,hk=F,fixtau=T)$estimates[,c(2,1,3)],2)
),
SCASiv=c(
RDbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="IV",skew=T,random=F,hk=T,fixtau=T)$estimates[,c(2,1,3)],3),
RRbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="IV",skew=T,random=F,hk=T,fixtau=T)$estimates[,c(2,1,3)],2),
OR=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="IV",skew=T,random=F,hk=F,fixtau=T)$estimates[,c(2,1,3)],2)
),
DL=c(
RDbin=fround(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="Inverse",contrast="RD")$estimates[,c(2,4,5)],3),
RRbin=fround(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="Inverse",contrast="RR")$estimates[,c(2,4,5)],2),
OR=fround(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="Inverse",contrast="OR")$estimates[,c(2,4,5)],2)
),
HKSJ=c(
RDbin=fround(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",hakn=T)$estimates[,c(2,4,5)],3),
RRbin=fround((diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",hakn=T)$estimates[,c(2,4,5)]),2),
OR=fround((diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",hakn=T)$estimates[,c(2,4,5)]),2)
),
TDAS=c(
RDbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="IV",random=T,hk=T,fixtau=T)$estimates[,c(2,1,3)],3),
RRbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",dist="bin",stratified=T,weights="IV",random=T,hk=T,fixtau=T,plotmax=12)$estimates[,c(2,1,3)],2),
OR=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="IV",random=T,hk=T,fixtau=T)$estimates[,c(2,1,3)],2)
)
)
mytab
xtable(mytab)



if(FALSE) {
x1=10
x2=0
n1=10
n2=10
alpha=0.05
stratified=T
random=F
bcf=F
delta=NULL
leftside=T
cc=F
skew=T
precis=5
dist="bin"
contrast="RR"

#Validation of ratesCI
x1=c(0,1); x2=c(0,0); n1=n2=15

x1=x2=0; n1=n2=5
x1=1; x2=0; n1=n2=15
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",dist="bin",skew=F,delta=-0.1,bcf=T,stratified=T) 
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",dist="bin",skew=F,bcf=T,stratified=T) 


x1<-c(322,14)
x2<-c(321,26)
n1<-c(379,31)
n2<-c(373,35)
system.time(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",dist="bin",skew=F,delta=-0.1,bcf=T,stratified=F,precis=5) )
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",dist="bin",skew=F,delta=-0.1,bcf=T,stratified=T,weights="IV") 
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",dist="bin",skew=T,delta=-0.1,bcf=T,stratified=T) 



x1 <- x2 <- 0
xs<-expand.grid(0:5,0:5)
x1 <- xs[,1]
x2 <- xs[,2]
n1 <- n2 <- 5
diffbinconf.score(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",dist="bin",skew=F,delta=0)  
diffbinconf.score(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",dist="bin",skew=T,delta=0)  

diffbinconf.score(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",dist="bin",skew=F,delta=1,bcf=T)  #pvals wrong
diffbinconf.score(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",dist="bin",skew=T,delta=1)  #??

RRconf.GN(1,0,5,5,skew=T)

result<-1
for(i in 1:length(x1)) {
c(result,diffbinconf.score(x1=x1[i],x2=x2[i],n1=n1,n2=n2,contrast="RR",dist="bin",skew=T,delta=0))  #?	
}
diffbinconf.score(c(0,1,5,0),c(0,1,5,5),5,5,skew=T,contrast="RR",bcf=T)
diffbinconf.score(5,0,5,5,skew=T,contrast="RR",bcf=T)
x1=1; x2=0; n1=n2=5

diffbinconf.score(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",dist="bin",skew=F,delta=0)
diffbinconf.score(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",dist="bin",skew=T,delta=0)




#binomial RD
#Laud&Dane
diffbinconf.score(x1=c(5,156,172,0,10),x2=c(0,160,179,0,0),n1=c(56,200,201,10,10),n2=c(29,200,201,10,20),dist="bin",contrast="RD",skew=F,delta=0)
diffbinconf.score(x1=c(5,156,172,0,10),x2=c(0,160,179,0,0),n1=c(56,200,201,10,10),n2=c(29,200,201,10,20),dist="bin",contrast="RD",skew=T)

cr<-crude.RD(x1=c(5,156,172,0,10),x2=c(0,160,179,0,0),n1=c(56,200,201,10,10),n2=c(29,200,201,10,20))
dimnames(cr)[[1]]

#Nurminen's dissertation, table II.2:
ratesCI(c(5,10),c(36,58),c(38,29),c(681,576),skew=T,stratified=T,weights="MH",plot=T,level=0.01)
ratesCI(c(5,10),c(36,58),c(38,29),c(681,576),skew=F,stratified=T,weights="MH",plot=T,level=0.01)
ratesCI(c(5,10),c(36,58),c(38,29),c(681,576),skew=T,stratified=F,plot=T,dist="poi")
ratesCI(c(5,10),c(36,58),c(38,29),c(681,576),skew=F,stratified=F,plot=T,dist="poi")
ratesCI(c(5,10),c(36,58),c(38,29),c(681,576),skew=T,stratified=F,plot=T,contrast="RR")
ratesCI(c(5,10),c(36,58),c(38,29),c(681,576),skew=F,stratified=F,plot=T,contrast="RR")

ratesCI(c(2,4),c(10,9),c(9,25),c(117,244),contrast="OR",skew=F,stratified=T,weights="MH",plot=T)




#Gart & Nam 1990 example 1-3
x1<-c(7,24,4)
x2<-c(2,1,5)
n1<-c(25,25,10)
n2<-c(25,25,40)
round(ratesCI(x1,x2,n1,n2,bcf=F,skew=F)$estimates,3)[,c(1,3)] #only the first one matches
round(ratesCI(x1,x2,n1,n2,bcf=F,skew=T),3)[,c(1,3)] #only the first one matches

#Stratified example in Gart & Nam Table 4
ratesCI(c(34,30),c(1,2),c(46,38),c(9,10),stratified=T,skew=F,bcf=F) #tick (0.374,0.748) with IV weights
ratesCI(c(34,30),c(1,2),c(46,38),c(9,10),stratified=T,skew=T,bcf=F) #quite different, shift to the right (0.382,0.760) G&N result shifts to left 
ratesCI(c(34,30),c(1,2),c(46,38),c(9,10),stratified=T,skew=T,bcf=F,weights="MH") #quite different, shift to the right (0.382,0.760) G&N result shifts to left 
ratesCI(64,3,84,19,skew=T,bcf=F) #note crude pooling shifts to the right (0.376,0.755) not sure that proves anything


#Fagerland (note his "MN" is actually Mee)
ratesCI(7,1,34,34,bcf=F,skew=F)

#Poisson RD


#Fagerland & Newcombe
ratesCI(24,53,73,65,contrast="OR")
ratesCI(29,11,55,11,contrast="OR",skew=T,bcf=T)
ratesCI(7,1,18,18,contrast="OR",cc=T)
ratesCI(24,53,73,65,contrast="RR")
ratesCI(29,11,55,11,contrast="RR")
ratesCI(7,1,18,18,contrast="RR")
ratesCI(7,1,18,18,contrast="RD")
ratesCI(7,1,18,18,contrast="RD",cc=T)


#Liu et al
ratesCI(11,25,309.5,295.3,dist="poi",skew=F) #tick
ratesCI(0,4,309.9,294.2,dist="poi",skew=F) #tick

#binomial RR
#Gart & Nam 1988 example 1
x1=8;x2=4;n1=15;n2=15
contrast="RR";dist="bin";skew=T;bcf=F
diffbinconf.score(x1=8,x2=4,n1=15,n2=15,contrast="RR",dist="bin",skew=F,bcf=F)  #tick (.815, 5.34)
diffbinconf.score(x1=8,x2=4,n1=15,n2=15,contrast="RR",dist="bin",skew=T,bcf=F)  #tick (.806, 6.15)
diffbinconf.score(x1=6,x2=6,n1=10,n2=20,contrast="RR",dist="bin",skew=F)  #tick
diffbinconf.score(x1=6,x2=6,n1=10,n2=20,contrast="RR",dist="bin",skew=T)  #tick

#example 2
x1=6;x2=6;n1=10;n2=20
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",dist="bin",skew=F,bcf=F,plot=T)  #tick (.844, 4.59)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="logRR",dist="bin",skew=F,bcf=F,plot=T)  #tick (.844, 4.59)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",dist="bin",skew=F,bcf=F,plot=T)  #tick (.844, 4.59)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",dist="bin",skew=T,bcf=F)  #tick (.822, 4.95)

#stratified example, table 6
x1 <- c(4,2,4,1); x2 <- c(5,3,10,3); n1 <- c(16,16,18,15); n2 <- c(79,87,90,82)
sum(x1)/sum(n1)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",dist="bin",skew=F,bcf=F,stratified=T,delta=1,weights="MH")  #hmm (1.29, 5.46) IV or (1.36, 5.08) MH
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",dist="bin",skew=T,bcf=F,stratified=T,delta=1,wt=(n1*x2/(n1+n2)))  # (1.22, 5.51) IV or (1.32, 5.12) MH


#Graham et al p.2076
diffbinconf.score(x1=41,x2=15,n1=28010,n2=19017,contrast="RR",dist="bin",skew=F,delta=1)
diffbinconf.score(x1=41,x2=15,n1=28010,n2=19017,contrast="RR",dist="bin",skew=T)  #doesnt quite match


#Fagerland
ratesCI(7,1,34,34,bcf=F,skew=F,contrast="RR",plot=T) #?
ratesCI(7,1,34,34,bcf=F,skew=T,contrast="RR")
ll <- ratesCI(7,1,34,34,bcf=F,skew=T,contrast="RD")$estimates[3]
(1/34+ll)/(1/34) 



#Fagerland & Newcombe
ratesCI(24,53,73,65,bcf=T,skew=F,contrast="RR")
ratesCI(29,11,55,11,bcf=T,skew=F,contrast="RR")



#Poisson RR


#validation against Graham et al p.2076
diffbinconf.score(x1=41,x2=15,n1=28010,n2=19017,contrast="RR",dist="poi",skew=F,delta=1)
ratesCI(x1=41,x2=15,n1=28010,n2=19017,contrast="RR",dist="poi",skew=F,delta=1)



#Binomial OR
ratesCI(x1=8,x2=4,n1=15,n2=15,contrast="OR",dist="bin",skew=T,bcf=F)   
ratesCI(x1=6,x2=6,n1=10,n2=20,contrast="OR",dist="bin",skew=T)  

#Confirming odd behaviour of LNCP for OR near p1=1
ratesCI(x1=49,x2=149,n1=50,n2=150,contrast="OR",skew=T,bcf=T)
pp1 <- 0.94; pp2 <- 0.5
truOR <- pp1*(1-pp2)/(pp2*(1-pp1))
xx1 <- rbinom(10000,50,pp1)
xx2 <- rbinom(10000,150,pp2)
cls <- ratesCI(x1=xx1,x2=xx2,n1=50,n2=150,contrast="OR",skew=T,bcf=T)$estimates
mean(cls[,1]>truOR)
hist(cls[,1])

ratesCI(x1=15,x2=1,n1=15,n2=15,contrast="OR",dist="bin",skew=F,bcf=F)  #why is UL=1??? - fixed by solving quadratic likelihood function when A=0. 
score(1,0,5,10,10,skew=F,random=F,stratified=F,contrast="OR",bcf=F)$zscore

ratesCI(19,0,20,20,contrast="OR")$estimates[1:3]
1/ratesCI(0,19,20,20,contrast="OR")$estimates[3:1]


diffbinconf.score(x1=41,x2=15,n1=28010,n2=19017,contrast="RD",dist="bin",skew=F,delta=0)
diffbinconf.score(x1=41,x2=15,n1=28010,n2=19017,contrast="RR",dist="bin",skew=F,delta=1)
diffbinconf.score(x1=41,x2=15,n1=28010,n2=19017,contrast="RD",dist="poi",skew=F,delta=0)
diffbinconf.score(x1=41,x2=15,n1=28010,n2=19017,contrast="RR",dist="poi",skew=F,delta=1)

diffbinconf.score(x1=41,x2=15,n1=28010,n2=19017,contrast="RD",dist="bin",skew=T)
diffbinconf.score(x1=41,x2=15,n1=28010,n2=19017,contrast="RR",dist="bin",skew=T)
diffbinconf.score(x1=41,x2=15,n1=28010,n2=19017,contrast="RD",dist="poi",skew=T)
diffbinconf.score(x1=41,x2=15,n1=28010,n2=19017,contrast="RR",dist="poi",skew=T)


#Checking for discontinuous score function
x1=c(2,1); x2=c(1,1); n1=n2=c(5,5)
x1=0; x2=5; n1=n2=10
#get tau estimate for fixtau version
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,skew=F,random=F,hk=F,bcf=T,contrast="OR",dist="bin",weights="IV")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,skew=F,random=T,hk=T,bcf=T,contrast="RR",dist="bin",weights="IV",fixq=F)
metabin(event.e=x1,n.e=n1,event.c=x2,n.c=n2,sm="RR",level=0.95,level.comb=0.95,method="Inverse",warn=T,addincr=F,incr=0.5,MH.exact=F,hakn=T,allstudies=F,backtransf=T)	
sc <- sapply(seq(0,15,0.01),function(x) score(x,x1,x2,n1,n2,skew=F,random=T,hk=T,stratified=T,contrast="RR",bcf=T,weights="IV",tau2=0,qhat=0.04258)$zscore) #-qnorm(alpha/2))
sc <- sapply(seq(0,2,0.01),function(x) score(x,x1,x2,n1,n2,skew=F,random=T,hk=T,stratified=T,contrast="RR",bcf=T,weights="IV",tau2=0.0)$rqhat) #-qnorm(alpha/2))
sc <- sapply(seq(0,15,0.01),function(x) score(x,x1,x2,n1,n2,skew=F,random=T,hk=T,stratified=T,contrast="RR",bcf=T,weights="IV")$zscore) #-qnorm(alpha/2))
sc <- sapply(seq(0,15,0.01),function(x) score(x,x1,x2,n1,n2,skew=F,random=T,hk=T,stratified=T,contrast="RR",bcf=T,tau2=0,weights="MH")$zscore) #-qnorm(alpha/2))
#sc <- sapply(seq(0,15,0.01),function(x) score(x,x1,x2,n1,n2,skew=F,random=T,hk=T,stratified=T,contrast="RR",bcf=T,tau2=NULL,weights="IV")$zscore) #-qnorm(alpha/2))
#sc <- sapply(1,function(x) score(x,x1,x2,n1,n2,skew=F,random=F,stratified=F,contrast="OR",bcf=F)$zscore) #-qnorm(alpha/2))
#scGN <- sapply(seq(0,15,0.1),function(x) score(x,x1,x2,n1,n2,skew=T,random=F,stratified=T,contrast="RR",bcf=F)$zscoreGN)#-qnorm(alpha/2))
plot(seq(0,15,0.01),sc,type="l",ylim=c(-20,20))
abline(h=c(-1,1)* qt(0.975,1))
abline(h=c(-1,1)* qnorm(0.975))
abline(v=c(0.89))

lines(seq(0,15,0.1),scGN)

score(theta=1.5,x1=x1,x2=x2,n1=n1,n2=n2,skew=T,random=F,stratified=T,contrast="RR",bcf=F)$zscore
score(1.5,x1,x2,n1,n2,skew=T,random=F,stratified=T,contrast="RR",bcf=F)$zscoreGN
theta=1.5
warnings()

plot(tan(seq(-0.99,0.99,0.01)*pi/2),sc)
abline(h=0)
abline(h=c(-1,1)*1.96)

z<-seq(-0.9,0.9,0.01)
plot(z,tan(z*pi/2))


sc <- sapply(seq(-1,1,0.01),function(x) score(x,x1,x2,n1,n2,skew=T,random=F,stratified=T,contrast="RD",bcf=F)$zscore-qnorm(alpha/2))
scGN <- sapply(seq(-1,1,0.01),function(x) score(x,x1,x2,n1,n2,skew=T,random=F,stratified=T,contrast="RD",bcf=F)$zscoreGN-qnorm(alpha/2))
plot(seq(-1,1,0.01),sc,type="l")
lines(seq(-1,1,0.01),scGN)


sc <- sapply(seq(0,2,0.01),function(x) score.RD(x,x1,x2,n1,n2,skew=F,random=F,stratified=T,contrast="RR"))
sc <- sapply(seq(0,2,0.01),function(x) score.RD(x,x1,x2,n1,n2,skew=F,random=F,stratified=T,contrast="RR"))
sc <- sapply(seq(0,2,0.01),function(x) score.RD(0.5,c(0,2),c(0,1),c(58,58),c(58,58),skew=F,random=T,stratified=T,contrast="RR"))


diffBinconf.BL(1,1,15,15)
diffbinconf.score(11,9,10,10)
diffBinconf.BL(c(1,2),c(1,2),c(15,15),c(15,15))
diffbinconf.score(c(0,0),c(0,0),c(10,10),c(10,10),delta=-0.2,stratified=F,skew=F)
diffbinconf.score(c(0,0),c(0,0),c(10,10),c(10,10),delta=-0.2,stratified=F,skew=F,dist="poi")
diffbinconf.score(c(0,0),c(0,0),c(10,10),c(10,10),delta=-0.2,stratified=T)

diffbinconf.score(c(0,0),c(0,0),c(10,10),c(10,10),stratified=T,skew=F)
diffbinconf.score(0,0,20,20,skew=F)

diffbinconf.score(c(0,0),c(0,0),c(10,10),c(10,10),stratified=T,skew=T)
diffbinconf.score(0,0,20,20,skew=T)

#Whitehead & Whitehead stratified example:
if(FALSE){
x1 <- c(0,59,0,13,60,5,1,25,43,1,2,0,10,18,32,20) #[-c(1,12)]
x2 <- c(0,88,5,22,109,20,6,36,52,3,1,0,21,34,48,39) #[-c(1,12)]
n1 <- c(508,3903,406,1721,8700,186,193,1048,233,68,45,58,49,534,416,419) #[-c(1,12)]
n2 <- c(504,3922,379,1706,8654,194,196,1004,219,63,42,58,48,529,424,465) #[-c(1,12)]

#excluding the zeros for direct comparison with W&W
x1 <- c(59,0,13,60,5,1,25,43,1,2,10,18,32,20)
x2 <- c(88,5,22,109,20,6,36,52,3,1,21,34,48,39)
n1 <- c(3903,406,1721,8700,186,193,1048,233,68,45,49,534,416,419)
n2 <- c(3922,379,1706,8654,194,196,1004,219,63,42,48,529,424,465)

#unequal stratum sizes
x1 <- c(0,59,0,13,60,5,1,25,43,1,2,0,10,18,32,20)
x2 <- c(0,88,5,22,109,20,6,36,52,3,1,0,21,34,48,39)
n1 <- c(508,3903,406,1721,8700,186,193,1048,233,68,45,58,49,534,416,419)
n2 <- c(504,3922,379,1706,8654,194,196,1004,219,63,42,58,48,529,424,465)



metabin(event.e=x1,n.e=n1,event.c=x2,n.c=n2,sm="RD",level=0.95,level.comb=0.95,method="Inverse",warn=T,addincr=F,MH.exact=F,hakn=T,allstudies=F)	#standard IV metaanalysis
metabin(event.e=x1,n.e=n1,event.c=x2,n.c=n2,sm="RD",level=0.95,level.comb=0.95,method="Inverse",warn=T,addincr=T,incr=(1/16),MH.exact=F,hakn=T,allstudies=F)	#HK version
metabin(event.e=x1,n.e=n1,event.c=x2,n.c=n2,sm="RR",level=0.95,level.comb=0.95,method="Inverse",warn=T,addincr=T,incr=0.5,MH.exact=F,hakn=F,allstudies=F,backtransf=F)	
metabin(event.e=x1,n.e=n1,event.c=x2,n.c=n2,sm="OR",level=0.95,level.comb=0.95,method="Inverse",warn=T,addincr=T,incr=0.5,MH.exact=F,hakn=T,allstudies=F,backtransf=F)	
(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="Inverse",contrast="RD",hakn=F,addincr=F)$estimates )  #addincr=T & allstudies=F for consistency with HK results
(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="Inverse",contrast="RD",hakn=T,addincr=F)$estimates )  #addincr=T & allstudies=F for consistency with HK results
(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="Inverse",contrast="RR",hakn=T)$estimates )  #addincr=T & allstudies=F for consistency with HK results
(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="Inverse",contrast="OR",hakn=F)$estimates )

ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=T,random=F,hk=F,bcf=T,contrast="RD",dist="bin",weights="IV")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=F,random=T,hk=T,bcf=T,contrast="RD",dist="bin",weights="IV")
log(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=T,random=F,bcf=T,contrast="RR",dist="bin",weights="MH")$estimates)
log(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=F,random=T,hk=T,bcf=T,contrast="RR",dist="bin",weights="IV")$estimates)
log(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=F,random=T,hk=T,bcf=T,contrast="OR",dist="bin",weights="IV")$estimates)

ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=T,contrast="RD",dist="bin",weights="IV")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=F,contrast="RD",dist="bin",weights="IV")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=T,contrast="RD",dist="bin",weights="MH")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=F,contrast="RD",dist="bin",weights="MH")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=T,contrast="RR",dist="bin",weights="MH")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=F,contrast="RR",dist="bin",weights="MH")



ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=F,delta=0,skew=T,contrast="RD",dist="bin",weights="MH")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=F,delta=0,skew=F,contrast="RD",dist="bin",weights="MH")


x1 <- c(1,1)
x2 <- c(0,0)
n1 <- c(508,3903)
n2 <- c(504,3922)

#Mehrotra & Railkar example 1
x1 <- c(26,25,11)
x2 <- c(21,11,8)
n1 <- c(453,174,70)
n2 <- c(464,165,69)

#Mehrotra & Railkar example 2
x1 <- c(48,48)
x2 <- c(53,57)
n1 <- c(100,60)
n2 <- c(100,60)
contr="OR"; dist="bin"
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,random=T,hk=T,bcf=T,contrast=contr,dist=dist,weights="IV",fixtau=T,plot=T,fixq=T)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,random=T,hk=T,bcf=T,contrast="OR",fixq=T,dist=dist,weights="IV",fixtau=T,plot=T)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,random=T,hk=T,bcf=T,contrast="RD",fixq=F,dist=dist,weights="IV",fixtau=T,plot=T)
diffBinconf.metabin(x1=x1,n1=n1,x2=x2,n2=n2,contrast=contr,dist=dist,level=0.95,weight="Inverse",hakn=T)	

#Chan & Wang Poisson example
x1 <- c(3,6,9)
x2 <- c(9,24,33)
n1 <- c(31323,26881,58203)
n2 <- c(30953,26783,57736)

#Tarone/Gart - NS heterog with CMH weights but sig with IV (93% weight given to stratum 1)
x1 <- c(4,38,119,221,259,310,226,65)[1:8]
x2 <- c(1,16,30,71,102,130,133,40)[1:8]
n1 <- c(181343,146207,121374,111353,83004,55932,29007,7538)[1:8]
n2 <- c(172675,123065,96216,92051,72159,54722,32195,8328)[1:8]

(x1/n1)/(x2/n2)
diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="Inverse",contrast="RR",hakn=TRUE) 
metabin(event.e=x1,n.e=n1,event.c=x2,n.c=n2,method="Inverse",sm="RR",hakn=FALSE)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="IV",random=F,hk=F,fixtau=T,fixvs=T,plot=T,plotmax=12,skew=T)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="IV",random=T,hk=T,fixtau=T,fixvs=T,plot=T,plotmax=12,skew=F)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="IV",random=T,hk=T,fixtau=F,fixvs=F,plot=T,plotmax=12,skew=F)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="IV",random=T,hk=T,fixtau=T,fixvs=T,plot=T,plotmax=12,skew=F,weightiter=20) #iterative weights makes little difference here

ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="MH",random=F,hk=F,fixtau=F,fixvs=F,plot=T,plotmax=12,skew=T)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="MH",random=T,hk=F,fixtau=F,fixvs=F,plot=T,plotmax=12,skew=T)
scr <-
 score(theta=1,x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,weights="IV",contrast="RR")$Stheta
sco <- score(theta=2.24,x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,weights="MH",contrast="OR")$Stheta




qqnorm(scr)
qqnorm(sco)

#Asia CAP/FOCUS meta-analysis
x1 <- c(244,139,151)
x2 <- c(200,121,112)
n1 <- c(301,154,189)
n2 <- c(297,157,156)
(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="Inverse",contrast="RD",hakn="T")$estimates )
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=F,leftside=T,random=F,bcf=T,contrast="RD",dist="bin",weights="MH")


x2=c(rep(0,12),0)
dev.off()

#CRASH meta-analysis - http://www.bmj.com/content/314/7098/1855.full
x1 <- c(9,16,16,26,35,114,8,44,34,33,4,19,38)
x2 <- c(13,22,16,13,36,38,9,47,7,21,4,21,49)
n1 <- c(17,55,67,49,81,201,50,81,72,68,12,133,175)
n2 <- c(18,55,28,27,83,74,50,80,16,62,12,136,195)
round((x1/n1+x2/n2)/2,2)



#CRASH meta-analysis updated with final results- Lancet
x1 <- c(9,16,16,26,35,114,8,44,34,33,4,19,38,13,1052)
x2 <- c(13,22,16,13,36,38,9,47,7,21,4,21,49,5,893)
n1 <- c(17,55,67,49,81,201,50,81,72,68,12,133,175,98,4985)
n2 <- c(18,55,28,27,83,74,50,80,16,62,12,136,195,54,4979)

contr="RR"; dist="bin"
mymeta <- metabin(event.e=x1,n.e=n1,event.c=x2,n.c=n2,sm=contr,backtransf=T,level=0.95,level.comb=0.95,method="Inverse",warn=T,addincr=F,incr=0.5,hakn=T,allstudies=F)	
funnel(mymeta)
diffBinconf.metabin(x1=x1,n1=n1,x2=x2,n2=n2,contrast=contr,dist=dist,level=0.95,weight="Inverse",hakn=F)	
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=F,delta=1,skew=F,leftside=T,random=F,hk=F,bcf=T,contrast=contr,dist=dist,weights="IV")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,leftside=T,random=F,hk=F,bcf=T,contrast=contr,dist=dist,weights="IV")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=T,leftside=T,random=F,hk=F,bcf=T,contrast=contr,dist=dist,weights="IV")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=T,leftside=T,random=F,hk=F,bcf=T,contrast=contr,dist=dist,weights="MH")
metabin(event.e=x1,n.e=n1,event.c=x2,n.c=n2,sm=contr,backtransf=T,level=0.95,level.comb=0.95,method="Inverse",warn=T,addincr=F,incr=0.5,hakn=T,allstudies=F)	
diffBinconf.metabin(x1=x1,n1=n1,x2=x2,n2=n2,contrast=contr,dist=dist,level=0.95,weight="Inverse",hakn=T)	
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,skew=F,weights="IV",random=T,hk=T,fixtau=T,contrast="OR",dist="bin")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,random=T,hk=T,bcf=T,contrast=contr,dist=dist,weights="IV",fixtau=T,plot=T)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,random=T,hk=T,bcf=T,contrast=contr,dist=dist,weights="IV",fixtau=T,plot=T,robust=T)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,random=T,hk=T,bcf=T,contrast=contr,dist=dist,weights="MH",plot=T)

#excluding Faupel (and Dearden?)
metabin(event.e=x1[-3],n.e=n1[-3],event.c=x2[-3],n.c=n2[-3],sm=contr,backtransf=T,level=0.95,level.comb=0.95,method="MH",warn=T,addincr=F,incr=0.5,hakn=F,allstudies=F)	
metabin(event.e=x1[-3],n.e=n1[-3],event.c=x2[-3],n.c=n2[-3],sm=contr,backtransf=T,level=0.95,level.comb=0.95,method="Inverse",warn=T,addincr=F,incr=0.5,hakn=T,allstudies=F)	
ratesCI(x1=x1[-3],x2=x2[-3],n1=n1[-3],n2=n2[-3],stratified=T,delta=1,skew=F,leftside=T,random=F,bcf=T,contrast=contr,dist="bin",weights="IV",plot=T)
ratesCI(x1=x1[-c(3,10)],x2=x2[-c(3,10)],n1=n1[-c(3,10)],n2=n2[-c(3,10)],stratified=T,delta=1,skew=F,leftside=T,random=F,bcf=T,contrast=contr,dist="bin",weights="IV",plot=T)
ratesCI(x1=x1[-3],x2=x2[-3],n1=n1[-3],n2=n2[-3],stratified=T,delta=1,skew=F,leftside=T,random=T,hk=T,bcf=T,contrast=contr,dist="bin",weights="IV",fixtau=T,plot=T)

#final CRASH results - death within 2 weeks
x1 <- 1052; n1 <- 4985; x2 <- 893; n2 <- 4979
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=F,delta=1,skew=T,bcf=T,contrast="RR",dist="bin")
metabin(event.e=x1,n.e=n1,event.c=x2,n.c=n2,sm="RR",backtransf=T,level=0.95,level.comb=0.95,method="Inverse",warn=T,addincr=F,incr=0.5,hakn=T,allstudies=F)	

#AZ reporting standards doc / Kevin Carroll position paper
#eg1
x1 <- c(0,6,2,7,0,0,0,0,3,1,1,20)
n1 <- c(59.5,100.9,31.1,118.1,92.8,53,47.8,121.3,128.9,90.9,52.4,948.6)
x2 <- c(1,2,1,0,0,1,0,0,3,0,0,3)
n2 <- c(35.3,36.1,25.7,48.6,73.4,18.1,5,34.7,120.9,19.5,20.9,37.2)

#eg2
x1 <- c(6,6,8,4,2,2,127,8,5,2)
n1 <- c(127.5,94.8,295.2,150.3,28.5,67.6,16180.1,2185,875,528)
x2 <- c(1,1,7,7,0,2,104,8,1,4)
n2 <- c(61.7,41.1,298,131.5,27.4,61.1,15891.2,1269,809,489)

#disconcerting infinite UL with MH or IV weights for RR with TDAS - maybe just a "few-stratum" effect
x1 <- c(0,11,2)
n1 <- c(10,100,50)
x2 <- c(0,36,3)
n2 <- c(10,100,50)

x1 <- c(23,2)
n1 <- c(100,50)
x2 <- c(36,3)
n2 <- c(100,50)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="IV",random=T,hk=T,fixtau=F,weightiter=2,plot=T,plotmax=12)
#ratesCI(x1=x2,x2=x1,n1=n2,n2=n1,contrast="RR",stratified=T,weights="IV",random=T,hk=T,fixtau=F,weightiter=2,plot=T,plotmax=12)
metabin(event.e=x1,n.e=n1,event.c=x2,n.c=n2,sm="OR",level=0.95,level.comb=0.95,method="Inverse",warn=T,addincr=F,incr=0.5,MH.exact=F,hakn=T,allstudies=F,backtransf=T)	

#try with more strata
ratesCI(x1=rep(x1,3),x2=rep(x2,3),n1=rep(n1,3),n2=rep(n2,3),contrast="RR",stratified=T,weights="IV",random=T,hk=T,fixtau=F,weightiter=2,plot=T,plotmax=12)
metabin(event.e=rep(x1,3),n.e=rep(n1,3),event.c=rep(x2,3),n.c=rep(n2,3),sm="RR",level=0.95,level.comb=0.95,method="Inverse",warn=T,addincr=F,incr=0.5,MH.exact=F,hakn=T,allstudies=F,backtransf=T)	
#pooled
ratesCI(x1=sum(x1),x2=sum(x2),n1=sum(n1),n2=sum(n2),contrast="RD",stratified=T,weights="IV",random=T,hk=T,fixtau=F,weightiter=2,plot=T,plotmax=12)

#what about RD? - small Q results in narrower CI than fixed effects, as per HKSJ
x1 <- c(0,10,2)
n1 <- c(10,100,50)
x2 <- c(0,11,3)
n2 <- c(10,100,50)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="IV",random=T,hk=T,fixtau=F,weightiter=2,plot=T,plotmax=12)
#ratesCI(x1=x2,x2=x1,n1=n2,n2=n1,contrast="RR",stratified=T,weights="IV",random=T,hk=T,fixtau=F,weightiter=2,plot=T,plotmax=12)
metabin(event.e=x1,n.e=n1,event.c=x2,n.c=n2,sm="RD",level=0.95,level.comb=0.95,method="Inverse",warn=T,addincr=F,incr=0.5,MH.exact=F,hakn=T,allstudies=F,backtransf=T)	



x1 <- c(30,15)
n1 <- c(100,50)
x2 <- c(35,14)
n2 <- c(100,50)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="IV",random=T,hk=T,fixtau=F,weightiter=2,plot=T,plotmax=2)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="IV",random=F,hk=T,fixtau=F,weightiter=2,plot=T,plotmax=12)
ratesCI(x1=x2,x2=x1,n1=n2,n2=n1,contrast="RR",stratified=T,weights="IV",random=T,hk=T,fixtau=F,weightiter=2,plot=T,plotmax=12)

ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,random=T,hk=T,bcf=T,contrast="RD",dist="bin",weights="IV",fixq=F,fixtau=F,plot=T) #random
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,random=T,hk=T,bcf=T,contrast="RR",dist="bin",weights="IV",fixvs=T,robust=F,plot=T) #random
ratesCI(x1=x2,x2=x1,n1=n2,n2=n1,stratified=T,delta=1,skew=F,random=T,hk=T,bcf=T,contrast="RR",dist="bin",weights="IV",fixvs=T,robust=F,plot=T) #random
metabin(event.e=x1,n.e=n1,event.c=x2,n.c=n2,sm="RR",backtransf=T,level=0.95,level.comb=0.95,method="MH",warn=T,addincr=F,incr=0.5,hakn=T,allstudies=F)	

x1 <- c()
x2 <- c()
n1 <- c()
n2 <- c()

#Collins Pre-eclampsia
#Another example of the unreliability of IV weights with a stratum with small p
labels <- c("Wesley","Flowers","Menzies","Fallis","Cuadros","Landesman","Kraus","Tervila","Campbell")
x1 <- c(14,21,14,6,12,138,15,6,65)
x2 <- c(14,17,24,18,35,175,20,2,40)
n1 <- c(131,385,57,38,1011,1370,506,108,153)
n2 <- c(136,134,48,40,760,1336,524,103,102)


#Eom pneumonia
labels <- c("Cheadle","Driks","Laggner","Reusser","Eddleston","Apte","Martin","Metz","Pickworth","Ryan","Ben-Menachem","Cloud","Maier","Prod'hom","Mustafa","Thomason","Cook","Hanisch","Moesgaard","O'Keefe","Yildizdas","Kantorova","Misra")
x1 <- c(13,16,2,7,10,13,0,12,5,7,13,7,14,21,9,27,114,10,11,14,20,7,2)
x2 <- c(3,7,0,8,3,9,4,15,6,8,6,4,10,10,3,30,98,12,19,10,17,5,5)
n1 <- c(100,69,16,19,30,16,56,84,44,56,100,81,51,80,16,80,596,57,77,49,42,71,45)
n2 <- c(100,61,16,21,30,18,61,79,39,58,100,38,47,83,15,80,604,57,87,47,42,75,47)


#CRASH meta-analysis updated with final results, excluding Faupel
x1 <- c(9,16,16,26,35,114,8,44,34,33,4,19,38,1,13,1052)[-c(3)]
x2 <- c(13,22,16,13,36,38,9,47,7,21,4,21,49,0,5,893)[-c(3)]
n1 <- c(17,55,67,49,81,201,50,81,72,68,12,133,175,5,98,4985)[-c(3)]
n2 <- c(18,55,28,27,83,74,50,80,16,62,12,136,195,5,54,4979)[-c(3)]

#CRASH meta-analysis updated with final results, excluding the 2 most extreme
x1 <- c(9,16,16,26,35,114,8,44,34,33,4,19,38,1,13,1052)[-c(3,14)]
x2 <- c(13,22,16,13,36,38,9,47,7,21,4,21,49,0,5,893)[-c(3,14)]
n1 <- c(17,55,67,49,81,201,50,81,72,68,12,133,175,5,98,4985)[-c(3,14)]
n2 <- c(18,55,28,27,83,74,50,80,16,62,12,136,195,5,54,4979)[-c(3,14)]
labels <- 1:length(x1)



#CRASH meta-analysis updated with final results, without CRASH itself- from Lancet paper
x1 <- c(9,16,16,26,35,114,8,44,34,33,1,4,13,19,38,0)
x2 <- c(13,22,16,13,36,38,9,47,7,21,0,4,5,21,49,0)
n1 <- c(17,55,67,49,81,201,50,81,72,68,5,12,98,133,175,30)
n2 <- c(18,55,28,27,83,74,50,80,16,62,5,12,54,136,195,30)
#CRASH meta-analysis updated with final results, plus some other trials- from Lancet paper
x1 <- c(9,16,16,26,35,114,8,44,34,33,1,4,13,19,38,0,1052)
x2 <- c(13,22,16,13,36,38,9,47,7,21,0,4,5,21,49,0,893)
n1 <- c(17,55,67,49,81,201,50,81,72,68,5,12,98,133,175,30,4985)
n2 <- c(18,55,28,27,83,74,50,80,16,62,5,12,54,136,195,30,4979)
#CRASH meta-analysis updated with final results, plus some other trials- minus Faupel & Chacon
x1 <- c(9,16,16,26,35,114,8,44,34,33,1,4,13,19,38,0,1052)[-c(3,11)]
x2 <- c(13,22,16,13,36,38,9,47,7,21,0,4,5,21,49,0,893)[-c(3,11)]
n1 <- c(17,55,67,49,81,201,50,81,72,68,5,12,98,133,175,30,4985)[-c(3,11)]
n2 <- c(18,55,28,27,83,74,50,80,16,62,5,12,54,136,195,30,4979)[-c(3,11)]
labels <- 1:length(x1)


#metabin smoking dataset
data(smoking)
attach(smoking)
attach(lungcancer)
x1 <- d.smokers
x2 <- d.nonsmokers
n1 <- py.smokers
n2 <- py.nonsmokers

#Helicopter ambulance Example from cochrane library???  Helicopter emergency medical services for adults with major trauma
x1 <- c(10,32,64,5,79,23,9,9,11,28,39,215,31,759,20,135,19,25,682,274,55,88,183,98,7813,128,407,546,32,64,28,759,19,130,10,9,215,88)[1:28]
x2 <- c(19,51,33,31,72,73,52,9,31,22,51,425,42,1845,39,1469,26,12,1874,540,72,283,308,26,17775,77,643,5765,51,33,22,1845,26,170,19,52,425,283)[1:28]
n1 <- c(150,104,347,42,288,565,77,32,15,92,122,2292,140,3017,104,3870,89,478,10049,2717,854,516,703,333,61909,446,2949,2090,104,347,92,3017,89,492,150,77,2292,516)[1:28]
n2 <- c(150,128,259,98,402,7277,308,20,34,92,306,14407,102,7295,172,3664,105,150,46695,7467,1310,1442,1346,166,161566,627,4467,22203,128,259,92,7295,105,576,150,308,14407,1442)[1:28]
labels <- 1:length(x1)

ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="IV",random=T,hk=T,fixtau=F,fixvs=F,weightiter=1,plot=T)
ratesCI(x1=x1[17],x2=x2[17],n1=n1[17],n2=n2[17],contrast="RR",stratified=T,weights="IV",random=T,hk=T,fixtau=T,fixvs=F,weightiter=1,plot=T)

U.RD <- ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="IV",random=T,hk=T,fixtau=F,fixvs=F)$estimates[,c(3)]
pbar <- (sum(x1)/sum(n1)+sum(x2)/sum(n2))/2
U.p1 <- pbar + U.RD/2
U.p2 <- pbar - U.RD/2
U.p1/U.p2
labels <- 1:length(x1)



p1.hat <- x1/n1
p2.hat <- x2/n2
p0.hat <- round((p1.hat + p2.hat)/2,3)
RD.hat <- (p1.hat - p2.hat)
RR.hat <- (p1.hat / p2.hat)
OR.hat <- p1.hat*(1-p2.hat)/(p2.hat*(1-p1.hat))

xtable(data.frame(labels,x1,n1,x2,n2,p0.hat,OR.hat,RR.hat,RD.hat,row.names=1),digits=c(0,0,0,0,0,3,3,3,3))

fround <- function(x,digits=6) paste("(",paste(format(round(x,digits=digits),nsmall=digits),collapse=", "),")",sep="")

data(cisapride)
    m.or <- metabin(event.cisa, n.cisa, event.plac, n.plac,
                    data=cisapride, sm="OR", method="Inverse",
                    studlab=study, addincr=TRUE,backtransf=FALSE)

m.rr.hakn <- metabin(event.cisa, n.cisa, event.plac, n.plac,
                     data=cisapride, sm="RR", method="Inverse",
                     studlab=study, addincr=TRUE,
                     hakn=TRUE)
# Results for log risk ratio - see Table VII in Hartung and Knapp (2001)
#
    m.rr <- metabin(event.cisa, n.cisa, event.plac, n.plac,
                    data=cisapride, sm="RR", method="Inverse",
                    studlab=study, addincr=TRUE)
res.rr <- rbind(data.frame(summary(m.rr)$fixed)[c("TE", "lower", "upper")],
                data.frame(summary(m.rr)$random)[c("TE", "lower", "upper")],
                data.frame(summary(m.rr.hakn)$random)[c("TE", "lower", "upper")])
#
row.names(res.rr) <- c("FE", "RE", "RE (HaKn)")
names(res.rr) <- c("Log risk ratio", "CI lower", "CI upper")
#
res.rr

x1 <- x1hk
x2 <- x2hk
n1 <- n1hk
n2 <- n2hk

mytab <-rbind(
#I2=c(
##RDbin=round(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD")$homog[,c(3)],0),
##RRbin=round(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR")$homog[,c(3)],0),
##OR=round(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR")$homog[,c(3)],0)
#RDbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="IV",skew=T,random=F,hk=F,fixtau=T)$Qtest["I2"],0),
#RRbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="IV",skew=T,random=F,hk=F,fixtau=T)$Qtest["I2"],0),
#OR=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="IV",skew=T,random=F,hk=F,fixtau=T)$Qtest["I2"],0)
#),
#pval=c(
##RDbin=round(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD")$homog[,c(4)],3),
##RRbin=round(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR")$homog[,c(4)],3),
##OR=round(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR")$homog[,c(4)],3)
#RDbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="IV",skew=T,random=F,hk=F,fixtau=T)$Qtest["pval.het"],3),
#RRbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="IV",skew=T,random=F,hk=F,fixtau=T)$Qtest["pval.het"],3),
#OR=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="IV",skew=T,random=F,hk=F,fixtau=T)$Qtest["pval.het"],3)
#),
ANiv=c(
RDbin=fround(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",weight="Inverse")$estimates[,c(1,3)],3),
RRbin=fround((diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",weight="Inverse")$estimates[,c(1,3)]),2),
OR=fround((diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",weight="Inverse")$estimates[,c(1,3)]),2)
),
ANmh=c(
RDbin=fround(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",weight="MH")$estimates[,c(1,3)],3),
RRbin=fround(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",weight="MH")$estimates[,c(1,3)],2),
OR=fround(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",weight="MH")$estimates[,c(1,3)],2)
),
MNmh=c(
RDbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="MH",skew=F,random=F,hk=T,fixtau=T)$estimates[,c(1,3)],3),
RRbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="MH",skew=F,random=F,hk=T,fixtau=T)$estimates[,c(1,3)],2),
OR=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="MH",skew=F,random=F,hk=T,fixtau=T)$estimates[,c(1,3)],2)
),
SCASmh=c(
RDbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="MH",skew=T,random=F,hk=T,fixtau=T)$estimates[,c(1,3)],3),
RRbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="MH",skew=T,random=F,hk=T,fixtau=T)$estimates[,c(1,3)],2),
OR=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="MH",skew=T,random=F,hk=F,fixtau=T)$estimates[,c(1,3)],2)
),
SCASiv=c(
RDbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="IV",skew=T,random=F,hk=T,fixtau=T)$estimates[,c(1,3)],3),
RRbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="IV",skew=T,random=F,hk=T,fixtau=T)$estimates[,c(1,3)],2),
OR=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="IV",skew=T,random=F,hk=F,fixtau=T)$estimates[,c(1,3)],2)
),
DL=c(
RDbin=fround(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD")$estimates[,c(4,5)],3),
RRbin=fround(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR")$estimates[,c(4,5)],2),
OR=fround(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR")$estimates[,c(4,5)],2)
),
HKSJ=c(
RDbin=fround(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",hakn=T)$estimates[,c(4,5)],3),
RRbin=fround((diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",hakn=T)$estimates[,c(4,5)]),2),
OR=fround((diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",hakn=T)$estimates[,c(4,5)]),2)
),
TDAS=c(
RDbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="IV",random=T,hk=T,fixtau=T)$estimates[,c(1,3)],3),
RRbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",dist="bin",stratified=T,weights="IV",random=T,hk=T,fixtau=T,plotmax=12)$estimates[,c(1,3)],2),
OR=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="IV",random=T,hk=T,fixtau=T)$estimates[,c(1,3)],2)
)
)
mytab



if(FALSE) {

#bizarre things happen with fixtau=F, even with iterative weights (e.g. CRASH excluding 2 extremes):
#but it seems to improve the results for the Tarone dataset.
TDASi=c(
RDbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="IV",random=T,hk=T,fixtau=F,weightiter=1,plot=T)$estimates[,c(1,3)],3),
RRbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="IV",random=T,hk=T,fixtau=F,weightiter=1,plot=T,plotmax=12)$estimates[,c(1,3)],2),
OR=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="IV",random=T,hk=T,fixtau=F,weightiter=1,plot=T)$estimates[,c(1,3)],2)
),
TDASi2=c(
RDbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="IV",random=T,hk=T,fixtau=F,weightiter=2,plot=T)$estimates[,c(1,3)],3),
RRbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="IV",random=T,hk=T,fixtau=F,weightiter=2,plot=T,plotmax=12)$estimates[,c(1,3)],2),
OR=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="IV",random=T,hk=T,fixtau=F,weightiter=2,plot=T)$estimates[,c(1,3)],2)
),
TDASi10=c(
RDbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="IV",random=T,hk=T,fixtau=F,weightiter=10,plot=T)$estimates[,c(1,3)],3),
RRbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="IV",random=T,hk=T,fixtau=F,weightiter=10,plot=T,plotmax=12)$estimates[,c(1,3)],2),
OR=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="IV",random=T,hk=T,fixtau=F,weightiter=10,plot=T)$estimates[,c(1,3)],2)
),
#TDASr=c(
#RDbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="IV",random=T,hk=T,fixtau=T,robust=T)$estimates[,c(1,3)],3),
#RRbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="IV",random=T,hk=T,fixtau=T,robust=T)$estimates[,c(1,3)],2),
#OR=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="IV",random=T,hk=T,fixtau=T,robust=T)$estimates[,c(1,3)],2)
#),
#TDAStrim=c(
#RDbin=fround(ratesCI(x1=x1[-c(3,11)],x2=x2[-c(3,11)],n1=n1[-c(3,11)],n2=n2[-c(3,11)],contrast="RD",stratified=T,weights="IV",random=T,hk=T,fixtau=T)$estimates[,c(1,3)],3),
#RRbin=fround(ratesCI(x1=x1[-c(3,11)],x2=x2[-c(3,11)],n1=n1[-c(3,11)],n2=n2[-c(3,11)],contrast="RR",stratified=T,weights="IV",random=T,hk=T,fixtau=T)$estimates[,c(1,3)],2),
#OR=fround(ratesCI(x1=x1[-c(3,11)],x2=x2[-c(3,11)],n1=n1[-c(3,11)],n2=n2[-c(3,11)],contrast="OR",stratified=T,weights="IV",random=T,hk=T,fixtau=T)$estimates[,c(1,3)],2)
#),
TDASmh=c(
RDbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="MH",random=T,hk=T,fixtau=T)$estimates[,c(1,3)],3),
RRbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="MH",random=T,hk=T,fixtau=T)$estimates[,c(1,3)],2),
OR=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="MH",random=T,hk=T,fixtau=T)$estimates[,c(1,3)],2)
)
}


mytab
##library(xtable)
xtable(mytab)

,
TDASmh=c(
RDbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="MH",random=T,hk=T,fixtau=T)$estimates[,c(1,3)],3),
RRbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="MH",random=T,hk=T,fixtau=T)$estimates[,c(1,3)],2),
OR=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="MH",random=T,hk=T,fixtau=T)$estimates[,c(1,3)],2)
),
TDASmhr=c(
RDbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="MH",random=T,hk=T,fixtau=T,robust=T)$estimates[,c(1,3)],3),
RRbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="MH",random=T,hk=T,fixtau=T,robust=T)$estimates[,c(1,3)],2),
OR=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="MH",random=T,hk=T,fixtau=T,robust=T)$estimates[,c(1,3)],2)
)


ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="IV",hk=T,fixtau=T)

ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,leftside=T,random=T,hk=T,bcf=T,contrast="RD",dist="bin",weights="IV",plot=T) #fixed
xyz <- ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,leftside=T,random=F,hk=T,bcf=T,contrast="RR",dist="bin",weights="IV") #fixed

qqnorm(xyz$Stheta)

ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,leftside=T,random=T,hk=T,bcf=T,contrast="OR",dist="bin",weights="IV",plot=T,robust=T) #random
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,leftside=T,random=T,hk=T,bcf=T,contrast="RD",dist="bin",weights="IV",plot=T) #random
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,leftside=T,random=T,hk=T,bcf=T,contrast="RR",dist="bin",weights="IV",simpletau=T) #random
ratesCI(x1=sum(x1),x2=sum(x2),n1=sum(n1),n2=sum(n2),stratified=T,delta=1,skew=F,leftside=T,random=F,hk=T,bcf=T,contrast="RR",dist="bin",weights="MH") #fixed
ratesCI(x1=x1[3],x2=x2[3],n1=n1[3],n2=n2[3],stratified=T,delta=1,skew=F,leftside=T,random=F,hk=T,bcf=T,contrast="OR",dist="bin",weights="MH") #fixed
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,leftside=T,random=F,hk=T,bcf=T,contrast="RR",dist="poi",weights="MH") #fixed
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,leftside=T,random=F,hk=T,bcf=T,contrast="RR",dist="poi",weights="IV") #random(IV)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,leftside=T,random=T,hk=T,bcf=T,contrast="RR",dist="poi",weights="MH") #random(MH)
metabin(event.e=x1,n.e=n1,event.c=x2,n.c=n2,sm="RR",backtransf=T,level=0.95,level.comb=0.95,method="MH",warn=T,addincr=F,incr=0.5,hakn=F,allstudies=F)	
metabin(event.e=x1,n.e=n1,event.c=x2,n.c=n2,sm="RR",backtransf=F,level=0.95,level.comb=0.95,method="Inverse",warn=T,addincr=F,incr=0.5,hakn=F,allstudies=F)	


metabin(event.e=x1,n.e=n1,event.c=x2,n.c=n2,sm="RR",backtransf=T,level=0.95,level.comb=0.95,method="MH",warn=T,addincr=F,incr=0.5,hakn=T,allstudies=F)	
metabin(event.e=x1,n.e=n1,event.c=x2,n.c=n2,sm="RD",backtransf=T,level=0.95,level.comb=0.95,method="MH",warn=T,addincr=F,incr=0.5,hakn=T,allstudies=F)	

(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="Inverse",level=0.95,hakn=T,contrast="RD")$estimates)
(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="Inverse",level=0.95,hakn=F,contrast="RD")$estimates)
(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="MH",level=0.95,hakn=F,contrast="RD")$estimates)
exp
(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="Inverse",0.95,hakn=F,contrast="OR")$estimates)
(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="Inverse",0.95,hakn=T,contrast="OR")$estimates)

ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,leftside=T,random=F,bcf=T,contrast="RD",dist="bin",weights="IV")

ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,leftside=T,random=F,hk=F,bcf=T,contrast="OR",dist="bin",weights="IV")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=T,leftside=T,random=F,hk=F,bcf=T,contrast="OR",dist="bin",weights="IV")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,leftside=T,random=T,hk=T,bcf=T,contrast="OR",dist="bin",weights="IV")



exp(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="MH",contrast="RD")$estimates )
exp(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="Inverse",contrast="OR",hakn=T)$estimates )
exp
(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="Inverse",contrast="RR",hakn=T)$estimates )
exp(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="MH",contrast="RR")$estimates )
(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="Inverse",contrast="RD",hakn=T)$estimates )
(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="MH",contrast="RD",hakn=T)$estimates )
diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="Inverse",contrast="RD")$homog
diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="MH",contrast="RD")$homog
exp(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="MH",contrast="RR")$estimates )
exp(diffBinconf.metabin(x1=x1,x2=x2,n1=n1,n2=n2,weight="Inverse",contrast="RR",hakn=T)$estimates )
diffBinconf.metabin.fixed(x1=x1,x2=x2,n1=n1,n2=n2,weight="MH") 
diffBinconf.meta.fixed(x1=x1,x2=x2,n1=n1,n2=n2,weight="invVar")
diffBinconf.meta.fixed(x1=x1,x2=x2,n1=n1,n2=n2,weight="MH")
ivwt <- 1/((x1/n1)*(1-(x1/n1))/n1 + (x2/n2)*(1-(x2/n2))/n2) #IV for RD
ivwt <- 1/(1/x1 + 1/x2 - 1/n1 - 1/n2)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=F,leftside=T,random=F,bcf=T,contrast="RD",dist="bin",weights="IV")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=F,leftside=T,random=F,bcf=T,contrast="RD",dist="bin",weights="MH")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=F,leftside=T,random=T,bcf=T,contrast="RD",dist="bin",weights="IV")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=F,leftside=T,random=F,bcf=T,contrast="RD",dist="bin",weights="MH")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=F,leftside=T,random=F,bcf=T,contrast="RD",dist="bin",wt=c(1,1,1,1,1,1,1,1))
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=F,leftside=T,random=T,bcf=T,contrast="RD",dist="bin",wt=c(1,1,1,1,1,1,1,1))
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=F,delta=0,skew=F,leftside=T,random=F,bcf=T,contrast="RD",dist="bin")

ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=T,leftside=T,random=F,bcf=T,contrast="RR",dist="bin",weights="IV")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=T,leftside=T,random=F,bcf=T,contrast="RR",dist="bin",weights="MH")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=T,leftside=T,random=T,bcf=T,contrast="RR",dist="bin",wt=(n1+n2))
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=T,leftside=T,random=T,hk=T,bcf=T,contrast="RR",dist="bin",weights="IV",fixtau=F)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=T,leftside=T,random=T,hk=T,bcf=T,contrast="RR",dist="bin",weights="IV",fixtau=T)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=T,leftside=T,random=T,hk=T,bcf=T,contrast="RR",dist="bin",weights="MH",fixtau=F)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=T,leftside=T,random=T,hk=T,bcf=T,contrast="RR",dist="bin",weights="MH",fixtau=T)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=T,leftside=T,random=T,hk=T,bcf=T,contrast="RR",dist="bin",weights="MH")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=T,leftside=T,random=T,hk=F,bcf=T,contrast="RR",dist="bin",weights="MH")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=T,leftside=T,random=T,hk=F,bcf=T,contrast="RR",dist="bin",weights="MH",simpletau=T)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=T,leftside=T,random=T,hk=F,bcf=T,contrast="RR",dist="bin",weights="IV")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=F,leftside=T,random=F,bcf=T,contrast="RR",dist="bin",wt=c(1,1,1,1,1,1,1,1))
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=F,leftside=T,random=T,bcf=T,contrast="RR",dist="bin",wt=c(1,1,1,1,1,1,1,1))

ratesCI(x1=sum(x1),x2=sum(x2),n1=sum(n1),n2=sum(n2),stratified=T,delta=1,skew=F,leftside=T,random=F,bcf=T,contrast="RR",dist="bin",weights="MH")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=F,delta=1,skew=F,leftside=T,random=F,bcf=T,contrast="RR",dist="bin")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,leftside=T,random=F,bcf=T,contrast="RR",dist="poi",weights="IV")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,leftside=T,random=F,bcf=T,contrast="RR",dist="poi",weights="MH")
ratesCI(x1=sum(x1),x2=sum(x2),n1=sum(n1),n2=sum(n2),stratified=T,delta=1,skew=F,leftside=T,random=F,bcf=T,contrast="OR",dist="bin",weights="MH")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,leftside=T,random=F,bcf=T,contrast="OR",dist="bin",weights="IV")

ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=T,leftside=T,random=T,hk=T,bcf=T,contrast="OR",dist="bin",weights="MH",simpletau=T)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=T,leftside=T,random=T,hk=T,bcf=T,contrast="OR",dist="bin",weights="MH",simpletau=F)


ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=T,leftside=T,random=T,bcf=T,contrast="OR",dist="bin",weights="IV")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=T,leftside=T,random=T,hk=F,bcf=T,contrast="OR",dist="bin",weights="MH")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=F,delta=1,skew=F,leftside=T,random=F,bcf=T,contrast="OR",dist="bin")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=T,leftside=T,random=F,bcf=T,contrast="RR",dist="bin")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=T,leftside=F,random=F,bcf=T,contrast="RD")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,leftside=F,random=F,bcf=T,contrast="OR")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0.5,skew=T,leftside=F,random=F,bcf=T,contrast="OR")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0.5,skew=F,leftside=F,random=F,bcf=T,contrast="RR")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0.5,skew=F,leftside=F,random=F,bcf=T,contrast="RD",weights="CMH")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0.00000,skew=F,leftside=F,random=T)
sapply(c(-0.1,-0.01,-0.00002,-0.00001,-0.0000000001,0,0.0000001,0.00001,0.00002,0.01,0.1), function(x) score(x, x1, x2, n1, n2,stratified=T,random=T,skew=F)$pval)

ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0.5,skew=F,leftside=F,random=F,bcf=T,contrast="RR")

ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=T,leftside=F,random=F,bcf=T,weights="CMH")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=F,leftside=F,dist="poi")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=T,leftside=F)

#fabricated data with unequal stratum sizes to see skewness correction
ratesCI(x1=3*x1[12:14],x2=x2[12:14],n1=3*n1[12:14],n2=n2[12:14],stratified=T,delta=0,skew=F,leftside=F,random=F,bcf=T)
ratesCI(x1=3*x1[12:14],x2=x2[12:14],n1=3*n1[12:14],n2=n2[12:14],stratified=T,delta=0,skew=T,leftside=F,random=F,bcf=T)

x1<-7; x2<-1; n1<-n2<-18
x1<-29; x2<-11; n1<-55; n2<-11

#Newcombe-ish example for publication
x1<-1; x2<-1; n1<-56; n2<-29

fround <- function(x,digits=6) paste("(",paste(format(round(x,digits=digits),nsmall=digits),collapse=", "),")",sep="")

#Fagerland-Newcombe example for publication
#output as latex table format
x1<-7; x2<-1; n1<-18; n2<-18;
x1<-53; x2<-24; n1<-65; n2<-73;

#selected examples from CRASH meta-analysis
x1<-16; x2<-16; n1<-67; n2<-28; #Faupel
x1<-114; x2<-38; n1<-201; n2<-74; #Pitts
x1<-33; x2<-21; n1<-68; n2<-62; #Dearden
x1<-26; x2<-13; n1<-49; n2<-27; #Cooper 1979

x1<-9; x2<-13; n1<-17; n2<-18; #Ransohoff 1972
x1<-34; x2<-7; n1<-72; n2<-16; #Giannotta 1984
x1<-34; x2<-4*7; n1<-72; n2<-4*16; #Giannotta 1984 as reported by CRASH
x1<-13; x2<-5; n1<-98; n2<-54; #Stubbs 1989

ratesCI(x1=34,x2=7,n1=72,n2=16,dist="bin",contrast="RD",skew=T,plot=T) #SCAS

ratesCI(x1=34,x2=7,n1=72,n2=16,dist="bin",contrast="RD",skew=F,plot=T) #MN






install.packages("xtable")
library(xtable)
fround <- function(x,digits=6) paste("(",paste(format(round(x,digits=digits),nsmall=digits),collapse=", "),")",sep="")
alpha<-0.05


ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0.5,skew=F,leftside=T,random=F,contrast="RR")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0.5,skew=F,leftside=T,random=T,contrast="RR")

(x1/n1)/(x2/n2)


score.RD(-0.3,x1,x2,n1,n2,stratified=T,bcf=T,contrast="RD",dist="bin",cc="N",alpha=0.05,skew=F,random=T)

#made up heterogeneity example:
x1 <- c(40,35,45)
x2 <- c(35,30,20)
n1 <- c(50,50,50)
n2 <- c(50,50,50)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=F,leftside=T,random=F)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=F,leftside=T,random=T)
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=0,skew=T,leftside=T,random=T)

ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,leftside=F,random=F,contrast="RR")
ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=F,leftside=F,random=T,contrast="RR")


ratesCI(0,0,50,50,contrast="RD",dist="bin",skew=T) #OK
ratesCI(0,0,1,1,contrast="RD",dist="poi",skew=F) #problem with upper - fixed by not allowing V=0
ratesCI(0,0,50,50,contrast="RR",dist="bin") #force upper=inf & MLE=NA - done
ratesCI(1,1,5,5,contrast="RR",dist="poi") #problem with MLE, should be NA
ratesCI(1,2,50,50,contrast="OR",dist="bin") #force upper=inf & MLE=NA - done

ratesCI(c(0,1),c(0,2),c(50,50),c(50,50),contrast="RR",dist="poi",skew=T) #force upper=inf & MLE=NA - done

x1 <- 0; x2 <- 0; n1 <- 1; n2 <- 1
myrange <- seq(-0.999,0.999,0.001)
myrange <- seq(0.01,10,0.1)
#sc <- sapply(myrange,function(x) score(x,x1,x2,n1,n2,skew=T,random=F,stratified=F,contrast="RD",bcf=T,dist="bin")$zscore)


x1 <- ss1[,99]
x2 <- ss2[,99]
n1 <- nn1
n2 <- nn2

pointest <- ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,stratified=T,delta=1,skew=T,leftside=T,random=T,hk=T,bcf=T,contrast="OR",dist="bin",weights="IV",fixtau=F)$estimates[2]
tfix <-  score(pointest,x1,x2,n1,n2,skew=T,random=T,hk=T,stratified=T,contrast="OR",bcf=T,dist="bin",weights="IV",tau2=NULL)$tau2
sc <- sapply(myrange,function(x) score(x,x1,x2,n1,n2,skew=T,random=T,hk=T,stratified=T,contrast="OR",bcf=T,dist="bin",weights="IV",tau2=tfix)$zscore)
#sc <- sapply(myrange,function(x) score(x,x1,x2,n1,n2,skew=T,random=F,hk=F,stratified=T,contrast="RR",bcf=T,dist="bin",weights="IV",tau2=tfix)$zscore)
plot(myrange,sc,type="l",ylim=c(-13,14))

nm <- sapply(myrange,function(x) score(x,x1,x2,n1,n2,skew=T,random=T,hk=T,stratified=T,contrast="RR",bcf=T,dist="bin",weights="IV",tau2=tfix)$numer)
vd <- sapply(myrange,function(x) score(x,x1,x2,n1,n2,skew=T,random=T,hk=T,stratified=T,contrast="RR",bcf=T,dist="bin",weights="IV",tau2=tfix)$Vdot)
plot(myrange,nm/sqrt(vd),type="l",ylim=c(-13,4))
qh <- sapply(myrange,function(x) score(x,x1,x2,n1,n2,skew=T,random=T,hk=T,stratified=T,contrast="RR",bcf=T,dist="bin",weights="IV",tau2=tfix,Vfix=Vfix)$qhat)
plot(myrange,qh,type="l",ylim=c(0,4))

plot(myrange,sqrt(vd),type="l",ylim=c(0,4))
#sc <- sapply(myrange,function(x) score(x,x1,x2,n1,n2,skew=T,random=T,hk=T,stratified=T,contrast="RR",bcf=T,dist="bin",weights="MH",tau2=NULL)$numer)
#sc <- sapply(myrange,function(x) score(x,x1,x2,n1,n2,skew=T,random=T,hk=T,stratified=T,contrast="RR",bcf=T,dist="bin",weights="MH",tau2=NULL)$zscore)
#sc <- sapply(myrange,function(x) score(x,x1,x2,n1,n2,skew=F,random=F,stratified=T,contrast="OR",bcf=T,dist="bin",weights="IV")$zscore)
abline(h=1.96*c(-1,1))
abline(h=qt(0.975,1)*c(-1,1))
abline(v=c(1.9,2.6))

score(0.005,x1,x2,n1,n2,skew=T,random=F,stratified=F,contrast="RD",bcf=F,dist="poi")

}







#code graveyard
#zlimv <- function(x1,x2,n1,n2,bcf=T,alpha=0.05,uplow="low",precis=4,skew=T,cc="N"){
bisect <- function(ftn,contrast, precis, max.iter = 100, d1=0, d2=1,uplow="low", logscale=F){
	if (uplow=="low") dir <- (-1) else dir <- (1)
	tiny <- (10^-(precis))/2
   #find min/max theta for which |score| < |threshold|
   threshold <- dir*qnorm(alpha/2)
   proot <- d1
   best <- rep(dir,length(proot))
   dp <- abs(proot - best)
   niter <- 1
   while(niter <= max.iter && any(dp>tiny)){
     dp <- 0.5*dp
     low2 <- pmax(-1,pmin(1,round(proot+dir*dp,10)))  #rounding avoids machine precision problem with, e.g. 7/10-6/10
     score <- ftn(low2)
     check <- (dir*score>dir*threshold) | is.na(score) #score=NA only happens when |p1-p2|=1 and |theta|=1, in which case proot==best anyway
     proot[check] <- low2[check]
     best[!check] <- low2[!check]
     niter <- niter+1 
   }
   proot	
}



#vectorized secant function including log scale calculation for RR - see Gart & Nam 1988
secant <- function(ftn, tol, max.iter = 100, d1=1, d2=1, logscale=F){
  d3 <- d2
  iter <- 0
  f2 <- ftn(d2)
#  while( (max(abs(f2)) > tol) && (max.iter > iter) && any(d1!=d2) ){
  while( any(abs(f2)>tol | d1==d2) && (max.iter > iter) ){
    if(logscale==TRUE) { d3 <- exp(log(d2) - f2 * (log(d2) - log(d1)) / (f2 - ftn(d1)))
#	} else d3 <- pmin(1,pmax(-1,d2 - f2 * (d2 - d1) / (f2 - ftn(d1))))
	} else d3 <- atan(tan(d2*pi/2) - f2 * (tan(d2*pi/2) - tan(d1*pi/2)) / (f2 - ftn(d1)))*2/pi
#	} else d3 <- d2 - f2 * (d2 - d1) / (f2 - ftn(d1))
    d1 <- d2
    d2[!is.na(d3)] <- d3[!is.na(d3)]
    f2 <- ftn(d2)
    iter <- iter + 1
  }
#  if(max((abs(f2))) > tol && iter >= max.iter d1==d2){

  if(any((abs(f2) > tol & d1!=d2))  && iter >= max.iter){
     print("Failed to Find root")
  } else return(d2)
}


	#identify crude estimate of CI for use in secant optimisation.  For stratified data, we use crude pooling here
	if(contrast=="RD") {
		crude.CI <- crude.RD(x1,x2,n1,n2,alpha=alpha,pooled=stratified)
	} else if(contrast=="RR") crude.CI <- crude.RR(x1,x2,n1,n2,alpha=alpha,pooled=stratified)
	#dimnames(crude.CI)[[1]]<-NULL

#NB should really add approximate Poisson CI method here, but BL may suffice
crude.RD <- function(x1,x2,n1,n2,alpha=0.05,cc=FALSE,pooled=FALSE,dist="bin") {
	#vectorised function for Brown-Li CI for first approximation in secant method
	if(pooled==TRUE) {x1 <- sum(x1); x2 <- sum(x2); n1 <- sum(n1); n2 <- sum(n2)}
	p1.hat <- (x1+0.5)/(n1+1)
	p2.hat <- (x2+0.5)/(n2+1)
	lcl <- pmax(-1,p1.hat-p2.hat - qnorm(1-alpha/2)*sqrt((p1.hat*(1-p1.hat)/(n1) + p2.hat*(1-p2.hat)/(n2))))
	ucl <- pmin(1,p1.hat-p2.hat + qnorm(1-alpha/2)*sqrt((p1.hat*(1-p1.hat)/(n1) + p2.hat*(1-p2.hat)/(n2))))
	cbind(lcl,ucl)
}
#Rothman Greenland method for Poisson RD
crude.RDID <- function(x1,x2,n1,n2,alpha=0.05) {
	if (length(n1)==1 & length(x1)>1) n1 <- rep(n1,length(x1)) #in case x1,x2 are input as vectors but n1,n2 are not
	if (length(n2)==1 & length(x2)>1) n2 <- rep(n2,length(x2))
	z <- qnorm(1-alpha/2)
	lcl <- (x1/n1-x2/n2)-z*sqrt(x1/n1^2+x2/n2^2)
	ucl <- (x1/n1-x2/n2)+z*sqrt(x1/n1^2+x2/n2^2)
	(cbind(lcl,ucl))
}	
#Rothman Greenland method for use as start point for secant iteration
crude.RR <- function(x1,x2,n1,n2,alpha=0.05,pooled=FALSE) {
	if(pooled==TRUE) {x1 <- sum(x1); x2 <- sum(x2); n1 <- sum(n1); n2 <- sum(n2)}
	z <- qnorm(1-alpha/2)
	x1[x1==0] <- 0.5
	x2[x2==0] <- 0.5
	lcl <- log(x1*n2/(x2*n1))-z*sqrt(1/x1+1/x2)
	ucl <- log(x1*n2/(x2*n1))+z*sqrt(1/x1+1/x2)
	exp(cbind(lcl,ucl))
}
	
}


