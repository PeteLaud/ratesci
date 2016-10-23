#Attempt to apply the MOVER method across different contrasts and distributions
# using Jeffreys intervals instead of Wilson

Jeff<-function(x,n,k=0.5,alpha=0.05,dist="bin",adj="FALSE",exact=F) {
	if(dist=="bin") {
	#	CI.lower <- qbeta(alpha/2, x + k, n - x + 1- k)
	#	CI.upper <- qbeta(1 - alpha/2, x + 1 - k, n - x + k)
		CI.lower <- qbeta(alpha/2, x + k, n - x + k)
		CI.upper <- qbeta(1 - alpha/2, x + k, n - x + k)
		CI.upper[x==n] <- 1
		if(adj==TRUE){ 
			CI.lower[x==n] <- (alpha/2)^(1/n)[x==n]
			CI.upper[x==0] <- 1-(alpha/2)^(1/n)[x==0]
		}
		if(exact==T) {
	#		CI.lower <- ifelse(x==0,0,qbeta(alpha/2, x , n - x + 1))
	#		CI.upper <- ifelse(x==n,1,qbeta(1 - alpha/2, x + 1, n - x))
			CI.lower <- ifelse(x==0,0,qbeta(alpha/2, x + (0.5-k), n - x + (0.5+k)))
			CI.upper <- ifelse(x==n,1,qbeta(1 - alpha/2, x + (0.5+k), n - x + (0.5-k)))
		}
	} else if(dist=="poi") {
		#From Li et al.
		CI.lower <- qgamma(alpha/2, x + k, n)
		CI.upper <- qgamma(1-alpha/2, x + k, n)
		if(exact==T) {
	#		CI.lower <- qchisq(alpha/2, 2*(x+(0.5-k)))/(2*n)
	#		CI.upper <- qchisq(1-alpha/2, 2*(x+(0.5+k)))/(2*n)
			CI.lower <- qgamma(alpha/2, x+(0.5-k), scale=1/n)
			CI.upper <- qgamma(1-alpha/2, (x+(0.5+k)), scale=1/n)
		}
		
	}
	CI.lower[x==0] <- 0
	CI <- cbind(Lower=CI.lower, Upper=CI.upper)
	CI
}
Jeff(14,400,dist="bin",exact=T)

Jeff(5,56,dist="poi",exact=T)
qchisq(0.025,2*5)/(2*56)
qgamma(0.025,5,scale=1/56)


quadroot <- function(a, b, c)
{
	# GET ROOTS OF A QUADRATIC EQUATION
	r1x <- ( - b + sqrt(b^2 - 4 * a * c))/(2 * a)
	r2x <- ( - b - sqrt(b^2 - 4 * a * c))/(2 * a)
	r1<-pmin(r1x,r2x)
	r2<-pmax(r1x,r2x)
	cbind(r1,r2)
}

x1=24; x2=53; n1=73;n2=65; contrast="OR"; dist="bin"

MOVER <- function(x1, x2, n1, n2,k=0.5, level = 0.95,dist="bin",contrast="RD",type="Jeff"){
	alpha <- 1-level
	z <- qnorm(1 - alpha/2)

	if (length(n1)==1 & length(x1)>1) n1=rep(n1,length(x1)) #in case x1,x2 are vectors but n1,n2 are not
	if (length(n2)==1 & length(x1)>1) n2=rep(n2,length(x1))
	if (length(x2)==1 & length(x1)>1) x2=rep(x2,length(x1))
	
	if(contrast=="OR") {
		#need to vectorize!
		special <- ((x2==n2) | (x1==n1))
		xx <- x2; x2[special] <- (n1-x1)[special]; x1[special] <- (n2-xx)[special]; 
		nx <- n1; n1[special] <- n2[special]; n2[special] <- nx[special]		
#		if(x2==n2 || x1==n1) {xx <- x2; x2 <- n1-x1; x1 <- n2-xx; nx <- n1; n1 <- n2; n2 <- nx }
#		if(x1==n1) {x1 <- n2-x2; nx <- n1; n1 <- n2; n2 <- nx; x2 <- 0}
	}
	
	p1hat<-x1/n1
	p2hat<-x2/n2

	if(contrast=="OR" && dist!="bin") {
		print("WARNING: Odds Ratio mst use dist='bin'")
		dist <- "bin"
	}

#	l<-Jeff(x1,n1,k=k,alpha=alpha,dist=dist)
#	u<-Jeff(x2,n2,k=k,alpha=alpha,dist=dist)
	if(type=="Jeff") {
		j1<-Jeff(x1,n1,k=k,alpha=alpha,dist=dist,adj=paste((contrast=="OR"))); 
		j2<-Jeff(x2,n2,k=k,alpha=alpha,dist=dist,adj=paste((contrast=="OR")));
	} else if(type=="exact") {
		j1<-Jeff(x1,n1,k=k,exact=T,alpha=alpha,dist=dist,adj=paste((contrast=="OR"))); 
		j2<-Jeff(x2,n2,k=k,exact=T,alpha=alpha,dist=dist,adj=paste((contrast=="OR")));
	} else {
		j1 <- quadroot(a=1 + z^2/(n1), b=-(2*p1hat + z^2/(n1)), c=p1hat^2)
		j2 <- quadroot(1 + z^2/(n2), -(2*p2hat + z^2/(n2)), p2hat^2)
	}
	l1 <- j1[,1]; u1 <- j1[,2]
	l2 <- j2[,1]; u2 <- j2[,2]


	if(contrast=="RD") {
		#From Newcombe 1998
		lower <- p1hat - p2hat - Re(sqrt(as.complex((p1hat - l1)^2 + (u2 - p2hat)^2)))
		upper <- p1hat - p2hat + Re(sqrt(as.complex((u1 - p1hat)^2 + (p2hat - l2)^2)))
	} else if(contrast=="OR") {
		#From Fagerland & Newcombe 2013
		q1hat <- p1hat/(1-p1hat)
		q2hat <- p2hat/(1-p2hat)
		L1 <- l1/(1-l1); U1 <- u1/(1-u1)
		L2 <- l2/(1-l2); U2 <- u2/(1-u2)
		lower <- pmax(0,(q1hat*q2hat - Re(sqrt(as.complex((q1hat*q2hat)^2 - L1*U2*(2*q1hat-L1)*(2*q2hat-U2))))) / (U2*(2*q2hat-U2)))
		upper <-        (q1hat*q2hat + Re(sqrt(as.complex((q1hat*q2hat)^2 - U1*L2*(2*q1hat-U1)*(2*q2hat-L2))))) / (L2*(2*q2hat-L2))
		upper[x2==0] <- Inf
		lower[(x1==0 & x2==n2) | (x1==n1 & x2==0)] <- 0
		upper[(x1==0 & x2==n2) | (x1==n1 & x2==0)] <- Inf
#		lower[x2==n2] <- 0
	} else if(contrast=="RR") {
		#From Donner & Zou 2012 / Li et al
		lower <- (p1hat*p2hat - Re(sqrt(as.complex((p1hat*p2hat)^2 - l1*(2*p2hat-u2)*(u2*(2*p1hat-l1)))))) / (u2*(2*p2hat-u2))
		upper <- (p1hat*p2hat + Re(sqrt(as.complex((p1hat*p2hat)^2 - u1*(2*p2hat-l2)*(l2*(2*p1hat-u1)))))) / (l2*(2*p2hat-l2))
		upper[x2==0] <- Inf
		lower[x1==0] <- 0
	}
	
	CI <- cbind(Lower =  lower, Upper = upper)
	CI
}

MOVER(5,0,56,29,type="exact",contrast="RD")
MOVER(5,0,56,29,type="exact",contrast="RD",k=0)
MOVER(5,0,56,29,type="Jeff",contrast="RD")
MOVER(5,0,56,29,type="Jeff",contrast="RD",k=0.25)



CI.BC <- function(x1,x2,n1,n2,level=0.95,k=0.5) {
	#conditional confidence interval for RR, based on Jeffreys interval for single proportion
	#(could also base on skewness-corrected score interval for a single proportion)
	#as per Barker & Caldwell paper on ratio of Poisson rates - may also apply for binomial with rare events
#	cip <- ratesCI(x1,x1,n1=(x1+x2),(x1+x2),level=level,contrast="p",dist="bin",stratified=stratified,wt=n1+n2,random=random,hk=hk) #NB dist=bin here even in the poisson case, because x1/(x1+x2)~bin
#	cij <- Jeff(x=x1,n=(x1+x2),k=k,alpha=(1-level),dist="bin")
	l1 <- qbeta((1-level)/2, x1 + k, x2 + k)
	u1 <- qbeta(1 - (1-level)/2, x1 + k, x2 + k)
	l1[x1==0] <- 0
	u1[x2==0] <- 1
	lower=(n2/n1)*l1/(1-l1)
	upper=(n2/n1)*u1/(1-u1)
#	lower=l1/(1-l1)
#	upper=u1/(1-u1)
	cbind(lower,upper)
}

CI.BC(0,1,45,15)



if(F) {

xs <- expand.grid(0:30,0:30)
x1 <- xs[,1]
x2 <- xs[,2]

BarkerCI(31,0,30,30,k=1)
BarkerCI(x1,x2,10,10,k=0.5)
MOVER(x1,x2,30,30,k=0.5)
MOVER(x1,x2,30,30,k=1)

Jeff(5,5,k=1)
k=1
qbeta(alpha/2, 5 + k, 0 + k)

BarkerCI(c(10,8),c(5,4),c(1000,1000),c(1000,1000),k=0.5)
BarkerCI(c(10,8),c(5,4),c(1000,1000),c(1000,1000),k=1)



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
