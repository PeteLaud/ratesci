# Confidence intervals applying the MOVER method across different contrasts and distributions
# using Jeffreys intervals (& other more general Beta and Gamma priors) instead of Wilson

JeffreysCI<-function(
  x,
  n,
  ai=0.5,
  bi=0.5,
  gamma=0, #continuity correction. gamma=0.5 gives Clopper-Pearson 'exact' (conservative) interval
  level=0.95,
  dist="bin",
  adj=FALSE,  #option for exact binomial interval for zero cell counts as suggested by Brown et al
  exact=F,
  ...
  ) {
	if(dist=="bin") {
	  CI.lower <- ifelse(x==0,0,qbeta((1-level)/2, x + (ai-gamma), n - x + (bi+gamma)))
	  CI.upper <- ifelse(x==n,1,qbeta(1 - (1-level)/2, x + (ai+gamma), n - x + (bi-gamma)))
		if(adj==TRUE){ 
			CI.lower[x==n] <- ((1-level)/2)^(1/n)[x==n]
			CI.upper[x==0] <- 1-((1-level)/2)^(1/n)[x==0]
		}
	} else if(dist=="poi") {
		#Jeffreys prior for Poisson rate uses gamma distribution, as defined in Li et al.
		CI.lower <- ifelse(x==0,0,qgamma(alpha/2, x+(ai-gamma), scale=1/n))
		CI.upper <- qgamma(1-alpha/2, (x+(ai+gamma)), scale=1/n)
	}
	CI <- cbind(Lower=CI.lower, Upper=CI.upper)
	CI
}

quadroot <- function(a, b, c)
{
	# GET ROOTS OF A QUADRATIC EQUATION
	r1x <- ( - b + sqrt(b^2 - 4 * a * c))/(2 * a)
	r2x <- ( - b - sqrt(b^2 - 4 * a * c))/(2 * a)
	r1<-pmin(r1x,r2x)
	r2<-pmax(r1x,r2x)
	cbind(r1,r2)
}

ratesCI.MOVER <- function(
  x1, 
  x2, 
  n1, 
  n2,
  a1=0.5, 
  b1=0.5, 
  a2=0.5, 
  b2=0.5, 
  cc=0,
  level = 0.95,
  dist="bin",
  contrast="RD",
  type="Jeff",  # "Jeff", "Exact", "Wilson"
  ...
  ){
	alpha <- 1-level
	z <- qnorm(1 - alpha/2)
	if(as.character(cc)=="TRUE") cc <- 0.5
	
	if (length(n1)==1 & length(x1)>1) n1=rep(n1,length(x1)) #in case x1,x2 are vectors but n1,n2 are not
	if (length(n2)==1 & length(x1)>1) n2=rep(n2,length(x1))

	if(contrast=="OR") { #special cases for OR handled as per Fagerland & Newcombe Table III
		special <- ((x2==n2) | (x1==n1))
		xx <- x2
		x2[special] <- (n1-x1)[special]
		x1[special] <- (n2-xx)[special]
		nx <- n1
		n1[special] <- n2[special]
		n2[special] <- nx[special]		
	}

	p1hat<-x1/n1
	p2hat<-x2/n2

	if(contrast=="OR" && dist!="bin") {
		print("WARNING: Odds Ratio must use dist='bin'")
		dist <- "bin"
	}

	if(type=="Jeff") { #MOVER-J, including optional "continuity correction" gamma
		j1<-JeffreysCI(x1,n1,ai=a1,bi=b1,gamma=cc,alpha=alpha,dist=dist,adj=paste((contrast=="OR"))); 
		j2<-JeffreysCI(x2,n2,ai=a2,bi=b2,gamma=cc,alpha=alpha,dist=dist,adj=paste((contrast=="OR")));
	} else if(type=="exact") { #MOVER-E based on Clopper-Pearson exact intervals
		j1<-JeffreysCI(x1,n1,ai=a1,bi=b1,gamma=0.5,alpha=alpha,dist=dist,adj=paste((contrast=="OR"))); 
		j2<-JeffreysCI(x2,n2,ai=a2,bi=b2,gamma=0.5,alpha=alpha,dist=dist,adj=paste((contrast=="OR")));
	} else { #or use Wilson intervals as per Newcombe 1988
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
	} else if(contrast=="RR") {
		#From Donner & Zou 2012 / Li et al
		lower <- (p1hat*p2hat - Re(sqrt(as.complex((p1hat*p2hat)^2 - l1*(2*p2hat-u2)*(u2*(2*p1hat-l1)))))) / (u2*(2*p2hat-u2))
		upper <- (p1hat*p2hat + Re(sqrt(as.complex((p1hat*p2hat)^2 - u1*(2*p2hat-l2)*(l2*(2*p1hat-u1)))))) / (l2*(2*p2hat-l2))
		upper[x2==0] <- Inf
	}
	CI <- cbind(Lower =  lower, Upper = upper)
	CI
}

ratesCI.MOVER(5,0,56,29,type="Wilson",contrast="RD") #Newcombe 1998
ratesCI.MOVER(5,0,56,29,type="exact",contrast="RD")
ratesCI.MOVER(5,0,56,29,type="Jeff",contrast="RD",cc=0.25)








if(F) {



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
