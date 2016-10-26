#Need to add output of per-stratum CIs for stratified method
#rm(list=ls())

ratesCI.score <- function(

	#Score-based confidence intervals for the rate difference (RD) or ratio (RR), for independent binomial or Poisson rates, 
	#including options for bias correction (from MN), skewness correction (from GN) and continuity correction 
	#Written by Pete Laud, March 2014
	#executable in either S-plus or R
	#see ratesCI.qc for validation code

	#refs:	MN = Miettinen & Nurminen, Statistics in Medicine 1985; 4:213-226
	# 		FM = Farrington & Manning, Statistics in Medicine 1990; 9:1447-1454
	# 		GN = Gart & Nam, Biometrics 1990; 46(3):637-643

	#This function is vectorised in x1,x2,n1,n2 
	
	x1,
	x2=NULL,  		#x1,x2: vectors of number of events in group 1 & group 2 respectively
	n1,
	n2=NULL,		#n1,n2: vectors of sample sizes (or exposure times for Poisson rates) in each treatment group
	dist="bin",		#dist: whether the data represent binomial ("bin") or Poisson ("poi") rates
	contrast="RD",	#contrast: comparative parameter ("RD"=risk/rate difference, "RR"=risk/rate ratio, "OR"=odds ratio) 
	level=0.95,		#level: confidence level (default 0.95)
	skew=TRUE,		#skew: skewness correction, generalised from Gart & Nam and Laud & Dane (TRUE/FALSE) (default=TRUE)
	bcf=TRUE,			#bcf: whether to include bias correction factor in the score denominator (applicable to dist="bin" only) (default=TRUE). 
					#		NB: bcf=FALSE option is really only included for legacy validation against previous published methods (i.e. Gart & Nam)
	cc=0,			#cc: continuity correction (0=none, TRUE=0.5=conservative, intermediate value for a compromise)  (default=0)
	delta=NULL,		#delta: non-inferiority margin, for an optional one-sided test of non-zero null hypothesis eg. H0: theta<=delta (default=NULL)
					#		NB: 1-sided p will be <0.025 iff 2-sided 95% CI excludes delta. NB: can also be used for a superiority test by setting delta=0
	precis=6,		#precis: required level of precision (number of dps)
	plot=FALSE,		#plot: whether to output plot of the score function
	plotmax=1000,	

	stratified=FALSE,
	weights="IV",	#weighting strategy: "MH" "IV" "MN"
	wt=NULL,		#optional user-specified weights
	tdas=FALSE,	
	#warn=TRUE,
	
	...)
	
{ 

#	would like to have inputs checked in a separate subroutine, but not sure how
#	checkinputs(x1=x1,x2=x2,n1=n1,n2=n2,dist=dist,contrast=contrast,level=level,bcf=bcf,skew=skew,cc=cc,delta=delta,
#	precis=precis,plot=plot,plotmax=plotmax,stratified=stratified,weights=weights,wt=wt,tdas=tdas,warn=warn)
	
	if(!is.numeric(c(x1, n1, x2, n2, delta))) {
		print("Non-numeric inputs!")
		stop()
	}
	if(any(c(x1, n1, x2, n2)<0)) {
		print("Negative inputs!")
		stop()
	}	
	if(!is.null(delta)) {	
		if(contrast=="RD") {
			if(dist=="bin" && (delta< -1 || delta>1)) {
				print("Impossible delta!")
				stop()			
			}
		}
	}
	if(as.character(cc)=="TRUE") cc <- 0.5
	#etc
		
	nstrat <- length(x1)
	#in case x1,x2 are input as vectors but n1,n2 are not
	if (length(n1)==1 && nstrat>1) n1 <- rep(n1,nstrat) 
	if (length(n2)==1 && nstrat>1) n2 <- rep(n2,nstrat)
	if (dist=="poi" && contrast=="OR") {
		print("Odds ratio not applicable to Poisson rates")
		stop()
	}
	#check for empty strata, and for x1<=n1, x2<=n2
	if (stratified) { #for stratified calculations, remove any strata with no observations 
		empty.strat <- ((n1==0|n2==0) | (contrast %in% c("RR","OR") & (x1==0 & x2==0)) | (contrast=="OR" & (x1==n1 & x2==n2)))
		x1 <- x1[!empty.strat]
		x2 <- x2[!empty.strat]
		n1 <- n1[!empty.strat]
		n2 <- n2[!empty.strat]
#		if(sum(empty.strat)>0) print("Warning: at least one stratum contributed no information and was removed")

	#for double-zero cells for RD, add 0.5 to avoid 100% weight with IV weights. Same for double-100% cells for RR & RD
	#NB this should only be necessary for weights="IV"
		if(weights=="IV") {
			zeroRD <-  (contrast %in% c("RD") & (x1==0 & x2==0)) 
			x1 <- x1+(zeroRD*0.5)
			x2 <- x2+(zeroRD*0.5)
			n1 <- n1+(zeroRD*1)
			n2 <- n2+(zeroRD*1)
			fullRD <-  (contrast %in% c("RD","RR") & (x1==n1 & x2==n2)) 
			x1 <- x1+(fullRD*0.5)
			x2 <- x2+(fullRD*0.5)
			n1 <- n1+(fullRD*1)
			n2 <- n2+(fullRD*1)
		}
	}	
	nstrat <- length(x1) #update nstrat after removing empty strata

# Put these warnings back later	
	if(stratified==TRUE && nstrat<=1) {stratified <- FALSE; tdas=FALSE 
		print("Warning: only one stratum!")
	}
	if(nstrat==0) {
		print("Warning: no data!"); 
		if(contrast %in% c("RR","OR")) {x1 <- x2 <- 0; n1 <- n2 <- 10}
	}
	#quick fix for x>n for Binomial case only (instead of outputting an error message)
	if(dist=="bin") {
		if(any(x1>n1+0.001)||any(x2>n2+0.001)) print("Warning: At least one X>n substituted with X=n")
		x1[x1>n1] <- n1; x2[x2>n2] <- n2
	}

	p1hat <- x1/n1; p2hat <- x2/n2
	
	myfun <- function(theta,randswitch=tdas,ccswitch=cc) {
		scoretheta(theta,x1=x1,x2=x2,n1=n1,n2=n2,bcf=bcf,contrast=contrast,dist=dist,stratified=stratified,wt=wt,weights=weights,
		tdas=randswitch,skew=skew,cc=ccswitch)$score  #NOTE point estimate is obtained without applying continuity correction - essentially the MLE is a 0% CI so alpha=1
	}

	point <- point.FE <- bisect(ftn=function(theta) myfun(theta,randswitch=FALSE,ccswitch=0)-0,contrast=contrast,dist=dist,precis=precis+1,uplow="low")  #fixed effects point estimate taken with no cc
#	point[scoretheta(x1/n1-x2/n2,x1,x2,n1,n2,skew=skew,cc=cc,dist=dist,contrast=contrast)==0] <- 0  #needs developing further - idea is if the score at the point estimate is zero, then that should be the MLE
	if(stratified==TRUE) point <- bisect(ftn=function(theta) myfun(theta,randswitch=tdas,ccswitch=0)-0,contrast=contrast,dist=dist,precis=precis+1,uplow="low") 	#random effects point estimate if required

	#fix some extreme cases with zero counts 
	if(stratified) {
		if(skew==FALSE & contrast %in% c("RR","OR")) point[sum(x2)>0 & sum(x1)==0] <- 0   
		if(skew==FALSE & contrast %in% c("RR","OR")) point[sum(x1)>0 & sum(x2)==0] <- Inf  
#		if(contrast %in% c("OR")) point[sum(x1)==sum(n1) & sum(x2)>0] <- Inf
	} else {
		if(contrast %in% c("RR","OR")) { ## & skew==FALSE)  
			point[x1==0 & x2==0] <- NA
			point[(x1>0 & x2==0) && skew==FALSE] <- Inf
		}
		if(contrast == "OR") { 
			point[x1==n1 & x2==n1] <- NA
			point[(x1==n1 & x2>0) && skew==FALSE] <- Inf
		}
		if(contrast == "RD") { 
			point[x1==n1 & x2==n1] <- 0
			point[x1==0 & x2==0] <- 0
		}
	}

	#identify certain quantities evaluated at the maximum likelihood point estimate for theta (incorporating random effects)
	at.MLE <- scoretheta(point,x1,x2,n1,n2,bcf=bcf,contrast=contrast,dist=dist,stratified=stratified,weights=weights,wt=wt,tdas=tdas,skew=skew,cc=cc)
	Stheta.MLE <- at.MLE$Stheta
	p1d.MLE <- at.MLE$p1d
	p2d.MLE <- at.MLE$p2d
	wt.MLE <- at.MLE$wt

	if (stratified==TRUE && nstrat>1) {
		#if stratified=T, perform homogeneity test	
		#NB heterogeneity test should be based on tdas=FALSE, but we also want to report the updated weights used when tdas=T
		#requires certain quantities evaluated at the fixed effects maximum likelihood point estimate for theta, for heterogeneity test
		at.FE <- scoretheta(point.FE,x1,x2,n1,n2,bcf=bcf,contrast=contrast,dist=dist,stratified=stratified,weights=weights,wt=wt,tdas=FALSE,skew=skew,cc=cc)
		wt.FE <- at.FE$wt
		tau2.FE <- at.FE$tau2
		Q.each <- at.FE$Q.i
		Q.FE <- at.FE$Q
		I2 <- max(0,100*(Q.FE-(nstrat-1))/Q.FE)
		pval.het <- 1-pchisq(Q.FE,nstrat-1)
		#Qualitative interaction test to be added...
		#

		p1hat.w <- sum(wt.MLE*x1/n1)/sum(wt.MLE) #as per M&N p218 (little r), but actually no longer needed for point estimate. 
		p2hat.w <- sum(wt.MLE*x2/n2)/sum(wt.MLE)
		p1d.w <- sum(wt.MLE*at.MLE$p1d)/sum(wt.MLE) #as per M&N p218 (big R)
		p2d.w <- sum(wt.MLE*at.MLE$p2d)/sum(wt.MLE)

		if(!is.null(wt)) weights <- "User-defined"
					
	} else {
		stratified <- FALSE; #if removing empty strata leaves only one stratum
		p1hat.w <- p1hat; p2hat.w <- p2hat
		p1d.w <- p1d.MLE; p2d.w <- p2d.MLE
		wt.MLE <- NULL
		Sdot <- Q.each <- NULL
	}
	
	#fix(?) some extreme cases with zero counts
#	p2d.w[sum(x2)==0] <- 0
	
	if(stratified==TRUE && tdas==TRUE) {qtnorm <- qt(1-(1-level)/2,nstrat-1)  #for t-distribution method
	} else {qtnorm <- qnorm(1-(1-level)/2)}

	lower <- bisect(ftn=function(theta) myfun(theta)-qtnorm,contrast=contrast,dist=dist,precis=precis+1,uplow="low")
	upper <- bisect(ftn=function(theta) myfun(theta)+qtnorm,contrast=contrast,dist=dist,precis=precis+1,uplow="up")

	if(plot==TRUE) {
		if(contrast=="RD") { 
			if (dist=="bin") { xlim <-c(max(-1,min(lower-(upper-lower)/2)),min(1,max(upper+(upper-lower)/2)))
			} else if (dist=="poi") xlim <- c(lower-(upper-lower)/2,upper+(upper-lower)/2) 
		} else xlim <- c(max(0.00001,min(0.5*lower)),min(plotmax,1.5*upper))
		myseq <- seq(xlim[1],xlim[2],length.out=400)
		if(stratified) dim1 <- 1 else dim1 <- nstrat
		sc <- array(sapply(myseq,function(x) myfun(x)),dim=c(dim1,length(myseq)))
		if(stratified==FALSE) {
			qnval <- qtnorm #qnorm(1-(1-level)/2)
			ylim=c(-2.5,2.5)*qnval
			for(i in 1:nstrat) {
				plot(myseq,sc[i,],type="l",ylim=ylim,xlab=contrast,yaxs="i",ylab="Score",col="blue",main=paste0("Score function for ",ifelse(dist=="bin","binomial","Poisson")," ",contrast,"\n",x1,"/",n1," vs ",x2,"/",n2),log=ifelse(contrast=="RD","","x"))
				text(x=c(lower[i],point[i],upper[i]),y=c(-1.5,-1.75,-2)*qnval,labels=formatC(c(lower[i],point[i],upper[i]),format="fg",4,flag="#"),pos=4,offset=-0.2,cex=0.8,xpd=TRUE)
				abline(h=c(-1,1)*qnval)	
				abline(h=0,lty=2)	
				lines(rep(lower[i],2),c(ylim[1],-1.5*qnval-0.3),lty=3)		
				lines(rep(lower[i],2),c(-1.5*qnval+0.4,qnval),lty=3)		
				lines(rep(upper[i],2),c(ylim[1],-2*qnval-0.3),lty=3)	
				lines(rep(upper[i],2),c(-2*qnval+0.4,-qnval),lty=3)	
				lines(rep(point[i],2),c(ylim[1],-1.75*qnval-0.3),lty=3)		
				lines(rep(point[i],2),c(-1.75*qnval+0.4,0),lty=3)		
			}
		} else {
			qtval <- qtnorm 
			ylim=c(-2.5,2.5)*qtval
			plot(myseq,sc[1,],type="l",ylim=ylim,xlab=contrast,ylab="Score",yaxs="i",col="blue",main=paste("Score function for",dist,contrast),log=ifelse(contrast=="RD","","x"))
			abline(h=c(-1,1)*qtval)
			abline(h=0,lty=2)
			lines(rep(lower,2),c(ylim[1],-1.5*qtval-0.3),lty=3)		
			lines(rep(lower,2),c(-1.5*qtval+0.4,qtval),lty=3)		
			lines(rep(upper,2),c(ylim[1],-2*qtval-0.3),lty=3)	
			lines(rep(upper,2),c(-2*qtval+0.4,-qtval),lty=3)	
			lines(rep(point,2),c(ylim[1],-1.75*qtval-0.3),lty=3)		
			lines(rep(point,2),c(-1.75*qtval+0.4,0),lty=3)		
			text(x=c(lower,point,upper),y=c(-1.5,-1.75,-2)*qtval,labels=formatC(c(lower,point,upper),format="fg",4,flag="#"),pos=4,offset=-0.2,cex=0.8,xpd=TRUE)
		}
		if(stratified) qqnorm(Stheta.MLE)
	}
	
	#fix some extreme cases with zero counts
	if(contrast %in% c("RR","OR")) {
		if(stratified==FALSE) {	
			lower[x1==0] <- 0
			upper[x2==0] <- Inf
		} else {
			lower[sum(x1)==0] <- 0
			lower[myfun(0)-qtnorm<0] <- 0
			upper[sum(x2)==0] <- Inf			
			point[sum(x2)==0 & sum(x1)==0] <- NA			
			upper[myfun(10^100)+qtnorm>0] <- Inf
	#		upper[sum(x1)>0 & sum(x2)==0] <- Inf
		}
	}
	if(contrast %in% c("OR")) {
		if(stratified==FALSE) {	
			lower[x2==n2] <- 0
			upper[x1==n1] <- Inf
		} else {
			lower[sum(x2)==sum(n2)] <- 0
			upper[sum(x1)==sum(n1)] <- Inf			
		}
	}
	estimates <- cbind(round(cbind(Lower=lower,MLE=point,Upper=upper),precis),level=level,p1hat=p1hat.w,p2hat=p2hat.w,p1d=p1d.w,p2d=p2d.w) #,Q,tau2,het.pval)

	#optionally add p-value for a test of null hypothesis: theta<=delta
	if(contrast=="RD") delta0 <- 0 else delta0 <- 1
	if(is.null(delta)) delta <- delta0
	scorezero <- scoretheta(delta0, x1, x2, n1, n2,stratified=stratified,wt=wt,weights=weights,tdas=tdas,bcf=bcf,contrast=contrast,dist=dist,skew=skew,cc=cc)
	scoredelta <- scoretheta(delta, x1, x2, n1, n2,stratified=stratified,wt=wt,weights=weights,tdas=tdas,bcf=bcf,contrast=contrast,dist=dist,skew=skew,cc=cc)
	pval.left <- scoredelta$pval 
	pval.right <- 1-pval.left
	chisq.zero <- scorezero$score^2
	pval2sided <- pchisq(chisq.zero,1,lower.tail=FALSE)
	if(tdas==TRUE) pval2sided <- pf(chisq.zero,1,nstrat-1,lower.tail=FALSE)
	pval <- cbind(chisq=chisq.zero,pval2sided,delta=delta,scoredelta=scoredelta$score,pval.left,pval.right)
	
	outlist <- list(estimates=estimates,pval=pval,nstrat=nstrat)
	if(stratified==TRUE) {
		Qtest <- c(Q=Q.FE,tau2=tau2.FE,pval.het=pval.het,I2=I2) 
		wtpct <- 100*wt.MLE/sum(wt.MLE)
		wt1pct <- 100*wt.FE/sum(wt.FE)
		outlist <- append(outlist,list(Qtest=Qtest,weights=weights,stratdata=cbind(p1hat=p1hat,p2hat=p2hat,Q.each=Q.each,wtpct.fixed=wt1pct,wtpct.rand=wtpct)))#,p1d=p1d.MLE,p2d=p2d.MLE,Stheta=Stheta.MLE)))
	}
	outlist <- append(outlist,list(call=match.call()))	
	return(outlist)
}



#vectorized limit-finding routine - turns out not to be any quicker but is neater
#the bisection method is just as efficient as the secant method suggested by G&N, and affords greater control over whether the final estimate has score<z
#the secant method is better for RR and for Poisson rates, where there is no upper bound for d, 
#however it is not guaranteed to converge
#New version not reliant on point estimate
#This could be modified to solve upper and lower limits simultaneously
bisect <- function(ftn,contrast,dist, precis, max.iter = 100, uplow="low", logscale=FALSE){
	tiny <- (10^-(precis))/2
	nstrat <- length(eval(ftn(1)))
	hi <- rep(1,nstrat)
	dir <- lo <- rep(-1,nstrat) #rep(dir,length(proot))
	dp <- 2 #abs(proot - best)
	niter <- 1
	while(niter <= max.iter && any(dp>tiny | is.na(hi))){
		dp <- 0.5*dp
		mid <- pmax(-1,pmin(1,round((hi+lo)/2,10)))  #rounding avoids machine precision problem with, e.g. 7/10-6/10
     	if(contrast=="RD" && dist=="bin") {scor <- ftn(mid)
     	} else if((contrast=="RD" && dist=="poi") ) {scor <- ftn(round(tan(pi*mid/2),10))    #avoid machine precision producing values outside [-1,1]
     	} else if(contrast %in% c("RR","OR")) {scor <- ftn(round(tan(pi*(mid+1)/4),10))    #avoid machine precision producing values outside [-1,1]
	 	} else if(contrast=="p") scor <- ftn((mid+1)/2)   #need to modify for poisson 
     	check <- (scor<0) #| is.na(scor) #scor=NA only happens when |p1-p2|=1 and |theta|=1 (in which case hi==lo anyway), or if p1=p2=0
     	hi[check] <- mid[check]
     	lo[!check] <- mid[!check]
     	niter <- niter+1 
   }
   if(uplow=="low") best <- lo else best <- hi
   if(contrast=="RD" && dist=="bin") { return(best) 
   } else if((contrast %in% c("RD") && dist=="poi") ) { return(tan(best*pi/2)) 
   } else if(contrast %in% c("RR","OR")) {return(tan((best+1)*pi/4))	
   } else if(contrast=="p") return((best+1)/2)
}

scoretheta <- function (
	#function to evaluate the score at a given value of theta, given the observed data
	#uses the MLE solution (and notation) given in F&M, extended in Laud2015
	#This function is vectorised in x1,x2
	theta,
	x1,
	x2=NULL,
	n1,
	n2=NULL,
	dist="bin",
	contrast="RD",
	bcf=TRUE,
	skew=FALSE,
	cc=FALSE,
	stratified=FALSE,
	wt=NULL,
	weights="IV",
	tdas=FALSE,	
	...)
	
{
	nstrat <- length(n1)
	lambda <- switch(as.character(bcf), "TRUE"=(n1+n2)/(n1+n2-1), "FALSE"=1)
	p1hat <- x1/n1; p2hat <- x2/n2 
	x <- x1+x2
	N <- n1+n2

 	#RMLE of p1|theta depends on whether theta=RD or RR, and on whether a binomial or Poisson distribution is assumed
 	#Binomial RD
	if(contrast=="RD") {
		Stheta <- (p1hat-p2hat)-theta
	 	if(dist=="bin") {
			t <- n2/n1  #NB: Notation from F&M p1448, NB t='theta'=N2/N1.  Their notation for the risk difference is s_0
			s0 <- theta
			a <- N						
	   	 	b <- (n1+2*n2)*theta - N - x	
	   	 	c <- (n2*theta - N - 2*x2)*theta + x
	   	 	d <- x2*theta*(1-theta)				
		    v <- (b/a/3)^3 - b*c/(6*a*a) + d/a/2 	
    		s <- sqrt( pmax(0,(b/a/3)^2 - c/a/3)) 	#NaNs produced when? - machine precision again?
			u <- ifelse(v>0,1,-1)*s 				
			w <- (pi+acos(pmax(-1,pmin(1,ifelse(u==0 & v==0,0,v/u^3)))))/3 #avoids machine precision errors passing impossible values (outside (-1,1)) to acos
	    	p2d <- pmin(1,pmax(0,round(2*u*cos(w) - b/a/3,10))) #again, machine precision errors can give values a tiny fraction either side of 0 or 1
	    	p1d <- pmin(1,pmax(0,p2d + theta))
    		V <- (pmax(0.00000000000000,(p1d*(1-p1d)/n1 + p2d*(1-p2d)/n2)*lambda ) ) #set to zero if machine precision error in cos function produces a negative V
			mu3 <- (p1d*(1-p1d)*(1-2*p1d)/(n1^2)-p2d*(1-p2d)*(1-2*p2d)/(n2^2))
		} else if (dist=="poi") {
			A <- N
			B <- N*theta-x
			C <- -x2*theta
			num <- (-B+sqrt(B^2-4*A*C))
			p2d <- ifelse(num==0,0,num/(2*A))
			p1d  <-  p2d + theta
			V <- (p1d/n1 + p2d/n2) 
			mu3 <- (p1d/(n1^2)-p2d/(n2^2))
		}
	} else if(contrast=="RR") {
		Stheta <- p1hat-p2hat*theta
		if(dist=="bin") {
			A <- N*theta
			B <- (-(n1*theta+x1+n2+x2*theta))
			C <- x
			num <- (-B-Re(sqrt(as.complex(B^2-4*A*C))))
			p2d <- ifelse(A==0,-C/B,ifelse(num==0,0,num/(2*A)))
			p1d <- p2d*theta
			V <- pmax(0.0000000001,(p1d*(1-p1d)/n1 + (theta^2)*p2d*(1-p2d)/n2)*lambda) #A bit of a cheat to avoid a problem
			V[is.na(V)] <- Inf
			mu3 <- (p1d*(1-p1d)*(1-2*p1d)/(n1^2)-(theta^3)*p2d*(1-p2d)*(1-2*p2d)/(n2^2))
		} else if (dist=="poi") {
			p2d <- (x1+x2)/(n1*theta+n2) #incidence density version from M&N p223
			p1d  <-  p2d*theta
			V <- pmax(0,(p1d/n1 + (theta^2)*p2d/n2))
			mu3 <- (p1d/(n1^2)-(theta^3)*p2d/(n2^2))
		}
	} else if(contrast=="OR") {
		if(dist=="bin") {   #warning: singularity at theta=1 - this might not matter, except when x1=0
			A <- n2*(theta-1)
			B <- n1*theta + n2 - x*(theta-1)
			C <- -x
			num <- (-B+sqrt(B^2-4*A*C))
			p2d <- ifelse(A==0,-C/B,ifelse(num==0,0,num/(2*A)))  #If A=0 then we solve a linear equation instead
			p1d <- p2d*theta/(1+p2d*(theta-1))
			Stheta <- (p1hat-p1d)/(p1d*(1-p1d))-(p2hat-p2d)/(p2d*(1-p2d))
	    	V <- pmax(0.0000000001,(1/(n1*p1d*(1-p1d)) + 1/(n2*p2d*(1-p2d)))*lambda ) #set to zero if machine precision error in cos function produces a negative V
   			mu3 <- (1-2*p1d)/((n1*p1d*(1-p1d))^2) - (1-2*p2d)/((n2*p2d*(1-p2d))^2)
		} else if (dist=="poi") {
			print("Odds ratio not applicable to Poisson rates")
		}	
	} else if(contrast=="p") {
		Stheta <- p1hat-theta
		if(dist=="bin") {   #warning: singularity at theta=1 - this might not matter, except when x1=0
		    V <- (pmax(0.0000000001,(theta*(1-theta)/n1 ) )) #set to zero if machine precision error in cos function produces a negative V
			mu3 <- (theta*(1-theta)*(1-2*theta)/(n1^2))
		} else if (dist=="poi") {
			V <- theta/n1
			mu3 <- theta/(n1^2)
		}	
		p2d <- NA
	}
	
	#continuity corrections
	corr <- 0
	if (cc>0 && contrast=="OR") corr <- cc*(1/(n1*p1d*(1-p1d)) + 1/(n2*p2d*(1-p2d))) #Cornfield, try cc=0.25
	if (cc>0 && contrast=="RR") corr <- cc*(1/(n1) + theta/(n2))  #try 0.125 or 0.25
	if (cc>0 && contrast=="RD") 	corr <- cc*(1/pmin(n1,n2)) #Hauck Anderson - halved with cc=0.25
	if (cc>0 && contrast=="RD" && stratified==TRUE) corr <- (3/16)*(sum(n1*n2/(n1+n2)))^(-1) #from MehrotraRailkar, also Zhao et at.
	
	if (stratified==TRUE) {
		pval <- NA
		if(is.null(wt)) {
			if(weights=="MH") { wt <- n1*n2/(n1+n2) #MH weights for RD, applied across other comparative parameters too (without theoretical justification for OR)
			} else if(weights=="IV") {if (all(V==0) || all(is.na(V))) wt <- rep(1,nstrat) else wt <- 1/V  #IVS: inverse variance weights updated wih V.tilde
			} else if(weights=="MN") {
				if(contrast=="RR") { #M&Ns iterative weights - quite similar to MH
					wtx <- (1/n1 + theta/n2)^(-1)
					p2ds <- sum(wtx*p2d)/sum(wtx)
					p1ds <- sum(wtx*p1d)/sum(wtx)
					wty <- ((1-p1ds)/(n1*(1-p2ds)) + theta/n2)^(-1)
					wty[p2ds==1] <- 0
					p2ds <- sum(wty*p2d)/sum(wty)
					p1ds <- sum(wty*p1d)/sum(wty)
					wt <- ((1-p1ds)/(n1*(1-p2ds)) + theta/n2)^(-1)
					wt[p2ds==1] <- 0
				} else if(contrast=="RD") { #M&Ns iterative weights - quite similar to MH
					#wtx <- n1*n2/(n1+n2)
					wtx <- (1/n1 + 1/n2)^(-1)
					p2ds <- sum(wtx*p2d)/sum(wtx)
					p1ds <- sum(wtx*p1d)/sum(wtx)
					wty <- ((p1ds*(1-p1ds)/(p2ds*(1-p2ds)))/n1 + 1/n2)^(-1)
					wty[p2ds==0] <- 0
					p2ds <- sum(wty*p2d)/sum(wty)
					p1ds <- sum(wty*p1d)/sum(wty)
					wt <- ((p1ds*(1-p1ds)/(p2ds*(1-p2ds)))/n1 + 1/n2)^(-1)
					#wt <- (1/n1 + (p2ds*(1-p2ds))/n2)^(-1) #equivalent, according to ...
					wt[p2ds==0] <- 0
				} else if(contrast=="OR") { #M&Ns weights are very similar in structure to IVS
					 wt <- n1*n2*((1-p1d)*p2d)^2 / (n1*p1d*(1-p1d) + n2*p2d*(1-p2d)) 
				}
			}
		} else weights <- "User-defined"

		Sdot <- sum(wt*Stheta)/sum(wt)
		if(weights=="IV") {Q.i <- wt*(( Stheta - Sdot )^2)  #This version for iterative weights?
		} else Q.i <- (( Stheta - Sdot )^2)/V
		Q <- sum(Q.i)
		W <- sum(wt)

		if(weights=="IV") { tau2 <- max(0,(Q-(nstrat-1)) / (W - (sum(wt^2)/W)) )  # published formula for IVS weights
		} else tau2 <- max(0,(Q-(nstrat-1*(sum(V*wt^2)/W))) / (sum(1/V)-1*(sum(wt^2)/W)) )  #only needed if want to output tau2. uncertain derivation - seems to be experimenting with wrat=1

		if(tdas==TRUE && weights=="IV" &&  !(all(V==0) || all(is.na(V))) ) { wt <- 1/(V+tau2) } #updated inverse variance weights, as per meta package - makes SCr less conservative? (CHECK) and very similar to HK

		Sdot <- sum(wt*Stheta)/sum(wt)
		VS <- sum( wt*(Stheta - Sdot)^2 )/((nstrat-1)*sum(wt))

		if(tdas==TRUE) t2 <- tau2 else t2 <- 0
		Vdot <- sum(((wt/sum(wt))^2)*(V + t2)) #from equation (15) of Miettinen&Nurminen, with the addition of between strata variance from Whitehead&Whitehead

		score1 <-  sum((wt/sum(wt))*(Stheta-corr))/pmax(0.00000001,sqrt(Vdot))
		scterm <- sum(((wt/sum(wt))^3)*mu3)/(6*Vdot^(3/2))
		A <- scterm
		B <- 1
		C <- -(score1 + scterm)
		num <- (-B+ Re(sqrt(as.complex(B^2-4*A*C))))
		score <- ifelse((skew==FALSE | scterm==0),score1,num/(2*A))
		pval <- pnorm(score)
		if(tdas==TRUE) {
			score <- sum((wt/sum(wt))*(Stheta-corr)) / sqrt(VS) #simplified version that ignores tau2 altogether for MH weights but seems to work. Is consistent with Sidik2006
			pval <- pt(score,nstrat-1)
		}	
		p2ds <- sum(wt*p2d/sum(wt))
		p1ds <- sum(wt*p1d/sum(wt))
	} 
	else if (stratified==FALSE) {
 		corr <- corr*sign(Stheta) 
		p1ds <- p1d
		p2ds <- p2d
		
		#Calculation of score & p-value involves solving z.p = Stheta/sqrt(V) - (z.p^2)*mu3/(6*V^(3/2)) + mu3/(6*V^(3/2))
		#Note that in the special case of mu3=0, this reduces to the skew=FALSE case.
		A <- mu3/(6*V^(3/2))
		B <- 1
		C <- -((Stheta-corr)/sqrt(V) + mu3/(6*V^(3/2)))
		num <- (-B + Re(sqrt(as.complex(B^2-4*A*C))))
		score <- ifelse((skew==FALSE | mu3==0),ifelse(Stheta==0,0,(Stheta-corr)/sqrt(V)),num/(2*A))
		pval <- pnorm(score)
	}

	outlist <- list(score=score,p1d=p1d,Stheta=Stheta,V=V,p2d=p2d,mu3=mu3,pval=pval)
	if(stratified) {outlist <- append(outlist,list(Sdot=Sdot,Vdot=Vdot,tau2=tau2,VS=VS,t2=t2,Q.i=Q.i,Q=Q,wt=wt,p1ds=p1ds,p2ds=p2ds))}
	return(outlist)
}


if(FALSE) {
#Hartung & Knapp stratified example:
x1 <- c(15,12,29,42,14,44,14,29,10,17,38,19,21)
x2 <- c(9,1,18,31,6,17,7,23,3,6,12,22,19)
n1 <- c(16,16,34,56,22,54,17,58,14,26,44,29,38)
n2 <- c(16,16,34,56,22,55,15,58,15,27,45,30,38)

ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=F,skew=TRUE,cc=F) $estimates[c(2,12),]

fround <- function(x,digits=6) paste(format(round(x[1],digits=digits),nsmall=digits)," (",paste(format(round(x[2:3],digits=digits),nsmall=digits),collapse=", "),")",sep="")
mytab <-rbind(
SCASmh=c(
RDbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="MH",skew=T,tdas=F)$estimates[,c(2,1,3)],3),
RRbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="MH",skew=T,tdas=F)$estimates[,c(2,1,3)],2),
OR=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="MH",skew=T,tdas=F)$estimates[,c(2,1,3)],2)
),
SCASiv=c(
RDbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="IV",skew=T,tdas=F)$estimates[,c(2,1,3)],3),
RRbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="IV",skew=T,tdas=F)$estimates[,c(2,1,3)],2),
OR=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="IV",skew=T,tdas=F)$estimates[,c(2,1,3)],2)
),
TDAS=c(
RDbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RD",stratified=T,weights="IV",tdas=T)$estimates[,c(2,1,3)],3),
RRbin=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="RR",stratified=T,weights="IV",tdas=T)$estimates[,c(2,1,3)],2),
OR=fround(ratesCI(x1=x1,x2=x2,n1=n1,n2=n2,contrast="OR",stratified=T,weights="IV",tdas=T)$estimates[,c(2,1,3)],2)
)
)
mytab

}
