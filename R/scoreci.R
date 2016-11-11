#' Score confidence intervals for comparisons of independent binomial or Poisson
#' rates.
#' 
#' Score-based confidence intervals for the rate (or risk) difference ("RD") or 
#' ratio ("RR") for independent binomial or Poisson rates, or for odds ratio 
#' ("OR", binomial only). Including options for bias correction (from Miettinen 
#' & Nurminen), skewness correction ("GNbc" method from Laud & Dane, developed
#' from Gart & Nam, and generalised as "SCAS" in forthcoming publication) and
#' continuity correction. Also includes intervals for a single proportion, i.e.
#' Wilson score method, with skewness correction, which has slightly better
#' coverage properties than the Jeffreys method. This function is vectorised in
#' x1, x2, n1, and n2.  Vector inputs may also be combined into a single
#' stratified analysis (e.g. meta-analysis), either using fixed effects, or the
#' more general "TDAS" method, which incorporates stratum variability using a
#' t-distribution score (inspired by to Hartung-Knapp-Sidik-Jonkman).
#' 
#' @param x1,x2 Numeric vectors of numbers of events in group 1 & group 2 
#'   respectively.
#' @param n1,n2 Numeric vectors of sample sizes (for binomial rates) or exposure
#'   times (for Poisson rates) in each group.
#' @param distrib Character string indicating distribution assumed for the input
#'   data: "bin" = binomial (default), "poi" = Poisson.
#' @param contrast Character string indicating the contrast of interest: "RD" = 
#'   rate difference (default), "RR" = rate ratio, "OR" = odds ratio, "p" =
#'   single proportion.
#' @param level Number specifying confidence level (between 0 and 1, default 
#'   0.95).
#' @param skew Logical (default TRUE) indicating whether to apply skewness
#'   correction (Laud 2016).
#' @param bcf Logical (default TRUE) indicating whether to apply bias correction in the score 
#'   denominator. Applicable to distrib = "bin" only. (NB: bcf = FALSE option is
#'   really only included for legacy validation against previous published 
#'   methods (i.e. Gart & Nam, Mee).
#' @param cc Number or logical (default FALSE) specifying (amount of) continuity
#'   correction.
#' @param delta Number to be used in a one-sided significance test (e.g. 
#'   non-inferiority margin). 1-sided p-value will be <0.025 iff 2-sided 95\% CI
#'   excludes delta. NB: can also be used for a superiority test by setting 
#'   delta=0.
#' @param precis Number (default 6) specifying precision to be used in
#'   optimisation subroutine (i.e. number of decimal places).
#' @param plot Logical (default FALSE) indicating whether to output plot of the
#'   score function
#' @param plotmax Numeric value indicating maximum value to be displayed on 
#'   x-axis of plots (useful for ratio contrasts which can be infinite).
#' @param stratified Logical (default FALSE) indicating whether to combine
#'   vector inputs into a single stratified analysis.
#' @param weighting String indicating which weighting method to use if 
#'   stratified = "TRUE":  "IVS" = Inverse Variance of Score (default), "MH" = 
#'   Mantel-Haenszel, "MN" = Miettinen-Nurminen iterative weights.
#' @param wt Numeric vector containing (optional) user-specified weights.
#' @param tdas Logical (default FALSE) indicating whether to use t-distribution
#'   method for stratified data (defined in Laud 2016).
#' @param ... Other arguments.
#' @importFrom stats pchisq pf pnorm pt qbeta qgamma qnorm qqnorm qt
#' @importFrom graphics abline lines text
#' @return A list containing the following components: \describe{ 
#'   \item{estimates}{a matrix containing estimates of the rates in each group 
#'   and of the requested contrast, with its confidence interval} \item{pval}{a 
#'   matrix containing details of the corresponding 2-sided significance test 
#'   against the null hypothesis that p_1 = p_2, and one-sided significance 
#'   tests agains the null hypothesis that theta >= or <= delta} 
#'   \item{call}{details of the function call} }
#'   If stratified = TRUE, the following outputs are added: \describe{
#'   \item{Qtest}{a vector of values descibing and testing heterogeneity}
#'   \item{weighting}{a string indicating the selected weighting method}
#'   \item{stratdata}{a matrix containing stratum estimates and weights}}
#' @examples  
#'   #Binomial RD, SCAS method:
#'   scoreci(x1 = c(12,19,5), n1 = c(16,29,56), x2 = c(1,22,0), n2 = c(16,30,29))
#'   
#'   #Binomial RD, MN method:
#'   scoreci(x1 = c(12,19,5), n1 = c(16,29,56), x2 = c(1,22,0), n2 = c(16,30,29), skew = FALSE)
#'   
#'   #Poisson RR, SCAS method:
#'   scoreci(x1 = 5, n1 = 56, x2 = 0, n2 = 29, distrib = "poi", contrast = "RR")
#'   
#'   #Poisson RR, MN method:
#'   scoreci(x1 = 5, n1 = 56, x2 = 0, n2 = 29, distrib = "poi", contrast = "RR", skew = FALSE)
#'   
#'   #Binomial rate, SCAS method:
#'   scoreci(x1 = c(5,0), n1 = c(56,29), contrast = "p")
#'   
#'   #Binomial rate, Wilson score method:
#'   scoreci(x1 = c(5,0), n1 = c(56,29), contrast = "p", skew = FALSE)
#'   
#'   #Poisson rate, SCAS method:
#'   scoreci(x1 = c(5,0), n1 = c(56,29), distrib = "poi", contrast = "p")
#'
#'   #Stratified example, using data from Hartung & Knapp:
#'   scoreci(x1 = c(15,12,29,42,14,44,14,29,10,17,38,19,21),
#'           x2 = c(9,1,18,31,6,17,7,23,3,6,12,22,19),
#'           n1 = c(16,16,34,56,22,54,17,58,14,26,44,29,38),
#'           n2 = c(16,16,34,56,22,55,15,58,15,27,45,30,38),
#'           stratified = TRUE)
#'           
#'   #TDAS example, using data from Hartung & Knapp:
#'   scoreci(x1 = c(15,12,29,42,14,44,14,29,10,17,38,19,21),
#'           x2 = c(9,1,18,31,6,17,7,23,3,6,12,22,19),
#'           n1 = c(16,16,34,56,22,54,17,58,14,26,44,29,38),
#'           n2 = c(16,16,34,56,22,55,15,58,15,27,45,30,38),
#'           stratified = TRUE, tdas = TRUE)
#'
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#' @references 
#'   Laud PJ. Equal-tailed confidence intervals for comparison of 
#'   rates: Submitted to Pharmaceutical Statistics for peer review.
#'   
#'   Laud PJ, Dane A. Confidence intervals for the difference between independent 
#'   binomial proportions: comparison using a graphical approach and
#'   moving averages. Pharmaceutical Statistics 2014; 13(5):294–308
#'   
#'   Miettinen OS, Nurminen M. Comparative analysis of two rates. Statistics in 
#'   Medicine 1985; 4:213-226.
#'   
#'   Farrington CP, Manning G. Test statistics and sample size formulae for 
#'   comparative binomial trials with null hypothesis of non-zero risk 
#'   difference or non-unity relative risk. Statistics in Medicine 1990; 
#'   9(12):1447-1454.
#'   
#'   Gart JJ, Nam Jm. Approximate interval estimation of the ratio of binomial 
#'   parameters: A review and corrections for skewness. Biometrics 1988; 
#'   44(2):323-338.
#'   
#'   Gart JJ, Nam Jm. Approximate interval estimation of the difference in 
#'   binomial parameters: correction for skewness and extension to multiple 
#'   tables. Biometrics 1990; 46(3):637-643.
#' @export
scoreci <- function(
	x1,
	n1,
	x2 = NULL,
	n2 = NULL,
	distrib = "bin",
	contrast = "RD",
	level = 0.95,
	skew = TRUE,
	bcf = TRUE,
	cc = 0,
	delta = NULL,
	precis = 6,
	plot = FALSE,	
	plotmax = 100,
	stratified = FALSE,
	weighting = "IVS",
	wt = NULL,
	tdas = FALSE,	
	...
	) { 
  if (!(tolower(substr(distrib, 1, 3)) %in% c("bin", "poi"))) {
    print("Distrib must be one of 'bin' or 'poi'")
    stop()
  }
  if (!(tolower(substr(contrast, 1, 2)) %in% c("rd", "rr", "or", "p"))) {
    print("Contrast must be one of 'RD', 'RR', 'OR' or 'p'")
    stop()
  }   
  if (contrast != "p" && (is.null(x2) || is.null(n2))) {
	  print("argument x2 or n2 missing")
	  stop()
	}
  if (!is.numeric(c(x1, n1, x2, n2, delta))) {
		print("Non-numeric inputs!")
		stop()
	}
	if (any(c(x1, n1, x2, n2) < 0)) {
		print("Negative inputs!")
		stop()
	}	
  if (distrib == "bin" && (any(x1 > n1 + 0.001) || any(x2 > n2 + 0.001))) {
    print("x1 > n1 or x2 > n2 not possible for distrib = 'bin'")
    stop()
  }
  if (!is.null(delta)) {	
		if (contrast == "RD") {
			if (distrib == "bin" && (delta < -1 || delta > 1)) {
				print("Impossible delta!")
				stop()			
			}
		} else if (contrast == "p") {
		  if (delta < 0 || (distrib == "bin" && delta > 1)) {
		    print("Impossible delta!")
		    stop()			
		  }
		}
	}
  if (distrib == "poi" && contrast == "OR") {
    print("Odds ratio not applicable to Poisson rates")
    stop()
  }
  if(as.character(cc) == "TRUE") cc <- 0.5

	nstrat <- length(x1)
	# in case x1,x2 are input as vectors but n1,n2 are not
	if (length(n1) < nstrat && nstrat > 1) n1 <- rep(n1, length.out = nstrat)
	if (contrast != "p" && length(n2) < nstrat && nstrat > 1) {
	  n2 <- rep(n2, length.out = nstrat)
	}
	
	#check for empty strata, and for x1<=n1, x2<=n2
	if (stratified) {
	  # for stratified calculations, remove any strata with no observations
	  empty.strat <- ((n1 == 0 | n2 == 0) |
	                    (contrast %in% c("RR", "OR") & (x1 == 0 & x2 == 0)) |
	                    (contrast == "OR" & (x1 == n1 & x2 == n2)))
	  x1 <- x1[!empty.strat]
	  x2 <- x2[!empty.strat]
	  n1 <- n1[!empty.strat]
	  n2 <- n2[!empty.strat]
	  if(sum(empty.strat)>0) {
	    print('Warning: at least one stratum contributed no information and was removed')
	  }
	  
	  # for double-zero cells for RD, add 0.5 to avoid 100% weight with IVS weights.
	  # Same for double-100% cells for RR & RD. NB this should only be necessary
	  # for weighting='IVS'
	  if (weighting == "IVS") {
	    zeroRD <- (contrast %in% c("RD") & (x1 == 0 & x2 == 0))
	    x1 <- x1 + (zeroRD * 0.5)
	    x2 <- x2 + (zeroRD * 0.5)
	    n1 <- n1 + (zeroRD * 1)
	    n2 <- n2 + (zeroRD * 1)
	    fullRD <- (contrast %in% c("RD", "RR") & (x1 == n1 & x2 == n2))
	    x1 <- x1 + (fullRD * 0.5)
	    x2 <- x2 + (fullRD * 0.5)
	    n1 <- n1 + (fullRD * 1)
	    n2 <- n2 + (fullRD * 1)
	  }
	}	
	nstrat <- length(x1) # update nstrat after removing empty strata

	#Warnings if removal of empty strata leave 0 or 1 stratum
	if (stratified == TRUE && nstrat <= 1) {
	  stratified <- FALSE
	  tdas = FALSE
	  print("Warning: only one stratum!")
	}
	if (nstrat == 0) {
	  print("Warning: no usable data!")
	  if (contrast %in% c("RR", "OR")) {
	    x1 <- x2 <- 0
	    n1 <- n2 <- 10
	  }
	}

	p1hat <- x1/n1
	p2hat <- x2/n2
	
	#wrapper function for scoretheta
	myfun <- function(theta, randswitch = tdas, ccswitch = cc) {
	  scoretheta(theta = theta, x1 = x1, x2 = x2, n1 = n1, n2 = n2, bcf = bcf,
	             contrast = contrast, distrib = distrib, stratified = stratified,
	             wt = wt, weighting = weighting, tdas = randswitch, skew = skew,
	             cc = ccswitch)$score  
	}

	#find point estimate for theta as the value where scoretheta = 0
	#fixed effects point estimate taken with no cc
	point.FE <- bisect(ftn = function(theta) 
	  myfun(theta, randswitch = FALSE, ccswitch = 0) - 0, contrast = contrast,
	  distrib = distrib, precis = precis + 1, uplow = "low")
	point <- point.FE
	  
	# random effects point estimate if required
	if (stratified == TRUE) {
	  point <- bisect(ftn = function(theta) 
	    myfun(theta, randswitch = tdas, ccswitch = 0) - 0, contrast = contrast, 
	    distrib = distrib, precis = precis + 1, uplow = "low")  
	}
	
	# fix some extreme cases with zero counts
	if (stratified) {
	  if (skew == FALSE & contrast %in% c("RR", "OR")) 
	    point[sum(x2) > 0 & sum(x1) == 0] <- 0
	  if (skew == FALSE & contrast %in% c("RR", "OR")) 
	    point[sum(x1) > 0 & sum(x2) == 0] <- Inf
	  # if(contrast %in% c('OR')) point[sum(x1)==sum(n1) & sum(x2)>0] <- Inf
	} else {
	  if (contrast %in% c("RR", "OR")) {
	    ## & skew==FALSE)
	    point[x1 == 0 & x2 == 0] <- NA
	    point[(x1 > 0 & x2 == 0) && skew == FALSE] <- Inf
	    if (distrib == "bin") point[x1 == n1 & x2 == n2] <- 1
	  }
	  if (contrast == "OR") {
	    point[x1 == n1 & x2 == n2] <- NA
	    point[(x1 == n1 & x2 > 0) && skew == FALSE] <- Inf
	  }
	  if (contrast == "RD") {
	    if (distrib == "bin") point[x1 == n1 & x2 == n2] <- 0
	    point[x1 == 0 & x2 == 0] <- 0
	  }
	}

	# identify certain quantities evaluated at the maximum likelihood point
	# estimate for theta (incorporating random effects)
	at.MLE <- scoretheta(theta = point, x1 = x1, x2 = x2, n1 = n1, n2 = n2, bcf = bcf,
	                     contrast = contrast,
	                     distrib = distrib, stratified = stratified,
	                     weighting = weighting, wt = wt, tdas = tdas, skew = skew,
	                     cc = cc)
	Stheta.MLE <- at.MLE$Stheta
	p1d.MLE <- at.MLE$p1d
	p2d.MLE <- at.MLE$p2d
	wt.MLE <- at.MLE$wt
	V.MLE <- at.MLE$V

	# if stratified=TRUE, options are available for assuming fixed effects (tdas=FALSE)
	# or random effects (tdas=T). The IVS weights are different for each version, which 
	# in turn can lead to a different point estimate, at which certain quantities are 
	# evaluated. However, fixed
	# effects estimates are needed for both, in particular for the heterogeneity test
	if (stratified == TRUE && nstrat > 1) {
	  at.FE <- scoretheta(theta = point.FE, x1 = x1, x2 = x2, n1 = n1, n2 = n2,
	                      bcf = bcf, contrast = contrast,
	                      distrib = distrib, stratified = stratified,
	                      weighting = weighting, wt = wt, tdas = FALSE,
	                      skew = skew, cc = cc)
	  wt.FE <- at.FE$wt
	  tau2.FE <- at.FE$tau2
	  Q.each <- at.FE$Q.i
	  Q.FE <- at.FE$Q
	  I2 <- max(0, 100 * (Q.FE - (nstrat - 1))/Q.FE)
	  pval.het <- 1 - pchisq(Q.FE, nstrat - 1)
	  # Qualitative interaction test to be added...
	  
	  # as per M&N p218 (little r), actually no longer needed for point estimate. 
	  p1hat.w <- sum(wt.MLE * x1/n1)/sum(wt.MLE)  
	  p2hat.w <- sum(wt.MLE * x2/n2)/sum(wt.MLE)
	  # as per M&N p218 (big R)
	  p1d.w <- sum(wt.MLE * at.MLE$p1d)/sum(wt.MLE) 
	  if (contrast != "p") {
	    p2d.w <- sum(wt.MLE * at.MLE$p2d)/sum(wt.MLE)
	  } else p2d.w <- NULL
	  
	  if (!is.null(wt)) weighting <- "User-defined"
	  
	} else {
	  # if removing empty strata leaves only one stratum treat as unstratified
	  #stratified <- FALSE #this is already done above
	  p1hat.w <- p1hat
	  p2hat.w <- p2hat
	  p1d.w <- p1d.MLE
	  if (contrast != "p") p2d.w <- p2d.MLE else p2d.w <- NULL
	  wt.MLE <- NULL
	  Sdot <- Q.each <- NULL
	}
	
	#z- or t- quantile required for specified significance level
	if (stratified == TRUE && tdas == TRUE) {
	  qtnorm <- qt(1 - (1 - level)/2, nstrat - 1)  #for t-distribution method
	} else {
	  qtnorm <- qnorm(1 - (1 - level)/2)
	}
	
	# Use bisection routine to locate lower and upper confidence limits
	lower <- bisect(ftn = function(theta)
	  myfun(theta) - qtnorm, contrast = contrast, distrib = distrib,
	  precis = precis + 1, uplow = "low")
	upper <- bisect(ftn = function(theta)
	  myfun(theta) + qtnorm, contrast = contrast, distrib = distrib,
	  precis = precis + 1, uplow = "up")

  # Optional plot of the score function.
	# Ideally this would be in a separate function, but it is unlikely to be used
	# much in practice - only included for code development and validation purposes.
	if (plot == TRUE) {
	  if (contrast == "RD") {
	    if (distrib == "bin") {
	      xlim <- c(max(-1, min(lower - (upper - lower)/2)),
	                min(1, max(upper + (upper - lower)/2)))
	    } else if (distrib == "poi") {
	      xlim <- c(lower - (upper - lower)/2, upper + (upper - lower)/2)
	    }
	  } else xlim <- c(max(0, min(0.5 * lower)), min(plotmax, max(1.5 * upper)))
	myseq <- seq(xlim[1], xlim[2], length.out = 400)
	  if (stratified) dim1 <- 1 else dim1 <- nstrat
	  sc <- array(sapply(myseq, function(x) myfun(x)), dim = c(dim1, length(myseq)))
	  if (stratified == FALSE) {
	    qnval <- qtnorm
	    ylim = c(-2.5, 2.5) * qnval
	    for (i in 1:nstrat) {
	      plot(myseq, sc[i, ], type = "l", xlim = xlim, ylim = ylim, xlab = contrast,
	           yaxs = "i", ylab = "Score", col = "blue", 
	           main = paste0("Score function for ",
	                         ifelse(distrib == "bin", "binomial", "Poisson"), " ",
	                         contrast, "\n", x1[i], "/", n1[i], 
	                         ifelse(contrast == "p", "", paste0(" vs ", x2[i], "/", n2[i])))
	           #log = ifelse(contrast == "RD", "", "x")
	           )
	      text(x = c(lower[i], point[i], upper[i]), y = c(-1.5, -1.75, -2) * qnval,
	           labels = formatC(c(lower[i], point[i], upper[i]), format = "fg", 4,
	                            flag = "#"),
	           pos = 4, offset = -0.2, cex = 0.8, xpd = TRUE)
	      abline(h = c(-1, 1) * qnval)
	      abline(h = 0, lty = 2)
	      lines(rep(lower[i], 2), c(ylim[1], -1.5 * qnval - 0.3), lty = 3)
	      lines(rep(lower[i], 2), c(-1.5 * qnval + 0.4, qnval), lty = 3)
	      lines(rep(upper[i], 2), c(ylim[1], -2 * qnval - 0.3), lty = 3)
	      lines(rep(upper[i], 2), c(-2 * qnval + 0.4, -qnval), lty = 3)
	      lines(rep(point[i], 2), c(ylim[1], -1.75 * qnval - 0.3), lty = 3)
	      lines(rep(point[i], 2), c(-1.75 * qnval + 0.4, 0), lty = 3)
	    }
	  } else {
	    qtval <- qtnorm
	    ylim = c(-2.5, 2.5) * qtval
	    plot(myseq, sc[1, ], type = "l", ylim = ylim, xlab = contrast,
	         ylab = "Score", yaxs = "i", col = "blue",
	         main = paste("Score function for", distrib, contrast)
	         #log = ifelse(contrast == "RD", "", "x")
	         )
	    abline(h = c(-1, 1) * qtval)
	    abline(h = 0, lty = 2)
	    lines(rep(lower, 2), c(ylim[1], -1.5 * qtval - 0.3), lty = 3)
	    lines(rep(lower, 2), c(-1.5 * qtval + 0.4, qtval), lty = 3)
	    lines(rep(upper, 2), c(ylim[1], -2 * qtval - 0.3), lty = 3)
	    lines(rep(upper, 2), c(-2 * qtval + 0.4, -qtval), lty = 3)
	    lines(rep(point, 2), c(ylim[1], -1.75 * qtval - 0.3), lty = 3)
	    lines(rep(point, 2), c(-1.75 * qtval + 0.4, 0), lty = 3)
	    text(x = c(lower, point, upper), y = c(-1.5, -1.75, -2) * qtval,
	         labels = formatC(c(lower, point, upper), format = "fg", 4, flag = "#"),
	         pos = 4, offset = -0.2, cex = 0.8, xpd = TRUE)
	  }
	  if (stratified) 
	    qqnorm(Stheta.MLE)
	}

	# fix some extreme cases with zero counts
	#if (contrast == "RR" && skew == FALSE) p2d.w[sum(x2) == 0] <- 0
	
	if (contrast %in% c("RR", "OR")) {
	  if (stratified == FALSE) {
	    lower[x1 == 0] <- 0
	    upper[x2 == 0] <- Inf
	  } else {
	    lower[sum(x1) == 0] <- 0
	    lower[myfun(0) - qtnorm < 0] <- 0
	    upper[sum(x2) == 0] <- Inf
	    point[sum(x2) == 0 & sum(x1) == 0] <- NA
	    upper[myfun(10^100) + qtnorm > 0] <- Inf
	    # upper[sum(x1)>0 & sum(x2)==0] <- Inf
	  }
	}
	if (contrast %in% c("OR")) {
	  if (stratified == FALSE) {
	    lower[x2 == n2] <- 0
	    upper[x1 == n1] <- Inf
	  } else {
	    lower[sum(x2) == sum(n2)] <- 0
	    upper[sum(x1) == sum(n1)] <- Inf
	  }
	}
	if (stratified == FALSE) {
  	inputs <- cbind(x1 = x1, n1 = n1)
	  if (contrast != "p") inputs <- cbind(inputs, x2 = x2, n2 = n2)
	} else inputs <- NULL
	
	estimates <- cbind(
	  round(cbind(Lower = lower, MLE = point, Upper = upper), precis),
	  level = level, inputs,
	  p1hat = p1hat.w, p2hat = p2hat.w, p1mle = p1d.w, p2mle = p2d.w)

	# optionally add p-value for a test of null hypothesis: theta<=delta
	# default value of delta depends on contrast
	if (contrast == "RD") {
	  delta0 <- 0 
	} else if (contrast == "p") {
	  delta0 <- 0.5
	} else delta0 <- 1
	if (is.null(delta)) delta <- delta0
	scorezero <- scoretheta(theta = delta0, x1 = x1, x2 = x2, n1 = n1, n2 = n2,
	                        stratified = stratified,
	                        wt = wt, weighting = weighting, tdas = tdas,
	                        bcf = bcf, contrast = contrast, distrib = distrib,
	                        skew = skew, cc = cc)
	scoredelta <- scoretheta(theta = delta, x1 = x1, x2 = x2, n1 = n1, n2 = n2,
	                         stratified = stratified,
	                         wt = wt, weighting = weighting, tdas = tdas,
	                         bcf = bcf, contrast = contrast, distrib = distrib,
	                         skew = skew, cc = cc)
	pval.left <- scoredelta$pval
	pval.right <- 1 - pval.left
	chisq.zero <- scorezero$score^2
	pval2sided <- pchisq(chisq.zero, 1, lower.tail = FALSE)
	if (tdas == TRUE) pval2sided <- pf(chisq.zero, 1, nstrat - 1, lower.tail = FALSE)
	pval <- cbind(chisq = chisq.zero, pval2sided, delta = delta,
	              scoredelta = scoredelta$score, pval.left, pval.right)
	
	outlist <- list(estimates = estimates, pval = pval) 
	if (stratified == TRUE) {
	  Qtest <- c(Q = Q.FE, tau2 = tau2.FE, pval.het = pval.het, I2 = I2)
	  wtpct <- 100 * wt.MLE/sum(wt.MLE)
	  wt1pct <- 100 * wt.FE/sum(wt.FE)
	  outlist <- append(outlist,
	    list(Qtest = Qtest, weighting = weighting, 
	    stratdata = cbind(x1j = x1, n1j = n1, x2j = x2, n2j = n2,
	                      p1hatj = p1hat, p2hatj = p2hat, Qj = Q.each,
	                      wtpct.fixed = wt1pct, wtpct.rand = wtpct))) 
#	    p1d=p1d.MLE,p2d=p2d.MLE,Stheta=Stheta.MLE,V.MLE)))
	}
	outlist <- append(outlist, list(call = c(distrib = distrib,
	                 contrast = contrast, level = level, skew = skew,
	                 bcf = bcf, cc = cc)))
	return(outlist)
}

#' Skewness-corrected asymptotic score ("SCAS") confidence intervals for
#' comparisons of independent binomial or Poisson rates.
#' 
#' Wrapper function for the SCAS method. Score-based confidence intervals for
#' the rate (or risk) difference ("RD") or ratio ("RR") for independent binomial
#' or Poisson rates, or for odds ratio ("OR", binomial only), or the single rate
#' ("p"). (This is the "GNbc" method from Laud & Dane, developed from Gart &
#' Nam, and generalised as "SCAS" in forthcoming publication) including optional
#' continuity correction.  This function is vectorised in x1, x2, n1, and n2.
#' Vector inputs may also be combined into a single stratified analysis (e.g.
#' meta-analysis). This method assumes the contrast is constant across strata
#' (fixed effects).  For a 'random-effects' method use tdasci (or scoreci with
#' tdas = TRUE).
#' 
#' @param x1,x2 Numeric vectors of numbers of events in group 1 & group 2 
#'   respectively.
#' @param n1,n2 Numeric vectors of sample sizes (for binomial rates) or exposure
#'   times (for Poisson rates) in each group.
#' @inheritParams scoreci
#' @export
scasci <- function(
  x1,
  n1,
  x2 = NULL,
  n2 = NULL,
  distrib = "bin",
  contrast = "RD",
  level = 0.95,
  cc = 0,
  delta = NULL,
  precis = 6,
  plot = FALSE,	
  plotmax = 100,
  stratified = FALSE,
  weighting = "IVS",
  wt = NULL,
  ...
) { 
  scoreci(
    x1 = x1,
    n1 = n1,
    x2 = x2,
    n2 = n2,
    distrib = distrib,
    contrast = contrast,
    level = level,
    cc = cc,
    delta = delta,
    precis = precis,
    plot = plot,	
    plotmax = plotmax,
    stratified = stratified,
    weighting = weighting,
    wt = wt,
    skew = TRUE,
    bcf = TRUE,
    ...
  ) 
}

#' t-distribution asymptotic score ("TDAS") confidence intervals for
#' comparisons of independent binomial or Poisson rates.
#' 
#' Wrapper function for the TDAS method. Score-based stratified confidence
#' intervals for the rate (or risk) difference ("RD") or ratio ("RR") for
#' independent binomial or Poisson rates, or for odds ratio ("OR", binomial
#' only), or the single rate ("p"). This function combines vector inputs into a
#' single stratified analysis (e.g. meta-analysis). The TDAS method incorporates
#' any stratum variability into the confidence interval.
#' 
#' @param x1,x2 Numeric vectors of numbers of events in group 1 & group 2 
#'   respectively.
#' @param n1,n2 Numeric vectors of sample sizes (for binomial rates) or exposure
#'   times (for Poisson rates) in each group.
#' @inheritParams scoreci
#' @export
tdasci <- function(
  x1,
  n1,
  x2 = NULL,
  n2 = NULL,
  distrib = "bin",
  contrast = "RD",
  level = 0.95,
  cc = 0,
  delta = NULL,
  precis = 6,
  plot = FALSE,	
  plotmax = 100,
  weighting = "IVS",
  wt = NULL,
  ...
) { 
  scoreci(
    x1 = x1,
    n1 = n1,
    x2 = x2,
    n2 = n2,
    distrib = distrib,
    contrast = contrast,
    level = level,
    cc = cc,
    delta = delta,
    precis = precis,
    plot = plot,	
    plotmax = plotmax,
    stratified = TRUE,
    weighting = weighting,
    wt = wt,
    skew = FALSE,
    bcf = TRUE,
    ...
  ) 
}

# vectorized limit-finding routine - turns out not to be any quicker but is
# neater. The bisection method is just as efficient as the secant method
# suggested by G&N, and affords greater control over whether the final estimate
# has score<z the secant method is better for RR and for Poisson rates, where
# there is no upper bound for d, however it is not guaranteed to converge New
# version not reliant on point estimate This could be modified to solve upper
# and lower limits simultaneously
# Internal function
bisect <- function(
  ftn,
  contrast,
  distrib,
  precis,
  max.iter = 100,
  uplow = "low"
  ) {
  tiny <- (10^-(precis))/2
  nstrat <- length(eval(ftn(1)))
  hi <- rep(1, nstrat)
  dir <- lo <- rep(-1, nstrat)  
  dp <- 2  
  niter <- 1
  while (niter <= max.iter && any(dp > tiny | is.na(hi))) {
    dp <- 0.5 * dp
    mid <- pmax(-1, pmin(1, round((hi + lo)/2, 10)))  
    # rounding avoids machine precision problem with, e.g. 7/10-6/10
    if (contrast == "RD" && distrib == "bin") {
      scor <- ftn(mid)
    } else if (contrast == "RD" && distrib == "poi") {
      scor <- ftn(round(tan(pi * mid/2), 10))  
      # avoid machine precision producing values outside [-1,1]
    } else if (contrast %in% c("RR", "OR") || 
               (contrast == "p" && distrib == "poi")) {
      scor <- ftn(round(tan(pi * (mid + 1)/4), 10))  
      # avoid machine precision producing values outside [-1,1]
    } else if (contrast == "p" && distrib == "bin") {
      scor <- ftn((mid + 1)/2) 
    }
    check <- (scor < 0)  | is.na(scor) #??scor=NA only happens when |p1-p2|=1 and |theta|=1 (in which case hi==lo anyway), or if p1=p2=0
    hi[check] <- mid[check]
    lo[!check] <- mid[!check]
    niter <- niter + 1
  }
  if (uplow == "low") 
    best <- lo else best <- hi
  if (contrast == "RD" && distrib == "bin") {
    return(best)
  } else if ((contrast %in% c("RD") && distrib == "poi")) {
    return(tan(best * pi/2))
  } else if (contrast %in% c("RR", "OR")|| 
             (contrast == "p" && distrib == "poi")) {
    return(tan((best + 1) * pi/4))
  } else if (contrast == "p" && distrib =="bin") 
    return((best + 1)/2)
}

# Internal function
scoretheta <- function (
	#function to evaluate the score at a given value of theta, given the observed data
	#uses the MLE solution (and notation) given in F&M, extended in Laud2015
	#This function is vectorised in x1,x2
	theta,
	x1,
	n1,
	x2 = NULL,
	n2 = NULL,
	distrib = "bin",
	contrast = "RD",
	bcf = TRUE,
	skew = TRUE,
	cc = FALSE,
	stratified = FALSE,
	wt = NULL,
	weighting = "IVS",
	tdas = FALSE,	
	...
	) {

  nstrat <- length(n1)
	lambda <- switch(as.character(bcf), "TRUE" = (n1 + n2)/(n1 + n2 - 1),
	                 "FALSE" = 1)
	p1hat <- x1/n1
	p2hat <- x2/n2 
	x <- x1 + x2
	N <- n1 + n2

 	#RMLE of p1|theta depends on whether theta=RD or RR, and on whether a binomial
 	#or Poisson distribution is assumed Binomial RD
	if (contrast == "RD") {
	  Stheta <- (p1hat - p2hat) - theta
	  if (distrib == "bin") {
	    s0 <- theta # NB: Farrington & Manning notation (p1448) uses theta for N2/N1, and s_0 for RD
	    a <- N
	    b <- (n1 + 2 * n2) * theta - N - x
	    c_ <- (n2 * theta - N - 2 * x2) * theta + x
	    d <- x2 * theta * (1 - theta)
	    v <- (b/a/3)^3 - b * c_/(6 * a * a) + d/a/2
	    s <- sqrt(pmax(0, (b/a/3)^2 - c_/a/3))  
	    # NaNs produced when? - machine precision again?
	    u <- ifelse(v > 0, 1, -1) * s
	    w <- (pi + acos(pmax(-1, pmin(1, ifelse(u == 0 & v == 0, 0, v/u^3)))))/3
	    # avoids machine precision errors passing impossible values (outside (-1,1))
	    # to acos
	    p2d <- pmin(1, pmax(0, round(2 * u * cos(w) - b/a/3, 10))) 
	    # again, machine precision errors can give values a tiny fraction either
	    # side of 0 or 1
	    p1d <- pmin(1, pmax(0, p2d + theta))
	    V <- pmax(0, (p1d * (1 - p1d)/n1 + p2d * (1 - p2d)/n2) * lambda)
	    # V set to zero if machine precision error in cos function produces a negative V
	    mu3 <- p1d * (1 - p1d) * (1 - 2 * p1d)/(n1^2) -
	      p2d * (1 - p2d) * (1 - 2 * p2d)/(n2^2)
	  } else if (distrib == "poi") {
	    A <- N
	    B <- N * theta - x
	    C_ <- -x2 * theta
	    num <- (-B + sqrt(B^2 - 4 * A * C_))
	    p2d <- ifelse(num == 0, 0, num/(2 * A))
	    p1d <- p2d + theta
	    V <- (p1d/n1 + p2d/n2)
	    mu3 <- (p1d/(n1^2) - p2d/(n2^2))
	  }
	} else if (contrast == "RR") {
	  Stheta <- p1hat - p2hat * theta
	  if (distrib == "bin") {
	    A <- N * theta
	    B <- (-(n1 * theta + x1 + n2 + x2 * theta))
	    C_ <- x
	    num <- (-B - Re(sqrt(as.complex(B^2 - 4 * A * C_))))
	    p2d <- ifelse(A == 0, -C_/B, ifelse(num == 0, 0, num/(2 * A)))
	    p1d <- p2d * theta
	    V <- pmax(0, (p1d * (1 - p1d)/n1 + (theta^2) * p2d * (1 - p2d)/n2) * lambda)
	    V[is.na(V)] <- Inf
	    mu3 <- (p1d * (1 - p1d) * (1 - 2 * p1d)/(n1^2) -
	              (theta^3) * p2d * (1 - p2d) * (1 - 2 * p2d)/(n2^2))
	  } else if (distrib == "poi") {
	    p2d <- (x1 + x2)/(n1 * theta + n2)  #incidence density version from M&N p223
	    p1d <- p2d * theta
	    V <- pmax(0, (p1d/n1 + (theta^2) * p2d/n2))
	    mu3 <- (p1d/(n1^2) - (theta^3) * p2d/(n2^2))
	  }
	} else if (contrast == "OR") {
	  if (distrib == "bin") {
	    A <- n2 * (theta - 1)
	    B <- n1 * theta + n2 - x * (theta - 1)
	    C_ <- -x
	    num <- (-B + sqrt(B^2 - 4 * A * C_))
	    p2d <- ifelse(A == 0, -C_/B, ifelse(num == 0, 0, num/(2 * A))) 
	    # If A=0 then we solve a linear equation instead
	    p1d <- p2d * theta/(1 + p2d * (theta - 1))
	    Stheta <- (p1hat - p1d)/(p1d * (1 - p1d)) - (p2hat - p2d)/(p2d * (1 - p2d))
	    V <- pmax(0, (1/(n1 * p1d * (1 - p1d)) + 1/(n2 * p2d * (1 - p2d))) * lambda)  
	    # set to zero if machine precision error in cos function produces a negative V
	    mu3 <- (1 - 2 * p1d)/((n1 * p1d * (1 - p1d))^2) -
	      (1 - 2 * p2d)/((n2 * p2d * (1 - p2d))^2)
	    Stheta[(x1 == 0 & x2 == 0) | (x1 == n1 & x2 == n2)] <- 0
	    mu3[(x1 == 0 & x2 == 0) | (x1 == n1 & x2 == n2)] <- 0
	  } else if (distrib == "poi") {
	    print("Odds ratio not applicable to Poisson rates")
	  }
	} else if (contrast == "p") {
	  Stheta <- p1hat - theta
	  if (distrib == "bin") {
	    V <- (pmax(0, (theta * (1 - theta)/n1)))  
	    # set to zero if machine precision error in cos function produces a negative V
	    mu3 <- (theta * (1 - theta) * (1 - 2 * theta)/(n1^2))
	  } else if (distrib == "poi") {
	    V <- theta/n1
	    mu3 <- theta/(n1^2)
	  }
	  p1d <- theta
	  p2d <- NA
	}
	
	# continuity corrections
	corr <- 0
	if(cc > 0) {
	  if (contrast == "OR") {
	    corr <- cc * (1/(n1 * p1d * (1 - p1d)) + 1/(n2 * p2d * (1 - p2d)))  
	    # cc=0.5 gives Cornfield correction. Try cc=0.25 instead
	  } else if (contrast == "RR") {
	    corr <- cc * (1/(n1) + theta/(n2))  #try 0.125 or 0.25
	  } else if (contrast == "RD") {
	    corr <- cc * (1/pmin(n1, n2))  #cc=0.5 gives Hauck-Anderson. Try cc=0.25
  	  if (stratified == TRUE) corr <- (3/16) * (sum(n1 * n2/(n1 + n2)))^(-1)  
  	  # from Mehrotra & Railkar, also Zhao et al.
	  } else if (contrast == "p") {
	    corr <- cc/n1
	  }  
	}
		
	if (stratified == TRUE) {
	  pval <- NA
	  if (is.null(wt)) {
	    if (weighting == "MH") {
	      wt <- n1 * n2/(n1 + n2)  
	      # MH weights for RD, applied across other comparative parameters too 
	      # (without theoretical justification for OR)
	    } else if (weighting == "IVS") {
	      # IVS: inverse variance weights updated wih V.tilde
	      if (all(V == 0) || all(is.na(V)) || all (V == Inf)) { 
	        wt <- rep(1, nstrat)
	      } else wt <- 1/V  
	    } else if (weighting == "MN") {
	      if (contrast == "RR") {
	        # M&Ns iterative weights - quite similar to MH
	        wtx <- (1/n1 + theta/n2)^(-1)
	        p2ds <- sum(wtx * p2d)/sum(wtx)
	        p1ds <- sum(wtx * p1d)/sum(wtx)
	        wty <- ((1 - p1ds)/(n1 * (1 - p2ds)) + theta/n2)^(-1)
	        wty[p2ds == 1] <- 0
	        p2ds <- sum(wty * p2d)/sum(wty)
	        p1ds <- sum(wty * p1d)/sum(wty)
	        wt <- ((1 - p1ds)/(n1 * (1 - p2ds)) + theta/n2)^(-1)
	        wt[p2ds == 1] <- 0
	      } else if (contrast == "RD") {
	        # M&Ns iterative weights - quite similar to MH wtx <- n1*n2/(n1+n2)
	        wtx <- (1/n1 + 1/n2)^(-1)
	        p2ds <- sum(wtx * p2d)/sum(wtx)
	        p1ds <- sum(wtx * p1d)/sum(wtx)
	        wty <- ((p1ds * (1 - p1ds)/(p2ds * (1 - p2ds)))/n1 + 1/n2)^(-1)
	        wty[p2ds == 0] <- 0
	        p2ds <- sum(wty * p2d)/sum(wty)
	        p1ds <- sum(wty * p1d)/sum(wty)
	        wt <- ((p1ds * (1 - p1ds)/(p2ds * (1 - p2ds)))/n1 + 1/n2)^(-1)
	        # wt <- (1/n1 + (p2ds*(1-p2ds))/n2)^(-1) #equivalent, according to ...
	        wt[p2ds == 0] <- 0
	      } else if (contrast == "OR") {
	        # M&Ns weights are very similar in structure to IVS
	        wt <- n1 * n2 * ((1 - p1d) * p2d)^2/
	          (n1 * p1d * (1 - p1d) + n2 * p2d * (1 - p2d))
	      }
	    }
	  } else weighting <- "User-defined"

		Sdot <- sum(wt * Stheta)/sum(wt)
		if (weighting == "IVS") {
		  Q.i <- wt * ((Stheta - Sdot)^2)  #This version for iterative weights?
		} else Q.i <- ((Stheta - Sdot)^2)/V
		Q <- sum(Q.i)
		W <- sum(wt)
		
		if (weighting == "IVS") {
		  tau2 <- max(0, (Q - (nstrat - 1)))/(W - (sum(wt^2)/W))  
		  # published formula for IVS weights
		} else {
		  tau2 <- max(0, (Q - (nstrat - 1 * (sum(V * wt^2)/W))))/
		    (sum(1/V) - 1 * (sum(wt^2)/W))  
		  # only needed if want to output tau2.
		}
		if (tdas == TRUE && weighting == "IVS" && !(all(V == 0) || all(is.na(V)))) {
		  wt <- 1/(V + tau2)
		}  
		
		Sdot <- sum(wt * Stheta)/sum(wt)
		VS <- sum(wt * (Stheta - Sdot)^2)/((nstrat - 1) * sum(wt))
		
		if (tdas == TRUE) t2 <- tau2 else t2 <- 0
		Vdot <- sum(((wt/sum(wt))^2) * (V + t2))  
		# from equation (15) of Miettinen&Nurminen, with the addition of between
		# strata variance from Whitehead&Whitehead
		
		score1 <- sum((wt/sum(wt)) * (Stheta - corr))/pmax(0, sqrt(Vdot))
		scterm <- sum(((wt/sum(wt))^3) * mu3)/(6 * Vdot^(3/2))
		A <- scterm
		B <- 1
		C_ <- -(score1 + scterm)
		num <- (-B + sqrt(max(0, B^2 - 4 * A * C_)))
		score <- ifelse((skew == FALSE | scterm == 0), score1, num/(2 * A))
		pval <- pnorm(score)
		if (tdas == TRUE) {
		  score <- sum((wt/sum(wt)) * (Stheta - corr))/sqrt(VS)  
		  pval <- pt(score, nstrat - 1)
		}
		p2ds <- sum(wt * p2d/sum(wt))
		p1ds <- sum(wt * p1d/sum(wt))
	} 
	else if (stratified==FALSE) {
 		corr <- corr*sign(Stheta) 
		p1ds <- p1d
		p2ds <- p2d
		
		# Calculation of score & p-value involves solving 
		# z.p = Stheta/sqrt(V) -(z.p^2)*mu3/(6*V^(3/2)) + mu3/(6*V^(3/2)) 
		# Note that in the special case of mu3=0, this reduces to the skew=FALSE case.
		A <- mu3/(6 * V^(3/2))
		B <- 1
		C_ <- -((Stheta - corr)/sqrt(V) + mu3/(6 * V^(3/2)))
		num <- (-B + sqrt(pmax(0, B^2 - 4 * A * C_)))
		
		score <- ifelse((skew == FALSE | mu3 == 0 | (distrib == "poi" & abs(mu3) <= 10e-16)),
		                ifelse(Stheta == 0, 0, (Stheta - corr)/sqrt(V)), num/(2 * A)
		                )
		pval <- pnorm(score)
	}

	outlist <- list(score = score, p1d = p1d, Stheta = Stheta, num=num, V = V,
	                p2d = p2d, mu3 = mu3, pval = pval)
	if (stratified) {
	  outlist <- append(outlist, list(Sdot = Sdot, Vdot = Vdot, tau2 = tau2,
	             VS = VS, t2 = t2, Q.i = Q.i, Q = Q, wt = wt, p1ds = p1ds, p2ds = p2ds))
	}
	return(outlist)
}