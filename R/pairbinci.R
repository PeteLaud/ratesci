#' Confidence intervals for comparisons of paired binomial rates.
#' 
#' Score-based confidence intervals for the rate (or risk) difference ('RD'),
#' rate ratio ('RR') or odds ratio ('OR'), for paired binomial data.
#' [For paired Poisson rates, use the tdasci function with distrib='poi',
#' and weighting='MH', with pairs as strata.].
#' This function applies the Tango and Tang methods for RD and RR respectively, 
#' as well as an experimental method using the stratified TDAS method with 
#' pairs as strata. 
#' For OR, intervals are produced based on transforming various intervals for 
#' the single proportion, inclusing SCAS, mid-p and Jeffreys.
#' 
#' @param x A numeric vector object specified as c(a,b,c,d) 
#'   where:
#'   a is the number of pairs with the event (e.g. success) under both
#'     conditions (e.g. treated/untreated, or case/control)
#'   b is the count of the number with the event on condition 1 only (=n12)
#'   c is the count of the number with the event on condition 2 only (=n21)
#'   d is the number of pairs with no event under both conditions
#'   (Note the order of a and d is only important for contrast="RR".)
#' @param contrast Character string indicating the contrast of interest: 
#'   "RD" = rate difference (default), "RR" = rate ratio, "OR" = odds ratio.
#' @param level Number specifying confidence level (between 0 and 1, default 
#'   0.95).
#' @param theta0 Number to be used in a one-sided significance test (e.g. 
#'   non-inferiority margin). 1-sided p-value will be <0.025 iff 2-sided 95\% CI
#'   excludes theta0. NB: can also be used for a superiority test by setting 
#'   theta0=0.
#' @param method.RD Character string indicating the confidence interval method 
#'   to be used for contrast="RD". "Score" = Tango asymptotic score (default), 
#'   "TDAS" = t-distribution asymptotic score (experimental method, seems to
#'   struggle with low numbers).
#' @param method.RR Character string indicating the confidence interval method
#'   to be used for contrast="RR". "Score" = Tang asymptotic score (default),
#'   "TDAS" t-distribution asymptotic score (experimental method, seems to
#'   struggle with low numbers).
#' @param method.OR Character string indicating the confidence interval method
#'   to be used for contrast="OR", all of which are based on transformation of
#'   an interval for a single proportion b/(b+c):
#'   "SCAS" = transformed skewness-corrected score (default),
#'   "Jeffreys" = transformed Jeffreys (to be added),
#'   "midp" = transformed mid-p,
#'   ("Wilson" = transformed Wilson score - not yet included, would be for
#'   reference only, not recommended).
#' @param precis Number (default 6) specifying precision (i.e. number of decimal
#'   places) to be used in optimisation subroutine for the confidence interval. 
#   ...more parameters to be added: cc? skew??
#' @importFrom stats uniroot pbinom ppois dpois
#' @examples  
#'   #Data example from Agresti-Min 2005
#'   pairbinci(x = c(53,16,8,9), contrast="RD", method.RD="Score")
#'   pairbinci(x = c(53,16,8,9), contrast="RD", method.RD="TDAS")
#'   pairbinci(x = c(53,16,8,9), contrast="RR", method.RR="Score")
#'   pairbinci(x = c(53,16,8,9), contrast="RR", method.RR="TDAS")
#'   pairbinci(x = c(53,16,8,9), contrast="OR", method.OR="SCAS")
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#' @references 
#'   Laud PJ. Equal-tailed confidence intervals for comparison of 
#'   rates. Pharmaceutical Statistics 2017; 16:334-348.
#'   
#'   Tango T. Equivalence test and confidence interval for the difference
#'   in proportions for the paired-sample design. Statistics in Medicine
#'   1998; 17:891-908
#'   
#'   Tango T. Improved confidence intervals for the difference between binomial 
#'   proportions based on paired data by Robert G. Newcombe, Statistics in 
#'   Medicine, 17, 2635–2650 (1998). Statistics in Medicine 1999;
#'   18(24):3511-3513
#'   
#'   Tang N-S, Tang M-L, Chan ISF. On tests of equivalence via non-unity
#'   relative risk for matched-pair design. Statistics in Medicine 2003;
#'   22:1217-1233
#'   
#'   Fagerland MW, Lydersen S, Laake P. Recommended tests and
#'   confidence intervals for paired binomial proportions. Statistics in
#'   Medicine 2014; 33(16):2850–2875
#'   
#'   Agresti A, Min Y. Simple improved confidence intervals for
#'   comparing matched proportions. Statistics in Medicine 2005;
#'   24:729-740
#'   
#' @export

pairbinci <- function(
  x,
  contrast = "RD",
  level = 0.95,
  method.RD = "Score",
  method.RR = "Score",
  method.OR = "SCAS",
  theta0 = NULL,
  precis = 6
) {
  if (!(tolower(substr(contrast, 1, 2)) %in% c("rd", "rr", "or"))) {
    print("Contrast must be one of 'RD', 'RR' or 'OR'")
    stop()
  }   
  if (!is.numeric(c(x))) {
    print("Non-numeric inputs!")
    stop()
  }
  
  #Convert the data into 2 columns of 0s and 1s for use in TDAS-based method
  #and output 2x2 table for validation
  x1i <- rep(c(1, 1, 0, 0), x)
  x2i <- rep(c(1, 0, 1, 0), x)
  xi <- table(x1i, x2i)
  
   if (contrast =="OR") {
    #special case for OR, use conditional method based on transforming the 
    #SCAS interval for a single proportion
     b <- x[2] 
     c <- x[3] 
     if (method.OR == "SCAS") { #could add transformed Wilson version for reference
       trans.th0 <- NULL
       if(!is.null(theta0)) trans.th0 <- theta0/(1 + theta0)
       OR.ci <- scasci(x1 = b, n1 = b + c, contrast = "p", distrib = "bin", 
                      level = level, theta0 = trans.th0)
       estimates <- rbind(c(OR.ci$estimates[, c(1:3)]/(1 - OR.ci$estimates[, c(1:3)]), 
                     OR.ci$estimates[, 4]))
       pval <- OR.ci$pval
       outlist <- list(xi, estimates = estimates, pval = pval)
     } else if (method.OR == "midp") {
       trans.ci <- midpci(x = b, n = b + c, level = level)
       estimates <- c(trans.ci/(1 - trans.ci))
       outlist <- list(xi, estimates = estimates)
     }
  } else {
    if((contrast =="RD" && method.RD == "TDAS") ||
       (contrast =="RR" && method.RR == "TDAS")) {
      #stratified TDAS method for paired data as suggested in Laud 2017
      n1i <- n2i <- rep(1, sum(x))
      out <- tdasci(x1 = x1i, n1 = n1i, x2 = x2i, n2 = n2i, weighting = "MH",
                    contrast = contrast, distrib = "bin", level = level, 
                    theta0 = theta0, warn = FALSE)
      outlist <- list(xi, estimates = out$estimates, pval = out$pval)
    }
    #Score methods by Tango (for RD) & Tang (for RR):
    if((contrast =="RD" && method.RD == "Score") ||
       (contrast =="RR" && method.RR == "Score")) {
      myfun <- function(theta) {
        scorepair(theta = theta, x = x, contrast = contrast)$score  
      }
      # Use bisection routine to locate lower and upper confidence limits
      qtnorm <- qnorm(1 - (1 - level)/2)
      MLE <- bisect(ftn = function(theta) myfun(theta) - 0, distrib = "bin",
                    contrast = contrast, precis = precis + 1, uplow = "low")
      lower <- bisect(ftn = function(theta) myfun(theta) - qtnorm, distrib = "bin", 
                      precis = precis + 1, contrast = contrast, uplow = "low")
      upper <- bisect(ftn = function(theta) myfun(theta) + qtnorm, distrib = "bin", 
                      precis = precis + 1, contrast=contrast, uplow = "up")
      estimates <- cbind(Lower = lower, MLE = MLE, Upper = upper, level = level)
      
      # optionally add p-value for a test of null hypothesis: theta<=theta0
      # default value of theta0 depends on contrast
      if (contrast == "RD") {
        theta00 <- 0 
      } else theta00 <- 1
      if (is.null(theta0)) {
        theta0 <- theta00
      } 
      scorezero <- scorepair(theta = theta00, x = x, contrast = contrast)
      scorenull <- scorepair(theta = theta0, x = x, contrast = contrast)
      pval.left <- scorenull$pval
      pval.right <- 1 - pval.left
      chisq.zero <- scorezero$score^2
      pval2sided <- pchisq(chisq.zero, 1, lower.tail = FALSE)
      pval <- cbind(chisq = chisq.zero, pval2sided, theta0 = theta0,
                    scorenull = scorenull$score, pval.left, pval.right)
      
      outlist <- list(xi, estimates = estimates, pval = pval)
    }
    
    #Placeholder:
    #Add code for MOVER Jeffreys methods from Tang2010 (and add SCAS version)? 
    #For both RD & RR this appears to be inferior to the Score methods
#    if((contrast =="RD" && method.RD == "MOVER") ||
#       (contrast =="RR" && method.RR == "MOVER")) {
#    }

  }
  return(outlist)
}


# Internal function
scorepair <- function (
  #function to evaluate the score at a given value of theta, given the observed
  #data for paired binomial RD and RR
  #uses the MLE solution (and notation) given in Fagerland 2014 from 
  #Tango (1998/1999) & Tang (2003)
  #This function is not vectorised
  theta,
  x,
  contrast = "RD",
  ...
) {
  
  N <- sum(x)
  p1hat <- (x[1] + x[2])/N
  p2hat <- (x[1] + x[3])/N 

  if (contrast == "RD") { #notation per Tango 1999 letter
    Stheta <- ((x[2] - x[3]) - N * theta)
    A <- 2 * N
    B <- -x[2] - x[3] + (2*N - x[2] + x[3]) * theta
    C_ <- -x[3] * theta * (1 - theta) 
    num <- (-B + Re(sqrt(as.complex(B^2 - 4 * A * C_))))
    p2d <- ifelse(num == 0, 0, num/(2 * A))
    V <- pmax(0, N*(2 * p2d + theta * (1 - theta)))
  }
  if (contrast == "RR") { #per Tang 2003
    Stheta <- ((x[2] + x[1]) - (x[3] + x[1]) * theta)
    A <- N * (1 + theta)
    B <- (x[1] + x[3]) * theta^2 - (x[1] + x[2] + 2*x[3])
    C_ <- x[3] * (1 - theta) * (x[1] + x[2] + x[3])/N
    num <- (-B + Re(sqrt(as.complex(B^2 - 4 * A * C_))))
    q21 <- ifelse(num == 0, 0, num/(2 * A))
    q11 <- ((x[1] + x[2] + x[3])/N - (1 + theta) * q21)/theta
    q12 <- (q21 + (theta - 1) * (x[1] + x[2] + x[3])/N)/theta
    q1 <- q11 + q12
    q2 <- q1/theta
    V <- pmax(0, N * (1 + theta) * q21 + (x[1] + x[2] + x[3]) * (theta - 1))
  }
  score <- ifelse(Stheta == 0, 0, Stheta/sqrt(V))
  pval <- pnorm(score)
  outlist <- list(score = score, pval = pval)
  return(outlist)
}

# Internal function for both Clopper-Pearson and mid-p, binomial or Poisson
midpci <- function(
  #function to calculate exact 'mid-p' confidence interval for a single 
  #binomial or Poisson rate x/n
  x, 
  n, 
  level = 0.95, 
  exact = FALSE, 
  distrib = 'bin',
  precis = 8
  ) {

  alpha <- 1 - level
  cc <- (exact==FALSE) * 0.5
  if(distrib =='bin') {
    lowroot <- function (p) pbinom(x - 1, n, p) + cc * dbinom(x, n, p) - (1 - alpha/2)
    uproot <- function (p) pbinom(x, n, p) - cc * dbinom(x, n, p) - alpha/2
  } else if(distrib =='poi') {
    lowroot <- function (p) ppois(x, p) + cc * dpois(x, p) - (1 - alpha/2)
    uproot <- function (p) ppois(x, p) - cc * dpois(x, p) - alpha/2
  }
  lower <- bisect(ftn = lowroot, precis = precis, uplow = "low", contrast = 'p', distrib = distrib)
  upper <- bisect(ftn = uproot, precis = precis, uplow = "up", contrast = 'p', distrib = distrib)
  return(cbind(Lower = lower, Upper = upper)/ifelse(distrib=='poi',n,1))
}


if(FALSE) {
n
  alpha <- 0.0000
  ps <- seq(0,0.1,0.001)
  root <- lowroot(ps)
  plot(ps,root,type='l')
  max(root)

    #
  midpci(0:10,rep(10,11),level=0.000001)
  scasci.nonit(0:10,rep(10,11),level=0.999)
  scasci(0:10,rep(10,11),contrast='p',level=0.999)
  # midpci(1,29,exact=T)
  #midpci(1,29,exact=F)
  #midpci(1,29,exact=F,distrib='poi')
  
  # Internal function - previous version, precision not high enough
  midpci0 <- function(x, n, level = 0.95) {
    #function to calculate exact 'mid-p' confidence interval for a single proportion x/n
    alpha <- 1 - level
    lowroot <- function(p) {
      pbinom(x - 1, n, p, lower.tail = FALSE) - 0.5 * dbinom(x, n, p) - alpha/2
    }
    uproot <- function(p) {
      pbinom(x, n, p) - 0.5 * dbinom(x, n, p) - alpha/2
    }
    if (x == 0) {
      p.L <- 0
    } else  {
      p.L <- uniroot(f = lowroot, interval = c(0, 1))$root
    }
    if (x == n){
      p.U <- 1
    } else  {
      p.U <- uniroot(f = uproot, interval = c(0, 1))$root
    }
    return(c(p.L, p.U))
  }
  
  
    #  if (contrast == "RD") { #per Tango 1998 - something's flipped
  #    Stheta <- ((x[2] - x[3]) + N*theta)
  #    A <- 2 * N
  #    B <- -x[2] - x[3] - (2*N - x[2] + x[3])*theta
  #    C_ <- x[3] * theta * (theta + 1) 
  #    num <- (-B + sqrt(B^2 - 4 * A * C_))
  #    p2d <- ifelse(num == 0, 0, num/(2 * A))
  #    V <- pmax(0, N*(2 * p2d - theta * (theta + 1)))
  #  }  
  
  x <- c(1,1,7,12)
  level <- 0.95
  contrast="RD"
  myfun <- function(theta) {
    scorepair(theta = theta, x = x, contrast = contrast)$score  
  }
  # Use bisection routine to locate lower and upper confidence limits
  qtnorm <- qnorm(1 - (1 - level)/2)
  lower <- bisect(ftn = function(theta)
    myfun(theta) - qtnorm, precis=7, contrast=contrast,
    uplow="low")
  upper <- bisect(ftn = function(theta)
    myfun(theta) + qtnorm, precis=7, contrast=contrast,
    uplow="up")
  c(lower,upper)
  
  #RD example
  thetas <- seq(-0.6,0,0.01)
  scores <- sapply(1:length(thetas),function(i) scorepair(thetas[i],c(1,1,7,12),contrast="RD")$score)
  plot(thetas,scores)
  abline(h=c(-1,1)*1.96)
  abline(v=c(-0.517,-0.026))
  #scoreci.mp(7,1,21,0.95) #check using PropCIs
  
  #RR example
  thetas <- seq(0,1,0.01)
  scores <- sapply(1:length(thetas),function(i) scorepair(thetas[i],c(1,1,7,12),contrast="RR")$score)
  plot(thetas,scores)
  abline(h=c(-1,1)*1.96)
  abline(v=c(0.065,0.907))
  
  #pairbinci(x=c(1,1,7,12))
  #pairbinci(x=c(1,1,7,12),contrast="RR", method.RR="Score")
  #pairbinci(x=c(1,1,7,12),contrast="RR", method.RR="NB")
  #pairbinci(x=c(1,1,7,12),tang=TRUE,contrast="RD",precis=8)
  #pairbinci(x=c(1,1,7,12),tang=TRUE,contrast="RR",precis=8)
  
}    