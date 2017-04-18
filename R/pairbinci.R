#' Skewness-corrected confidence intervals for comparisons of paired binomial
#' rates.
#' 
#' Score-based confidence intervals for the rate (or risk) difference ("RD"), rate 
#' ratio ("RR") or for odds ratio ("OR"), for paired binomial data. 
#' [For paired Poisson rates, use the tdasci function with distrib="poi",
#' with pairs as strata.].
#' This function applies the stratified TDAS method for RD and RR. For OR, an 
#' interval is produced based on transforming a SCAS interval for the single 
#' proportion.
#' 
#' @param x A numeric vector object specified as c(a,b,c,d) 
#'   where:
#'   a is the number of pairs with the event (e.g. success) under both conditions 
#'     (e.g. treated/untreated, or case/control)
#'   b is the count of the number with the event on condition 1 only 
#'   c is the count of the number with the event on condition 2 only
#'   d is the number of pairs with no event under both conditions
#'   (Note the order of a and d is only important for contrast="RR".)
#'   Note for data in columns of success/failure, use the tdasci function instead
#'   (not recommended for contrast="OR")
#' @param contrast Character string indicating the contrast of interest: 
#'   "RD" = rate difference (default), "RR" = rate ratio, "OR" = odds ratio.
#' @param level Number specifying confidence level (between 0 and 1, default 
#'   0.95).
#' @param delta Number to be used in a one-sided significance test (e.g. 
#'   non-inferiority margin). 1-sided p-value will be <0.025 iff 2-sided 95\% CI
#'   excludes delta. NB: can also be used for a superiority test by setting 
#'   delta=0.
#   ...more parameters to be added: cc?, precis
#' @examples  
#'   #Data example from Agresti-Min 2005
#'   pairbinci(c(53,16,8,9),contrast="RD")
#'   pairbinci(c(53,16,8,9),contrast="RR")
#'   pairbinci(c(53,16,8,9),contrast="OR")
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#' @references 
#'   Laud PJ. Equal-tailed confidence intervals for comparison of 
#'   rates. Pharmaceutical Statistics [in press].
#'   
#'   Fagerland MW, Lydersen S, Laake P. Recommended tests and
#'   confidence intervals for paired binomial proportions. Statistics in
#'   Medicine 2014; 33(16):2850–2875
#' @export

pairbinci <- function(
  x,
  contrast="RD",
  level=0.95,
  delta=NULL
  ) {
  if (!(tolower(substr(contrast, 1, 2)) %in% c("rd", "rr", "or"))) {
    print("Contrast must be one of 'RD', 'RR' or 'OR'")
    stop()
  }   
  if (!is.numeric(c(x))) {
    print("Non-numeric inputs!")
    stop()
  }
  #Convert the data into 2 columns of 0s and 1s
  #and output 2x2 table for validation
  x1i <- rep(c(1,1,0,0),x)
  x2i <- rep(c(1,0,1,0),x)
  xi <- table(x1i,x2i)
  
  if (contrast =="OR") {
    #special case for OR, use conditional method based on transforming the 
    #SCAS interval for a single proportion
    b <- x[2] 
    c <- x[3] 
    trans.del <- NULL
    if(!is.null(delta)) trans.del <- delta/(1+delta)
    OR.ci <- scasci(b,b+c,contrast="p",level=level,delta=trans.del)
    estimates <- c(OR.ci$estimates[,c(1:3)]/(1-OR.ci$estimates[,c(1:3)]),OR.ci$estimates[,4])
    pval <- OR.ci$pval
    outlist <- list(xi,estimates=estimates,pval=pval)
  } else {
    #for paired RD and RR, use stratified TDAS method
    n1i <- n2i <- rep(1,sum(x))
    out <- tdasci(x1=x1i,n1=n1i,x2=x2i,n2=n2i,weighting="MH",contrast=contrast,level=level,delta=delta,warn=FALSE)
    outlist <- list(xi,estimates=out$estimates,pval=out$pval)
  }
  outlist
}

