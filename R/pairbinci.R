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
#' the single proportion, including SCAS, mid-p and Jeffreys.
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
#' @param method_RD Character string indicating the confidence interval method
#'   to be used for contrast="RD". "Score" = Tango asymptotic score (default),
#'   "TDAS" = t-distribution asymptotic score (experimental method, seems to
#'   struggle with low numbers).
#' @param method_RR Character string indicating the confidence interval method
#'   to be used for contrast="RR". "Score" = Tang asymptotic score (default),
#'   "TDAS" t-distribution asymptotic score (experimental method, seems to
#'   struggle with low numbers).
#' @param method_OR Character string indicating the confidence interval method
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
#'   pairbinci(x = c(53,16,8,9), contrast="RD", method_RD="Score")
#'   pairbinci(x = c(53,16,8,9), contrast="RD", method_RD="TDAS")
#'   pairbinci(x = c(53,16,8,9), contrast="RR", method_RR="Score")
#'   pairbinci(x = c(53,16,8,9), contrast="RR", method_RR="TDAS")
#'   pairbinci(x = c(53,16,8,9), contrast="OR", method_OR="SCAS")
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
  method_RD = "Score",
  method_RR = "Score",
  method_OR = "SCAS",
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

   if (contrast == "OR") {
    #special case for OR, use conditional method based on transforming the
    #SCAS interval for a single proportion
     b <- x[2]
     c <- x[3]
     if (method_OR == "SCAS") {
       #could add transformed Wilson version for reference
       trans_th0 <- NULL
       if (!is.null(theta0)) trans_th0 <- theta0/(1 + theta0)
       OR_ci <- scasci(x1 = b, n1 = b + c, contrast = "p", distrib = "bin",
                      level = level, theta0 = trans_th0)
       estimates <- rbind(
         c(OR_ci$estimates[, c(1:3)]/(1 - OR_ci$estimates[, c(1:3)]),
          OR_ci$estimates[, 4]))
       pval <- OR_ci$pval
       outlist <- list(xi, estimates = estimates, pval = pval)
     } else if (method_OR == "midp") {
       trans_ci <- exactci(x = b, n = b + c, midp = TRUE, level = level)
       estimates <- c(trans_ci/(1 - trans_ci))
       outlist <- list(xi, estimates = estimates)
     }
  } else {
    if ((contrast =="RD" && method_RD == "TDAS") ||
       (contrast =="RR" && method_RR == "TDAS")) {
      #stratified TDAS method for paired data as suggested in Laud 2017
      n1i <- n2i <- rep(1, sum(x))
      out <- tdasci(x1 = x1i, n1 = n1i, x2 = x2i, n2 = n2i, weighting = "MH",
                    contrast = contrast, distrib = "bin", level = level,
                    theta0 = theta0, warn = FALSE)
      outlist <- list(xi, estimates = out$estimates, pval = out$pval)
    }
    #Score methods by Tango (for RD) & Tang (for RR):
    if ((contrast =="RD" && method_RD == "Score") ||
       (contrast =="RR" && method_RR == "Score")) {
      myfun <- function(theta) {
        scorepair(theta = theta, x = x, contrast = contrast)$score
      }
      # Use bisection routine to locate lower and upper confidence limits
      qtnorm <- qnorm(1 - (1 - level)/2)
      MLE <- bisect(ftn = function(theta) myfun(theta) - 0, distrib = "bin",
                    contrast = contrast, precis = precis + 1, uplow = "low")
      lower <- bisect(ftn = function(theta) myfun(theta) - qtnorm,
                    distrib = "bin", precis = precis + 1, contrast = contrast,
                    uplow = "low")
      upper <- bisect(ftn = function(theta) myfun(theta) + qtnorm,
                    distrib = "bin", precis = precis + 1, contrast=contrast,
                    uplow = "up")
      estimates <- cbind(Lower = lower, MLE = MLE, Upper = upper,
                         level = level)

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
      pval_left <- scorenull$pval
      pval_right <- 1 - pval_left
      chisq_zero <- scorezero$score^2
      pval2sided <- pchisq(chisq_zero, 1, lower.tail = FALSE)
      pval <- cbind(chisq = chisq_zero, pval2sided, theta0 = theta0,
                    scorenull = scorenull$score, pval_left, pval_right)

      outlist <- list(xi, estimates = estimates, pval = pval)
    }

    #Placeholder:
    #Add code for MOVER Jeffreys methods from Tang2010 (and add SCAS version)?
    #For both RD & RR this appears to be inferior to the Score methods
#    if ((contrast =="RD" && method_RD == "MOVER") ||
#       (contrast =="RR" && method_RR == "MOVER")) {
#    }

  }
  return(outlist)
}


# Internal function
scorepair <- function(
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
  if (contrast == "RD") {
    #notation per Tango 1999 letter
    Stheta <- ((x[2] - x[3]) - N * theta)
    A <- 2 * N
    B <- -x[2] - x[3] + (2*N - x[2] + x[3]) * theta
    C_ <- -x[3] * theta * (1 - theta)
    num <- (-B + Re(sqrt(as.complex(B^2 - 4 * A * C_))))
    p2d <- ifelse(num == 0, 0, num/(2 * A))
    V <- pmax(0, N*(2 * p2d + theta * (1 - theta)))
  }
  if (contrast == "RR") {
    #per Tang 2003
    Stheta <- ((x[2] + x[1]) - (x[3] + x[1]) * theta)
    A <- N * (1 + theta)
    B <- (x[1] + x[3]) * theta^2 - (x[1] + x[2] + 2 * x[3])
    C_ <- x[3] * (1 - theta) * (x[1] + x[2] + x[3])/N
    num <- (-B + Re(sqrt(as.complex(B^2 - 4 * A * C_))))
    q21 <- ifelse(num == 0, 0, num/(2 * A))
    V <- pmax(0, N * (1 + theta) * q21 + (x[1] + x[2] + x[3]) * (theta - 1))
  }
  score <- ifelse(Stheta == 0, 0, Stheta/sqrt(V))
  pval <- pnorm(score)
  outlist <- list(score = score, pval = pval)
  return(outlist)
}
