#' MOVER confidence intervals for comparisons of independent binomial or Poisson
#' rates.
#'
#' Confidence intervals applying the MOVER method ("Method of Variance Estimates
#' Recovery", developed from the Newcombe method for binomial RD) across
#' different contrasts (RD, RR, OR) and distributions (binomial, Poisson) using
#' equal-tailed Jeffreys intervals instead of the Wilson score method for the
#' event rates.  Also allows more general Beta and Gamma priors for an
#' approximate Bayesian confidence interval incorporating prior beliefs about
#' the group event rates.
#'
#' @param x1,x2 Numeric vectors of numbers of events in group 1 & group 2
#'   respectively.
#' @param n1,n2 Numeric vectors of sample sizes (for binomial rates) or exposure
#'   times (for Poisson rates) in each group.
#' @param a1,b1,a2,b2 Numbers defining the Beta(ai,bi) prior distributions for
#'   each group (default ai = bi = 0.5 for Jeffreys method). Gamma priors for
#'   Poisson rates require only a1, a2.
#' @param cc Number or logical specifying (amount of) continuity correction
#'   (default FALSE). Numeric value is taken as the gamma parameter in Laud 2017,
#'   Appendix S2 (default 0.5 if cc=TRUE). Forced equal to 0.5 if type="exact".
#'   NB: cc currently not implemented for type="wilson".
#' @param contrast Character string indicating the contrast of interest: "RD" =
#'   rate difference (default), "RR" = rate ratio, "OR" = odds ratio.
#'   contrast="p" gives an interval for the single proportion x1/n1.
#' @param type Character string indicating the method used for the intervals for
#'   the individual group rates.
#'   "jeff" = Jeffreys equal-tailed intervals (default);
#'   "exact" = Clopper-Pearson/Garwood exact intervals (note this does NOT
#'   result in a strictly conservative interval for the contrast, except for
#'   contrast='p'. The scoreci function with cc=TRUE is recommended as a
#'   superior approximation of 'exact' methods);
#'   "midp" = mid-p intervals (to be added);
#'   "SCAS" = SCAS non-iterative intervals (to be added);
#'   "wilson" = Wilson score intervals (as per Newcombe 1998).
#'   NB: "wilson" option is included only for legacy validation against previous
#'   published method by Newcombe. It is not recommended, as type="jeff"
#'   achieves much better coverage properties.
#' @param adj Logical (default FALSE) indicating whether to apply the boundary
#'   adjustment for Jeffreys intervals recommended on p108 of Brown et al.
#'   (type = "jeff" only: set to FALSE if using informative priors)
#' @param ... Additional arguments.
#' @inheritParams jeffreysci
#' @importFrom stats pchisq pf pnorm pt qbeta qgamma qnorm qqnorm qt
#' @importFrom graphics abline lines text
#' @return A matrix containing the confidence interval for the requested contrast
#' @examples
#' # Binomial RD, MOVER-J method:
#' moverci(x1 = 5, n1 = 56, x2 = 0, n2 = 29)
#'
#' # Binomial RD, Newcombe method:
#' moverci(x1 = 5, n1 = 56, x2 = 0, n2 = 29, type = "wilson")
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#' @references
#'   Laud PJ. Equal-tailed confidence intervals for comparison of
#'   rates. Pharmaceutical Statistics 2017; 16:334-348.
#'
#'   Newcombe RG. Interval estimation for the difference between independent
#'   proportions: comparison of eleven methods. Statistics in Medicine 1998;
# ;   17(8):873-890.
#'
#'   Donner A, Zou G. Closed-form confidence intervals for functions of the
#'   normal mean and standard deviation. Statistical Methods in Medical Research
# ;   2012; 21(4):347-359.
#'
#'   Fagerland MW, Newcombe RG. Confidence intervals for odds ratio and relative
#'   risk based on the inverse hyperbolic sine transformation. Statistics in
#'   Medicine 2013; 32(16):2823-2836.
#'
#'   Li HQ, Tang ML, Wong WK. Confidence intervals for ratio of two Poisson
#'   rates using the method of variance estimates recovery. Computational
#'   Statistics 2014; 29(3-4):869-889.
#'
#' @export
moverci <- function(x1,
                    n1,
                    x2 = NULL,
                    n2 = NULL,
                    a1 = 0.5,
                    b1 = 0.5,
                    a2 = 0.5,
                    b2 = 0.5,
                    cc = FALSE,
                    level = 0.95,
                    distrib = "bin",
                    contrast = "RD",
                    type = "jeff",
                    adj = FALSE,
                    ...) {
  if (!(tolower(substr(type, 1, 4)) %in% c("jeff", "wils", "exac"))) {
    print("Type must be one of 'jeffreys', 'wilson' or 'exact'")
    stop()
  }
  if (!(tolower(substr(distrib, 1, 3)) %in% c("bin", "poi"))) {
    print("Distrib must be one of 'bin' or 'poi'")
    stop()
  }
  if (!(tolower(substr(contrast, 1, 2)) %in% c("rd", "rr", "or", "p"))) {
    print("Contrast must be one of 'RD', 'RR', 'OR' or 'p'")
    stop()
  }
  if (contrast != "p" && (is.null(x2) || is.null(n2))) {
    print("Argument x2 or n2 missing")
    stop()
  }
  if (!is.numeric(c(x1, n1, x2, n2))) {
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
  if (distrib == "poi" && contrast == "OR") {
    print("Odds ratio not applicable to Poisson rates")
    stop()
  }
  if (adj == TRUE & !(tolower(substr(type, 1, 4)) == c("jeff"))) {
    print("adj only applicable for type = 'jeff'")
    stop()
  }
  if (as.character(cc) == "TRUE") cc <- 0.5

  alpha <- 1 - level
  z <- qnorm(1 - alpha / 2)

  lenx <- length(x1)
  # in case x1,x2 are input as vectors but n1,n2 are not
  if (length(n1) < lenx && lenx > 1) n1 <- rep(n1, length.out = lenx)
  if (length(n2) < lenx && lenx > 1) n2 <- rep(n2, length.out = lenx)

  if (contrast == "OR" && type == "wilson") {
    # special cases for OR handled as per Fagerland & Newcombe Table III
    # - not needed for MOVER-J now that xi/ni is not used for p.hat
    # but still needed for legacy Wilson method
    special <- (x2 == n2) | (x1 == n1)
    xx <- x2
    x2[special] <- (n1 - x1)[special]
    x1[special] <- (n2 - xx)[special]
    nx <- n1
    n1[special] <- n2[special]
    n2[special] <- nx[special]
  }

  if (type %in% c("wilson", "p")) {
    p1hat <- x1 / n1
    if (type != "p") p2hat <- x2 / n2
  }

  if (type %in% c("jeff", "exact")) {
    if (type == "exact") cc <- 0.5
    # MOVER-J, including optional 'continuity correction'
    j1 <- jeffreysci(x1, n1,
      ai = a1, bi = b1, cc = cc, level = level,
      distrib = distrib, adj = adj
    )
    if (contrast != "p") {
      j2 <- jeffreysci(x2, n2,
        ai = a2, bi = b2, cc = cc, level = level,
        distrib = distrib, adj = adj
      )
    } else {
      j2 <- NULL
    }
    p1hat <- j1[, 3]
    p2hat <- j2[, 3]
  } else if (type == "wilson") {
    # or use Wilson intervals as per Newcombe 1998
    # (NB could add cc here for completeness)
    j1 <- quadroot(
      a = 1 + z^2 / n1, b = -(2 * p1hat + z^2 / n1),
      c_ = p1hat^2
    )
    if (contrast != "p") {
      j2 <- quadroot(
        a = 1 + z^2 / n2, b = -(2 * p2hat + z^2 / n2),
        c_ = p2hat^2
      )
    } else {
      j2 <- NULL
    }
  }
  l1 <- j1[, 1]
  u1 <- j1[, 2]
  l2 <- j2[, 1]
  u2 <- j2[, 2]

  if (contrast == "p") {
    lower <- l1
    upper <- u1
    est <- p1hat
  } else if (contrast == "RD") {
    # From Newcombe (1998), p876 "Method 10"
    lower <- p1hat - p2hat - sqrt(pmax(0, (p1hat - l1)^2 + (u2 - p2hat)^2))
    upper <- p1hat - p2hat + sqrt(pmax(0, (u1 - p1hat)^2 + (p2hat - l2)^2))
    est <- p1hat - p2hat
    if (adj == TRUE) {
      # This is optional for informative priors
      lower[(x1 == 0 & x2 == n2)] <- -1
      upper[(x1 == n1 & x2 == 0)] <- 1
    }
  } else if (contrast == "OR") {
    # From Fagerland & Newcombe (2013), p2828 "Method 4"
    q1hat <- p1hat / (1 - p1hat)
    q2hat <- p2hat / (1 - p2hat)
    L1 <- l1 / (1 - l1)
    U1 <- u1 / (1 - u1)
    L2 <- l2 / (1 - l2)
    U2 <- u2 / (1 - u2)
    lower <- pmax(0, (q1hat * q2hat - sqrt(pmax(0, (q1hat * q2hat)^2 -
      L1 * U2 * (2 * q1hat - L1) * (2 * q2hat - U2)))) /
      (U2 * (2 * q2hat - U2)))
    upper <- (q1hat * q2hat + sqrt(pmax(0, (q1hat * q2hat)^2 -
      U1 * L2 * (2 * q1hat - U1) * (2 * q2hat - L2)))) /
      (L2 * (2 * q2hat - L2))
    est <- p1hat * (1 - p2hat) / (p2hat * (1 - p1hat))
    if (adj == TRUE) {
      # This is optional for informative priors
      upper[x2 == 0] <- Inf
      lower[(x1 == 0 & x2 == n2) | (x1 == n1 & x2 == 0)] <- 0
      upper[(x1 == 0 & x2 == n2) | (x1 == n1 & x2 == 0)] <- Inf
    }
    if (cc > 0) lower[u2 == 1] <- 0 # Avoids LCL=NaN
  } else if (contrast == "RR") {
    # From Donner & Zou (2012), p351
    # or Li et al. (2014), p873
    lower <- (p1hat * p2hat - sqrt(pmax(0, (p1hat * p2hat)^2 -
      l1 * (2 * p2hat - u2) * (u2 * (2 * p1hat - l1))))) /
      (u2 * (2 * p2hat - u2))
    upper <- (p1hat * p2hat + sqrt(pmax(0, (p1hat * p2hat)^2 -
      u1 * (2 * p2hat - l2) * (l2 * (2 * p1hat - u1))))) /
      (l2 * (2 * p2hat - l2))
    est <- p1hat / p2hat
    if (adj == TRUE) {
      # This is optional for informative priors
      upper[x2 == 0] <- Inf
      lower[(x1 == 0)] <- 0
    }
  }
  CI <- cbind(Lower = lower, Estimate = est, Upper = upper)
  CI
}

#' Jeffreys and other approximate Bayesian confidence intervals for a single
#' binomial or Poisson rate.
#'
#' Generalised approximate Bayesian confidence intervals based on a Beta (for
#' binomial rates) or Gamma (for Poisson rates) conjugate priors. Encompassing
#' the Jeffreys method (with Beta(0.5, 0.5) or Gamma(0.5) respectively), as well
#' as any user-specified prior distribution. Clopper-Pearson method (as
#' quantiles of a Beta distribution as described in Brown et al. 2001) also
#' included by way of a "continuity correction" parameter.
#'
#' @param x Numeric vector of number of events.
#' @param n Numeric vector of sample sizes (for binomial rates) or exposure
#'   times (for Poisson rates).
#' @param ai,bi Numbers defining the Beta prior distribution (default ai = bi =
#'   0.5 for Jeffreys interval). Gamma prior for Poisson rates requires only ai.
#' @param cc Number or logical specifying (amount of) "continuity correction".
#'   cc = 0 (default) gives Jeffreys interval, cc = 0.5 gives the
#'   Clopper-Pearson interval (or Garwood for Poisson). A value between 0 and
#'   0.5 allows a compromise between proximate and conservative coverage.
#' @param level Number specifying confidence level (between 0 and 1, default
#'   0.95).
#' @param distrib Character string indicating distribution assumed for the input
#'   data: "bin" = binomial (default), "poi" = Poisson.
#' @param adj Logical (default TRUE) indicating whether to apply the boundary
#'   adjustment recommended on p108 of Brown et al. (set to FALSE if informative
#'   priors are used)
#' @param ... Other arguments.
#' @importFrom stats qbeta qgamma qnorm
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#' @examples
#' # Jeffreys method:
#' jeffreysci(x = 5, n = 56)
#' @references
#'   Laud PJ. Equal-tailed confidence intervals for comparison of
#'   rates. Pharmaceutical Statistics 2017; 16:334-348.
#'
#'   Brown LD, Cai TT, DasGupta A. Interval estimation for a binomial
#'   proportion. Statistical Science 2001; 16(2):101-133
#' @export
jeffreysci <- function(x,
                       n,
                       ai = 0.5,
                       bi = 0.5,
                       cc = 0,
                       level = 0.95,
                       distrib = "bin",
                       adj = TRUE,
                       ...) {
  if (!is.numeric(c(x, n))) {
    print("Non-numeric inputs!")
    stop()
  }
  if (any(c(x, n) < 0)) {
    print("Negative inputs!")
    stop()
  }
  if (distrib == "bin" && (any(x > n + 0.001))) {
    print("x > n not possible for distrib = 'bin'")
    stop()
  }
  if (as.character(cc) == "TRUE") cc <- 0.5

  alpha <- 1 - level
  if (distrib == "bin") {
    CI_lower <- qbeta(alpha / 2, x + (ai - cc), n - x + (bi + cc))
    est <- qbeta(0.5, x + (ai), n - x + (bi)) # Obtain phat as the median
    CI_upper <- qbeta(1 - alpha / 2, x + (ai + cc), n - x + (bi - cc))
    if (adj == TRUE) {
      # adjustment at boundary values
      CI_lower[x == 0] <- 0
      CI_upper[x == n] <- 1
      est[x == 0] <- 0
      est[x == n] <- 1
    }
  } else if (distrib == "poi") {
    # Jeffreys prior for Poisson rate uses gamma distribution,
    # as defined in Li et al. with "continuity correction" from Laud 2017.
    CI_lower <- qgamma(alpha / 2, x + (ai - cc), scale = 1 / n)
    est <- qgamma(0.5, x + (ai), scale = 1 / n)
    if (adj == TRUE) {
      CI_lower[x == 0] <- 0
      est[x == 0] <- 0
    }
    CI_upper <- qgamma(1 - alpha / 2, (x + (ai + cc)), scale = 1 / n)
  }
  CI <- cbind(Lower = CI_lower, Upper = CI_upper, est = est)
  CI
}

#' Approximate Bayesian ("MOVER-B") confidence intervals for
#' comparisons of independent binomial or Poisson rates.
#'
#' Wrapper function for the MOVER-B methods.  Approximate Bayesian confidence
#' intervals for the rate (or risk) difference ("RD") or ratio ("RR") for
#' independent binomial or Poisson rates, or for odds ratio ("OR", binomial
#' only). (developed from Newcombe, Donner & Zou, Li et al, and Fagerland &
#' Newcombe, and generalised as "MOVER-B" in forthcoming publication) including
#' special case "MOVER-J" using non-informative priors with optional continuity
#' correction.  This function is vectorised in x1, x2, n1, and n2.
#'
#' @param x1,x2 Numeric vectors of numbers of events in group 1 & group 2
#'   respectively.
#' @param n1,n2 Numeric vectors of sample sizes (for binomial rates) or exposure
#'   times (for Poisson rates) in each group.
#' @param a1,b1,a2,b2 Numbers defining the Beta(ai,bi) prior distributions for
#'   each group (default ai = bi = 0.5 for Jeffreys uninformative priors). Gamma
#'   priors for Poisson rates require only a1, a2.
#' @inheritParams moverci
#' @export
moverbci <- function(x1,
                     n1,
                     x2,
                     n2,
                     a1 = 0.5,
                     b1 = 0.5,
                     a2 = 0.5,
                     b2 = 0.5,
                     distrib = "bin",
                     contrast = "RD",
                     level = 0.95,
                     cc = 0,
                     ...) {
  moverci(
    x1 = x1,
    n1 = n1,
    x2 = x2,
    n2 = n2,
    a1 = a1,
    b1 = b1,
    a2 = a2,
    b2 = b2,
    distrib = distrib,
    contrast = contrast,
    level = level,
    cc = cc,
    type = "jeff",
    adj = FALSE,
    ...
  )
}

# Internal function
quadroot <- function(a, b, c_) {
  # GET ROOTS OF A QUADRATIC EQUATION
  r1x <- (-b + sqrt(b^2 - 4 * a * c_)) / (2 * a)
  r2x <- (-b - sqrt(b^2 - 4 * a * c_)) / (2 * a)
  r1 <- pmin(r1x, r2x)
  r2 <- pmax(r1x, r2x)
  cbind(r1, r2)
}
