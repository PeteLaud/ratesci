

#' MOVER confidence intervals for comparison of independent binomial or Poisson 
#' rates.
#' 
#' Confidence intervals applying the MOVER method across different contrasts and distributions
#' using Jeffreys intervals (& other more general Beta and Gamma priors) instead of Wilson
#' 
#' @param x1,x2 Numeric vectors of numbers of events in group 1 & group 2 
#'   respectively.
#' @param n1,n2 Numeric vectors of sample sizes (or exposure times for Poisson 
#'   rates) in each group.
#' @param a1,b1 Numbers defining the Beta prior distribution for group 1 (0.5 
#'   for Jeffreys prior).
#' @param a2,b2 Numbers defining the Beta prior distribution for group 2 (0.5 
#'   for Jeffreys prior).
#' @param cc Number or logical specifying (amount of) continuity correction.
#' @param level Number specifying confidence level.
#' @param dist Character string indicating distribution of data: "bin"=binomial,
#'   "poi"=Poisson.
#' @param contrast Character string indicating the contrast required: ("RD", 
#'   "RR", or "OR").
#' @param type Character string indicating the method used for the intervals for
#'   the individual group rates. "Jeff"=Jeffreys equal-tailed intervals,
#'   "exact"=Clopper-Pearson exact intervals (also obtained using "Jeff" with
#'   cc=0.5), "Wilson"=Wilson score intervals (as per Newcombe 1998).
#' @param ... Additional arguments.
#' @export
MOVERCI <- function(
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
  type="Jeff",
  ...
){
  alpha <- 1 - level
  z <- qnorm(1 - alpha / 2)
  if (as.character(cc) == "TRUE") cc <- 0.5
  
  #in case x1,x2 are vectors but n1,n2 are not
  if (length(n1) == 1 & length(x1) > 1) n1 <- rep(n1, length(x1))
  if (length(n2) == 1 & length(x1) > 1) n2 <- rep(n2, length(x1))
  
  if (contrast == "OR") {
    # special cases for OR handled as per Fagerland & Newcombe Table III
    special <- ( (x2 == n2) | (x1 == n1))
    xx <- x2
    x2[special] <- (n1 - x1)[special]
    x1[special] <- (n2 - xx)[special]
    nx <- n1
    n1[special] <- n2[special]
    n2[special] <- nx[special]
  }
  
  p1hat <- x1 / n1
  p2hat <- x2 / n2
  
  if (contrast == "OR" && dist != "bin") {
    print("WARNING: Odds Ratio must use dist='bin'")
    dist <- "bin"
  }
  
  if (type == "Jeff") {
    # MOVER-J, including optional 'continuity correction' g
    j1 <- JeffreysCI(x1, n1, ai = a1, bi = b1, g = cc, alpha = alpha,
                     dist = dist, adj = paste( (contrast == "OR") ))
    j2 <- JeffreysCI(x2, n2, ai = a2, bi = b2, g = cc, alpha = alpha,
                     dist = dist, adj = paste( (contrast == "OR") ))
  } else if (type == "exact") {
    # MOVER-E based on Clopper-Pearson exact intervals
    j1 <- JeffreysCI(x1, n1, ai = a1, bi = b1, g = 0.5, alpha = alpha,
                     dist = dist, adj = paste( (contrast == "OR") ))
    j2 <- JeffreysCI(x2, n2, ai = a2, bi = b2, g = 0.5, alpha = alpha,
                     dist = dist, adj = paste( (contrast == "OR") ))
  } else {
    # or use Wilson intervals as per Newcombe 1988 (NB could add cc here for completeness)
    j1 <- quadroot(a = 1 + z ^ 2 / n1, b = - (2 * p1hat + z ^ 2 / n1),
                   c = p1hat ^ 2)
    j2 <- quadroot(a = 1 + z ^ 2 / n2, b = - (2 * p2hat + z ^ 2 / n2),
                   c = p2hat ^ 2)
  }
  l1 <- j1[, 1]
  u1 <- j1[, 2]
  l2 <- j2[, 1]
  u2 <- j2[, 2]
  
  if (contrast == "RD") {
    # From Newcombe 1998
    lower <- p1hat - p2hat -
      Re(sqrt(as.complex( (p1hat - l1) ^ 2 + (u2 - p2hat) ^ 2)))
    upper <- p1hat - p2hat +
      Re(sqrt(as.complex( (u1 - p1hat) ^ 2 + (p2hat - l2) ^ 2)))
  } else if (contrast == "OR") {
    # From Fagerland & Newcombe 2013
    q1hat <- p1hat / (1 - p1hat)
    q2hat <- p2hat / (1 - p2hat)
    L1 <- l1 / (1 - l1)
    U1 <- u1 / (1 - u1)
    L2 <- l2 / (1 - l2)
    U2 <- u2 / (1 - u2)
    lower <- pmax(0, (q1hat * q2hat -
                        Re(sqrt(as.complex( (q1hat * q2hat) ^ 2 -
                                              L1 * U2 * (2 * q1hat - L1) * (2 * q2hat - U2))))) /
                    (U2 * (2 * q2hat - U2)))
    upper <- (q1hat * q2hat +
                Re(sqrt(as.complex( (q1hat * q2hat) ^ 2 - U1 * L2 *
                                      (2 * q1hat - U1) * (2 * q2hat - L2))))) /
      (L2 * (2 * q2hat - L2))
    upper[x2 == 0] <- Inf
    lower[(x1 == 0 & x2 == n2) | (x1 == n1 & x2 == 0)] <- 0
    upper[(x1 == 0 & x2 == n2) | (x1 == n1 & x2 == 0)] <- Inf
  } else if (contrast == "RR") {
    # From Donner & Zou 2012 / Li et al
    lower <- (p1hat * p2hat -
                Re(sqrt(as.complex( (p1hat * p2hat) ^ 2 - l1 * (2 * p2hat - u2) *
                                      (u2 * (2 * p1hat - l1)))))) / (u2 * (2 * p2hat - u2))
    upper <- (p1hat * p2hat +
                Re(sqrt(as.complex( (p1hat * p2hat) ^ 2 - u1 * (2 * p2hat - l2) *
                                      (l2 * (2 * p1hat - u1)))))) / (l2 * (2 * p2hat - l2))
    upper[x2 == 0] <- Inf
  }
  CI <- cbind(Lower = lower, Upper = upper)
  CI
}

# cc=continuity correction: cc=0.5 gives Clopper-Pearson 'exact' (conservative) interval
# "adj" is the option to use exact binomial interval for zero cell counts as suggested by Brown et al
JeffreysCI <- function(
  x,
  n,
  ai=0.5,
  bi=0.5,
  cc=0,
  level=0.95,
  dist="bin",
  adj=FALSE,
  exact=F,
  ...
  ) {
  alpha <- 1 - level
  if (dist == "bin") {
    CI.lower <- qbeta( alpha / 2, x + (ai - cc), n - x + (bi + cc))
    CI.lower[x == 0] <- 0
    CI.upper <- qbeta(1 - alpha / 2, x + (ai + cc), n - x + (bi - cc))
    CI.upper[x == n] <- 1
    if (adj == TRUE) {
      CI.lower[x == n] <- ( (1 - level) / 2) ^ ( 1 / n )[x == n]
      CI.upper[x == 0] <- 1 - ( (1 - level) / 2) ^ ( 1 / n )[x == 0]
    }
  } else if (dist == "poi") {
    # Jeffreys prior for Poisson rate uses gamma distribution,
    # as defined in Li et al.
    CI.lower <- qgamma(alpha / 2, x + (ai - cc), scale = 1 / n)
    CI.lower[x == 0] <-  0
    CI.upper <- qgamma(1 - alpha / 2, (x + (ai + cc)), scale = 1 / n)
  }
  CI <- cbind(Lower = CI.lower, Upper = CI.upper)
  CI
}

quadroot <- function(a, b, c) {
	# GET ROOTS OF A QUADRATIC EQUATION
	r1x <- ( - b + sqrt(b ^ 2 - 4 * a * c) ) / (2 * a)
	r2x <- ( - b - sqrt(b ^ 2 - 4 * a * c) ) / (2 * a)
	r1 <- pmin(r1x, r2x)
	r2 <- pmax(r1x, r2x)
	cbind(r1, r2)
}

