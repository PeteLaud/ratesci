
# Internal function for Clopper-Pearson/Garwood and mid-p, binomial or Poisson
exactci <- function( # function to calculate exact 'exact' confidence interval for a single
                    # binomial or Poisson rate x/n
                    x,
                    n,
                    level = 0.95,
                    midp = TRUE,
                    distrib = "bin",
                    precis = 8) {
  alpha <- 1 - level
  if (as.character(midp) == "TRUE") midp <- 0.5
  if (distrib == "bin") {
    lowroot <- function(p) {
      pbinom(x - 1, n, p) + midp * dbinom(x, n, p) -
        (1 - alpha / 2)
    }
    uproot <- function(p) pbinom(x, n, p) - midp * dbinom(x, n, p) - alpha / 2
  } else if (distrib == "poi") {
    lowroot <- function(p) ppois(x, p) + midp * dpois(x, p) - (1 - alpha / 2)
    uproot <- function(p) ppois(x, p) - midp * dpois(x, p) - alpha / 2
  }
  lower <- bisect(
    ftn = lowroot, precis = precis, uplow = "low",
    contrast = "p", distrib = distrib
  )
  upper <- bisect(
    ftn = uproot, precis = precis, uplow = "up", contrast = "p",
    distrib = distrib
  )
  return(cbind(Lower = lower, Upper = upper) / ifelse(distrib == "poi", n, 1))
}


#' Skewness-corrected asymptotic score ("SCAS") confidence intervals for
#' single binomial or Poisson rate using closed-form calculations.
#' This function is vectorised in x, n.
#'
#' @param x Numeric vector of number of events.
#' @param n Numeric vector of sample sizes (for binomial rates) or exposure
#'   times (for Poisson rates).
#' @inheritParams scoreci
#' @export
scaspci <- function(x,
                    n,
                    distrib = "bin",
                    level = 0.95,
                    cc = FALSE,
                    ...) {
#  x <- Rmpfr::mpfr(x, 120)
#  n <- Rmpfr::mpfr(n, 120)
#  level <- Rmpfr::mpfr(level, 120)
  if (as.character(cc) == "TRUE") cc <- 0.5
  z <- qnorm(1 - (1 - level) / 2)
  if (distrib == "poi") {
    Du <- (x + cc) / n - (z^2 - 1) / (6 * n)
#    Dl <- Rmpfr::pmax(0, (x - cc) / n - (z^2 - 1) / (6 * n))
    Dl <- pmax(0, (x - cc) / n - (z^2 - 1) / (6 * n))
    A <- 1
    Bu <- -2 * Du - z^2 / n
    Bl <- -2 * Dl - z^2 / n
    Cu <- Du^2
    Cl <- Dl^2
    D0 <- -1 / (6 * n) - x / n
    B0 <- 2 * D0
    A0 <- 1
    C0 <- D0^2
  } else if (distrib == "bin") {
    E <- (z^2 - 1) / (3 * n) - 1
    # Alteration to published formula,
    # to deal with non-nested intervals when level > 0.99
#    Du <- Rmpfr::pmax(0, (n - x - cc) / n - (z^2 - 1) / (6 * n))
#    Dl <- Rmpfr::pmax(0, (x - cc) / n - (z^2 - 1) / (6 * n))
    Du <- pmax(0, (n - x - cc) / n - (z^2 - 1) / (6 * n))
    Dl <- pmax(0, (x - cc) / n - (z^2 - 1) / (6 * n))
    A <- z^2 / n + E^2
    Bu <- 2 * E * Du - z^2 / n
    Bl <- 2 * E * Dl - z^2 / n
    Cu <- Du^2
    Cl <- Dl^2
    E0 <- 1 + 1 / (3 * n)
    D0 <- -1 / (6 * n) - x / n
    A0 <- E0^2
    B0 <- 2 * E0 * D0
    C0 <- D0^2
  }

  CI <- (cbind(
#    Lower = Rmpfr::asNumeric((-Bl - sqrt(Rmpfr::pmax(0, Bl^2 - 4 * A * Cl))) / (2 * A)),
#    MLE = Rmpfr::asNumeric((-B0 - sqrt(Rmpfr::pmax(0, (B0^2 - 4 * A0 * C0)))) / (2 * A0)),
    Lower = ((-Bl - sqrt(pmax(0, Bl^2 - 4 * A * Cl))) / (2 * A)),
    MLE = ((-B0 - sqrt(pmax(0, (B0^2 - 4 * A0 * C0)))) / (2 * A0)),
    Upper = if (distrib == "bin") {
#      Rmpfr::asNumeric(1 - (-Bu - sqrt(Rmpfr::pmax(0, Bu^2 - 4 * A * Cu))) / (2 * A))
      (1 - (-Bu - sqrt(pmax(0, Bu^2 - 4 * A * Cu))) / (2 * A))
    } else {
#      Rmpfr::asNumeric((-Bu + sqrt(Rmpfr::pmax(0, Bu^2 - 4 * A * Cu))) / (2 * A))
      ((-Bu + sqrt(pmax(0, Bu^2 - 4 * A * Cu))) / (2 * A))
    }
  ))
  return((CI))
}

#' Selected confidence intervals for the single binomial or Poisson rate.
#'
#' Confidence intervals for the single binomial or Poisson rate. Including
#' SCAS or Jeffreys intervals, with or without continuity correction, and
#' 'exact' Clopper-Pearson/Garwood or mid-p intervals.
#'
#' @param x Numeric vector of number of events.
#' @param n Numeric vector of sample size (for binomial rate) or exposure
#'   times (for Poisson rate).
#' @param distrib Character string indicating distribution assumed for the input
#'   data: "bin" = binomial (default), "poi" = Poisson.
#' @param level Number specifying confidence level (between 0 and 1, default
#'   0.95).
#' @param cc Number or logical (default FALSE) specifying continuity
#'   correction.
#' @param precis Number (default 6) specifying precision (i.e. number of decimal
#'   places) to be used in optimisation subroutine for the SCAS and exact methods.
#' @return A list containing, for each method, a matrix containing lower and upper
#'   confidence limits for each value of x and n. Methods shown depend on the cc
#'   parameter, which specifies whether the continuity correction is applied to
#'   the SCAS and Jeffreys methods. The corresponding 'exact' method is
#'   Clopper-Pearson/Garwood if cc == TRUE and mid-p if cc == FALSE.
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#' @export
rateci <- function(x,
                   n,
                   distrib = "bin",
                   level = 0.95,
                   cc = FALSE,
                   precis = 6) {
  # in case x is input as a vector but n is not
  if (length(n) < length(x) && length(x) > 1) {
    n <- rep(n, length.out = length(x))
  }
  if (as.character(cc) == "TRUE") cc <- 0.5

  ci_scas <- scaspci(
    x = x,
    n = n,
    distrib = distrib,
    level = level,
    cc = cc
  )[, c(1, 3)]
  ci_jeff <- jeffreysci(
    x = x,
    n = n,
    ai = 0.5,
    bi = 0.5,
    cc = cc,
    level = level,
    distrib = distrib,
    adj = TRUE
  )[, c(1, 2)]
  ci_exact <- exactci(
    x = x,
    n = n,
    level = level,
    midp = 0.5 - cc,
    distrib = distrib,
    precis = precis
  )
  if (cc == 0) {
    return(list(scas = ci_scas, jeff = ci_jeff, midp = ci_exact))
  } else if (cc == 0.5) {
    if (distrib == "bin") {
      return(list(scas_cc = ci_scas, jeff_cc = ci_jeff, cp = ci_exact))
    } else {
      return(list(scas_cc = ci_scas, jeff_cc = ci_jeff, garwood = ci_exact))
    }
  } else {
    return(list(scas_cc = ci_scas, jeff_cc = ci_jeff))
    # exact method not applicable if using a compromise value of cc
  }
}
