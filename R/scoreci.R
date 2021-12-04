#' Score confidence intervals for comparisons of independent binomial or Poisson
#' rates.
#'
#' Score-based confidence intervals for the rate (or risk) difference ("RD") or
#' ratio ("RR") for independent binomial or Poisson rates, or for odds ratio
#' ("OR", binomial only). Including options for bias correction (from Miettinen
#' & Nurminen), skewness correction ("GNbc" method from Laud & Dane, developed
#' from Gart & Nam, and generalised as "SCAS" in Laud 2017) and continuity
#' correction (for strictly conservative coverage).
#' Also includes score intervals for a single binomial proportion or Poisson
#' rate. Based on the Wilson score interval, when corrected for skewness,
#' coverage is almost identical to the mid-p method, or Clopper-Pearson
#' when also continuity-corrected.
#' Hypothesis tests for superiority or non-inferiority are provided using the
#' same score, to ensure consistency between test and CI.
#' This function is vectorised in x1, x2, n1, and n2.  Vector inputs may also be
#' combined into a single stratified analysis (e.g. meta-analysis), either using
#' fixed effects, or the more general random effects "TDAS" method, which
#' incorporates stratum variability using a t-distribution score (inspired by
#' Hartung-Knapp-Sidik-Jonkman).
#' For fixed-effects analysis of stratified datasets, with weighting = "MH" for
#' RD or RR, or weighting = "IVS" for OR, omitting the skewness correction
#' produces the CMH test, together with a coherent confidence interval for the
#' required contrast.
#'
#' @param x1,x2 Numeric vectors of numbers of events in group 1 & group 2
#'   respectively.
#' @param n1,n2 Numeric vectors of sample sizes (for binomial rates) or exposure
#'   times (for Poisson rates) in each group.
#' @param distrib Character string indicating distribution assumed for the input
#'   data: "bin" = binomial (default), "poi" = Poisson.
#' @param contrast Character string indicating the contrast of interest: "RD" =
#'   rate difference (default), "RR" = rate ratio, "OR" = odds ratio.
#'   contrast="p" gives an interval for the single proportion or rate x1/n1.
#' @param level Number specifying confidence level (between 0 and 1, default
#'   0.95).
#' @param skew Logical (default TRUE) indicating whether to apply skewness
#'   correction (for the SCAS method recommended in Laud 2017) or not (for
#'   the Miettinen-Nurminen method).
#' @param simpleskew Logical (default FALSE) indicating whether to use the
#'   "simplified" skewness correction instead of the quadratic solution.
#'   See Laud 2021 for details. NOTE: this version of the score is only
#'   suitable for obtaining confidence limits, not p-values.
#' @param ORbias Logical (default is TRUE for contrast="OR", otherwise
#'   NULL) indicating whether to apply additional bias correction for OR derived
#'   from Gart 1985. (Corrigendum to Laud 2017, published May 2018).
#'   Only applies if contrast is "OR".
#' @param bcf Logical (default TRUE) indicating whether to apply bias correction
#'   in the score denominator. Applicable to distrib = "bin" only. (NB: bcf =
#'   FALSE option is really only included for legacy validation against previous
#'   published methods (i.e. Gart & Nam, Mee, or standard Chi-squared test).
#'   Ignored for contrast = "p".
#' @param cc Number or logical (default FALSE) specifying (amount of) continuity
#'   correction. Numeric value is taken as the gamma parameter in Laud 2017,
#'   Appendix S2 (default 0.5 if cc = TRUE).
#'   IMPORTANT NOTES:
#'   1) This is a 'continuity correction' aimed at approximating strictly
#'   conservative coverage, NOT for dealing with zero cell counts. Such
#'   'sparse data adjustments' are not needed in the score method,
#'   except to deal with double-zero cells for RD (& double-100% cells for
#'   binomial RD & RR) with IVS/INV weights.
#'   2) The continuity corrections provided here have not been fully tested for
#'   stratified methods.
#' @param theta0 Number to be used in a one-sided significance test (e.g.
#'   non-inferiority margin). 1-sided p-value will be <0.025 iff 2-sided 95\% CI
#'   excludes theta0. If bcf = FALSE and skew = FALSE this gives a
#'   Farrington-Manning test.
#'   By default, a two-sided test against theta0 = 0 (for RD) or 1 (for RR/OR)
#'   is also output: if bcf = FALSE and skew = FALSE this is the same as
#'   Pearson's Chi-squared test.
#' @param precis Number (default 6) specifying precision (i.e. number of decimal
#'   places) to be used in optimisation subroutine for the confidence interval.
#' @param plot Logical (default FALSE) indicating whether to output plot of the
#'   score function
#' @param plotmax Numeric value indicating maximum value to be displayed on
#'   x-axis of plots (useful for ratio contrasts which can be infinite).
#' @param xlim pair of values indicating range of values to be plotted.
#' @param ylim pair of values indicating range of values to be plotted.
#' @param stratified Logical (default FALSE) indicating whether to combine
#'   vector inputs into a single stratified analysis.
#'   IMPORTANT NOTE: The mechanism for stratified calculations is enabled for
#'   contrast = "p", but the performance of the resulting intervals has not
#'   been fully evaluated.
#' @param weighting String indicating which weighting method to use if
#'   stratified = "TRUE":
#'   "IVS" = Inverse Variance of Score (see Laud 2017 for details),
#'   "INV" = Inverse Variance (bcf omitted, default for contrast = "OR"),
#'   "MH" = Mantel-Haenszel (default for contrast = "RD" or "RR"),
#'   "MN" = Miettinen-Nurminen iterative weights.
#'   For CI consistent with a CMH test, select skew = FALSE and use
#'   MH weighting for RD/RR and IVS for OR.
#' @param wt Numeric vector containing (optional) user-specified weights.
#' @param sda Sparse data adjustment to avoid zero variance when x1 + x2 = 0:
#'           Only applied when stratified = TRUE.
#'           Default 0.5 for RD with IVS/INV weights.
#'           Not required for RR/OR, default is to remove double-zero strata
#'           instead.
#' @param fda Full data adjustment to avoid zero variance when x1 + x2 = n1 + n2:
#'           Only applied when stratified = TRUE.
#'           Default 0.5 for RD & RR with IVS/INV weights.
#'           Not required for OR, default is to remove affected strata.
#' @param dropzeros Logical (default FALSE) indicating whether to drop
#'   uninformative strata for RR/OR, even when the choice of weights would allow
#'   them to be retained for a fixed effects analysis.
#'   Has no effect on estimates, just the heterogeneity test.
#' @param RRtang Logical indicating whether to use Tang's score for RR:
#'   Stheta = (p1hat - p2hat * theta) / p2d (see Tang 2020).
#'   Default TRUE for stratified = TRUE, with weighting = "IVS" or "INV".
#'   Forced to FALSE for stratified = TRUE, with fixed weighting.
#'   Experimental for distrib = "poi".
#' @param hetplot Logical (default FALSE) indicating whether to output plots for
#'   evaluating heterogeneity of stratified datasets.
#' @param MNtol Numeric value indicating convergence tolerance to be used in
#'   iteration with weighting = "MN".
#' @param tdas (deprecated: parameter renamed to random)
#' @param random Logical (default FALSE) indicating whether to perform random
#'   effects meta-analysis for stratified data, using the t-distribution (TDAS)
#'   method for stratified data (defined in Laud 2017).
#'   NOTE: If random = TRUE, then skew = TRUE only affects the per-stratum
#'   estimates.
#' @param prediction Logical (default FALSE) indicating whether to produce
#'   a prediction interval (work in progress).
#' @param warn Logical (default TRUE) giving the option to suppress warnings.
#' @param ... Other arguments.
#' @importFrom stats pchisq pf pnorm pt qbeta qgamma qnorm qqnorm qt dbinom
#' @importFrom graphics abline lines text
#' @return A list containing the following components: \describe{
#'   \item{estimates}{a matrix containing estimates of the rates in each group
#'   and of the requested contrast, with its confidence interval} \item{pval}{a
#'   matrix containing details of the corresponding 2-sided significance test
#'   against the null hypothesis that p_1 = p_2, and one-sided significance
#'   tests against the null hypothesis that theta >= or <= theta0}
#'   \item{call}{details of the function call} } If stratified = TRUE, the
#'   following outputs are added: \describe{ \item{Qtest}{a vector of values
#'   describing and testing heterogeneity} \item{weighting}{a string indicating
#'   the selected weighting method} \item{stratdata}{a matrix containing stratum
#'   estimates and weights}}
#' @examples
#' # Binomial RD, SCAS method:
#' scoreci(
#'   x1 = c(12, 19, 5), n1 = c(16, 29, 56),
#'   x2 = c(1, 22, 0), n2 = c(16, 30, 29)
#' )
#'
#' # Binomial RD, MN method:
#' scoreci(
#'   x1 = c(12, 19, 5), n1 = c(16, 29, 56),
#'   x2 = c(1, 22, 0), n2 = c(16, 30, 29), skew = FALSE
#' )
#'
#' # Poisson RR, SCAS method:
#' scoreci(x1 = 5, n1 = 56, x2 = 0, n2 = 29, distrib = "poi", contrast = "RR")
#'
#' # Poisson RR, MN method:
#' scoreci(
#'   x1 = 5, n1 = 56, x2 = 0, n2 = 29, distrib = "poi",
#'   contrast = "RR", skew = FALSE
#' )
#'
#' # Binomial rate, SCAS method:
#' scoreci(x1 = c(5, 0), n1 = c(56, 29), contrast = "p")
#'
#' # Binomial rate, Wilson score method:
#' scoreci(x1 = c(5, 0), n1 = c(56, 29), contrast = "p", skew = FALSE)
#'
#' # Poisson rate, SCAS method:
#' scoreci(x1 = c(5, 0), n1 = c(56, 29), distrib = "poi", contrast = "p")
#'
#' # Stratified example, using data from Hartung & Knapp:
#' scoreci(
#'   x1 = c(15, 12, 29, 42, 14, 44, 14, 29, 10, 17, 38, 19, 21),
#'   x2 = c(9, 1, 18, 31, 6, 17, 7, 23, 3, 6, 12, 22, 19),
#'   n1 = c(16, 16, 34, 56, 22, 54, 17, 58, 14, 26, 44, 29, 38),
#'   n2 = c(16, 16, 34, 56, 22, 55, 15, 58, 15, 27, 45, 30, 38),
#'   stratified = TRUE
#' )
#'
#' # TDAS example, using data from Hartung & Knapp:
#' scoreci(
#'   x1 = c(15, 12, 29, 42, 14, 44, 14, 29, 10, 17, 38, 19, 21),
#'   x2 = c(9, 1, 18, 31, 6, 17, 7, 23, 3, 6, 12, 22, 19),
#'   n1 = c(16, 16, 34, 56, 22, 54, 17, 58, 14, 26, 44, 29, 38),
#'   n2 = c(16, 16, 34, 56, 22, 55, 15, 58, 15, 27, 45, 30, 38),
#'   stratified = TRUE, random = TRUE
#' )
#'
#' # Stratified example, with extremely rare instance of non-calculable skewness
#' # correction seen on plot of score function:
#' scoreci(
#'   x1 = c(1, 16), n1 = c(20, 40), x2 = c(0, 139), n2 = c(80, 160),
#'   contrast = "RD", skew = TRUE, simpleskew = FALSE,
#'   distrib = "bin", stratified = TRUE, plot = TRUE, weighting = "IVS"
#' )
#'
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#' @references Laud PJ. Equal-tailed confidence intervals for comparison of
#'   rates. Pharmaceutical Statistics 2017; 16:334-348.
#'
#'   Laud PJ. Corrigendum: Equal-tailed confidence intervals for comparison of
#'   rates. Pharmaceutical Statistics 2018; 17:290-293.
#'
#'   Laud PJ, Dane A. Confidence intervals for the difference between
#'   independent binomial proportions: comparison using a graphical approach and
#'   moving averages. Pharmaceutical Statistics 2014; 13(5):294–308.
#'
#'   Miettinen OS, Nurminen M. Comparative analysis of two rates. Statistics in
#'   Medicine 1985; 4:213-226.
#'
#'   Farrington CP, Manning G. Test statistics and sample size formulae for
#'   comparative binomial trials with null hypothesis of non-zero risk
#'   difference or non-unity relative risk. Statistics in Medicine 1990;
#'   9(12):1447-1454.
#'
#'   Gart JJ. Analysis of the common odds ratio: corrections for bias and
#'   skewness. Bulletin of the International Statistical Institute 1985,
#'   45th session, book 1, 175-176.
#'
#'   Gart JJ, Nam Jm. Approximate interval estimation of the ratio of binomial
#'   parameters: a review and corrections for skewness. Biometrics 1988;
#'   44(2):323-338.
#'
#'   Gart JJ, Nam Jm. Approximate interval estimation of the difference in
#'   binomial parameters: correction for skewness and extension to multiple
#'   tables. Biometrics 1990; 46(3):637-643.
#'
#'   Tang Y. Score confidence intervals and sample sizes for stratified
#'   comparisons of binomial proportions. Statistics in Medicine 2020;
#'   39:3427–3457.
#'
#' @export
scoreci <- function(x1,
                    n1,
                    x2 = NULL,
                    n2 = NULL,
                    distrib = "bin",
                    contrast = "RD",
                    level = 0.95,
                    skew = TRUE,
                    simpleskew = FALSE,
                    ORbias = TRUE,
                    RRtang = NULL,
                    bcf = TRUE,
                    cc = FALSE,
                    theta0 = NULL,
                    precis = 6,
                    plot = FALSE,
                    plotmax = 100,
                    hetplot = FALSE,
                    xlim = NULL,
                    ylim = NULL,
                    stratified = FALSE,
                    weighting = NULL,
                    MNtol = 1E-8,
                    wt = NULL,
                    sda = NULL,
                    fda = NULL,
                    dropzeros = FALSE,
                    tdas = NULL,
                    random = FALSE,
                    prediction = FALSE,
                    warn = TRUE,
                    ...) {
  if (!is.null(tdas)) {
    warning(
      "argument tdas is deprecated; please use random instead.",
      call. = FALSE
    )
    random <- tdas
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
    print("argument x2 or n2 missing")
    stop()
  }
  if ((contrast != "p") && length(x1) != length(x2)) {
    print("x1 and x2 must be the same length")
    stop()
  }
  if ((contrast != "p") && !is.numeric(c(x1, n1, x2, n2, theta0))) {
    print("Non-numeric inputs!")
    stop()
  }
  if ((contrast == "p") && !is.numeric(c(x1, n1, theta0))) {
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
  if (contrast != "OR") {
    ORbias <- NULL
  }
  if (stratified == TRUE && is.null(weighting)) {
    weighting <- switch(contrast,
                        "RD" = "MH",
                        "RR" = "MH",
                        "OR" = "INV",
                        "p" = "IVS"
    )
  }
  if (is.null(sda)) {
    if (stratified == TRUE && weighting %in% c("IVS", "INV") &&
        contrast %in% c("RD")) {
      sda <- 0.5
    } else {
      sda <- 0
    }
  }
  if (is.null(fda)) {
    if (stratified == TRUE && weighting %in% c("IVS", "INV") &&
        contrast %in% c("RD", "RR")) {
      fda <- 0.5
    } else {
      fda <- 0
    }
  }
  if (!is.null(theta0)) {
    if (contrast == "RD") {
      if (distrib == "bin" && (theta0 < -1 || theta0 > 1)) {
        print("Impossible theta0!")
        stop()
      }
    } else if (contrast == "p") {
      if (theta0 < 0 || (distrib == "bin" && theta0 > 1)) {
        print("Impossible theta0!")
        stop()
      }
    }
  }
  if (distrib == "poi" && contrast == "OR") {
    print("Odds ratio not applicable to Poisson rates")
    stop()
  }
  if (contrast == "p" && weighting == "MN" && stratified == TRUE) {
    print("MN weights not applicable to the single rate")
    stop()
  }
  if (cc != FALSE && stratified == TRUE && contrast != "OR") {
    if (warn == TRUE) {
      print(paste(
        "Warning: Continuity correction is experimental for",
        "stratified RD and RR"
      ))
    }
  }
  # Tang RR score intended only for IVS/INV weighting -
  # Tang p3431 does not use it for MH weights.
  if (contrast != "RR" && !is.null(RRtang) && warn == TRUE) {
    print(paste(
      "Warning: RRtang argument has no effect for contrast =", contrast
    ))
    RRtang = FALSE
  } else if (contrast == "RR") {
    if (stratified == TRUE && !(weighting %in% c("IVS", "INV"))) {
      if (!is.null(RRtang)) {
      if (warn == TRUE && RRtang == TRUE) {
        print(paste(
          "Warning: RRtang set to FALSE - option designed for inverse variance weighting only"
        ))
      }
      }
    RRtang <- FALSE
    } else if (is.null(RRtang)) {
      RRtang <- TRUE
    }
  } else {
    RRtang <- FALSE
  }

  if (as.character(cc) == "TRUE") cc <- 0.5

  nstrat <- length(x1)
  # in case x1,x2 are input as vectors but n1,n2 are not
  if (length(n1) < nstrat && nstrat > 1) n1 <- rep(n1, length.out = nstrat)
  if (contrast != "p" && length(n2) < nstrat && nstrat > 1) {
    n2 <- rep(n2, length.out = nstrat)
  }

  # check for empty strata, and for uninformative strata
  if (stratified == TRUE) {
    # for stratified calculations, remove any strata with no observations
    if (contrast == "p") {
      empty_strat <- (n1 == 0)
    } else {
      empty_strat <- (n1 == 0 | n2 == 0) |
        # Also drop uninformative strata for RR & OR
        # Note for RR & OR with IVS/INV weighting such strata have zero weight
        # and for RR with fixed weights and RRtang = FALSE they dont contribute
        # for a fixed effects analysis,
        # so could remain to avoid ITT concerns - but note this has
        # implications for heterogeneity test and random effects method.
        ((x1 == 0 & x2 == 0) & !(sda > 0) & contrast == "RR" &
          weighting %in% c("INV", "IVS") & RRtang == FALSE) |
        #       below not needed because RRtang is forced to FALSE anyway
        #        ((x1 == 0 & x2 == 0) & !(sda > 0) & contrast == "RR" &
        #          !(weighting %in% c("INV", "IVS")) & RRtang == TRUE) |
        ((x1 == 0 & x2 == 0) & !(sda > 0) & contrast == "OR" &
          !(weighting %in% c("INV", "IVS"))) |
        ((x1 == n1 & x2 == n2) & !(fda > 0) & contrast == "OR" &
          !(weighting %in% c("INV", "IVS")))

      if (random == TRUE || dropzeros == TRUE) {
        # Exclude more uninformative strata using dropzeros argument,
        # to prevent them affecting the heterogeneity test.
        # Similarly if applying the TDAS method with 'random' argument.
        # This might be unnecessary for some random effects analyses, such as
        # RR with RRtang = TRUE and INV/IVS weighting and OR with INV/IVS
        # weighting, but needs further research.
        empty_strat <- (n1 == 0 | n2 == 0) |
          ((x1 == 0 & x2 == 0) & !(sda > 0) & (contrast %in% c("RR", "OR"))) |
          ((x1 == n1 & x2 == n2) & !(fda > 0) & (contrast == "OR"))
      }
    }

    x1 <- x1[!empty_strat]
    n1 <- n1[!empty_strat]
    if (contrast != "p") {
      x2 <- x2[!empty_strat]
      n2 <- n2[!empty_strat]
    }
    if (warn == TRUE && sum(empty_strat) > 0) {
      print(paste(
        "Note: at least one stratum contributed no information",
        "and was removed"
      ))
    }

    nstrat <- length(x1) # update nstrat after removing empty strata
    sda <- rep_len(sda, length.out = nstrat)
    fda <- rep_len(fda, length.out = nstrat)

    # Sparse/Full data adjustment:
    # for double-zero cells for RD, add sda (default 0.5) to avoid 100% weight
    # with IVS weights. Add fda (default 0.5) for double-100% cells for binomial RD
    # & RR (yes, really - variance becomes zero at theta=1).
    # NB this should only be necessary for weighting = 'IVS' or 'INV',
    # (but some might prefer to specify instead of excluding empty strata.)
    # Default is sda=fda=0.5 for RD, fda=0.5 for RR, sda=0 for OR
    if (stratified == TRUE) {
      # && weighting %in% c("IVS", "INV")) {
      #  if (contrast %in% c("OR") && weighting == "MH") {
      # ??OR struggles with zero event counts with MH weights
      #    zero_rd <- (x1 == 0 | x2 == 0 | x1 == n1 | x2 == n2) #& !empty_strat
      #  } else zero_rd <- 0
      # NB users have option to apply sda for RR/OR instead of dropping strata
      if (contrast %in% c("RD", "RR", "OR")) {
        zero_rd <- (x1 == 0 & x2 == 0) # & !empty_strat
      } else {
        zero_rd <- 0
      }
      x1 <- x1 + (zero_rd * sda)
      x2 <- x2 + (zero_rd * sda)
      n1 <- n1 + (zero_rd * 2 * sda)
      n2 <- n2 + (zero_rd * 2 * sda)
      if (contrast %in% c("RD", "RR", "OR") & distrib == "bin") {
        full_rd <- (x1 == n1 & x2 == n2) # & !empty_strat
      } else {
        full_rd <- 0
      }
      x1 <- x1 + (full_rd * fda)
      x2 <- x2 + (full_rd * fda)
      n1 <- n1 + (full_rd * 2 * fda)
      n2 <- n2 + (full_rd * 2 * fda)
    }
  }

  # Warnings if removal of empty strata leave 0 or 1 stratum
  if (stratified == TRUE && nstrat == 1) {
    stratified <- FALSE
    random <- FALSE
    if (warn == TRUE) print("Note: only one stratum!")
  }
  if (nstrat == 0) {
    if (warn == TRUE) print("Warning: no usable data!")
    stratified <- FALSE
    if (contrast %in% c("RR", "OR")) {
      x1 <- x2 <- 0
      n1 <- n2 <- 10
    }
  }

  p1hat <- x1 / n1
  p2hat <- x2 / n2

  # wrapper function for scoretheta
  myfun <- function(theta,
                    randswitch = random,
                    ccswitch = cc,
                    stratswitch = stratified,
                    skewswitch = skew,
                    ssswitch = simpleskew,
                    lev = level,
                    predswitch = FALSE) {
    scoretheta(
      theta = theta, x1 = x1, x2 = x2, n1 = n1, n2 = n2, bcf = bcf,
      contrast = contrast, distrib = distrib, stratified = stratswitch,
      wt = wt, weighting = weighting, MNtol = MNtol,
      random = randswitch, prediction = predswitch, skew = skewswitch,
      simpleskew = ssswitch, ORbias = ORbias, RRtang = RRtang,
      cc = ccswitch, level = lev
    )$score
  }

  # wrapper function for plotting discriminants
  mydsct <- function(theta,
                     randswitch = random,
                     ccswitch = cc,
                     stratswitch = stratified,
                     predswitch = FALSE) {
    myout <- scoretheta(
      theta = theta, x1 = x1, x2 = x2, n1 = n1, n2 = n2,
      bcf = bcf, contrast = contrast, distrib = distrib,
      stratified = stratswitch, wt = wt, weighting = weighting,
      MNtol = MNtol, random = randswitch, RRtang = RRtang,
      prediction = predswitch, skew = skew, ORbias = ORbias,
      cc = ccswitch, simpleskew = simpleskew, level = level
    )$dsct
  }

  # find point estimate for theta as the value where scoretheta = 0
  # fixed effects point estimate taken with no cc
  point_FE <- bisect(
    ftn = function(theta) {
      myfun(theta, randswitch = FALSE, ccswitch = 0, lev = 0) - 0
    }, contrast = contrast, distrib = distrib,
    precis = precis + 1, uplow = "low"
  )
  point_FE[myfun(theta = 1, randswitch = random, ccswitch = 0, lev = 0) == 0] <- 1
  point <- point_FE

  # random effects point estimate if required
  if (stratified == TRUE) {
    point <- bisect(
      ftn = function(theta) {
        myfun(theta, randswitch = random, ccswitch = 0, lev = 0) - 0
      }, contrast = contrast, distrib = distrib,
      precis = precis + 1, uplow = "low"
    )
  }
  point[myfun(theta = 1, randswitch = random, ccswitch = 0, lev = 0) == 0] <- 1

  # fix some extreme cases with zero counts
  if (stratified == TRUE) {
    if (skew == FALSE & contrast %in% c("RR", "OR")) {
      point[sum(x2) > 0 & sum(x1) == 0] <- 0
    }
    if (skew == FALSE & contrast %in% c("RR", "OR")) {
      point[sum(x1) > 0 & sum(x2) == 0] <- Inf
    }
  } else {
    if (contrast %in% c("RR", "OR")) {
      point[x1 == 0 & x2 == 0] <- NA
      point[(x1 > 0 & x2 == 0) & skew == FALSE] <- Inf
      if (distrib == "bin") point[x1 == n1 & x2 == n2] <- 1
    }
    if (contrast == "OR") {
      point[x1 == n1 & x2 == n2] <- NA
      point[(x1 == n1 & x2 > 0) & skew == FALSE] <- Inf
    }
    if (contrast == "RD") {
      if (distrib == "bin") point[x1 == n1 & x2 == n2] <- 0
      point[x1 == 0 & x2 == 0] <- 0
    }
  }

  # identify certain quantities evaluated at the maximum likelihood point
  # estimate for theta (incorporating random effects)
  at_MLE <- scoretheta(
    theta = ifelse(is.na(point), 1, point), x1 = x1, x2 = x2,
    n1 = n1, n2 = n2,
    bcf = bcf, contrast = contrast,
    distrib = distrib, stratified = stratified,
    weighting = weighting, MNtol = MNtol, wt = wt,
    random = random, skew = skew, ORbias = ORbias,
    RRtang = RRtang, cc = cc, simpleskew = simpleskew, level = level
  )
  ##  Stheta_MLE <- at_MLE$Stheta
  p1d_MLE <- at_MLE$p1d
  p2d_MLE <- at_MLE$p2d
  wt_MLE <- at_MLE$wt
  ##  V_MLE <- at_MLE$V

  # if stratified = TRUE, options are available for assuming fixed effects
  # (random = FALSE) or random effects (random = TRUE). The IVS weights are
  # different for each version, which in turn can lead to a different point
  # estimate, at which certain quantities are evaluated. However, fixed effects
  # estimates are needed for both, in particular for the heterogeneity test
  if (stratified == TRUE && nstrat > 1) {
    at_FE <- scoretheta(
      theta = point_FE, x1 = x1, x2 = x2, n1 = n1, n2 = n2,
      bcf = bcf, contrast = contrast,
      distrib = distrib, stratified = stratified,
      weighting = weighting, MNtol = MNtol, wt = wt,
      random = FALSE, skew = skew, ORbias = ORbias,
      RRtang = RRtang, cc = cc, simpleskew = simpleskew, level = level
    )
    Stheta_FE <- at_FE$Stheta
    wt_FE <- at_FE$wt
    V_FE <- at_FE$V
    tau2_FE <- at_FE$tau2
    Q_each <- at_FE$Q_j
    Q_FE <- at_FE$Q
    I2 <- max(0, 100 * (Q_FE - (nstrat - 1)) / Q_FE)
    pval_het <- 1 - pchisq(Q_FE, nstrat - 1)
    # Qualitative interaction test is calculated further down

    # as per M&N p218 (little r), actually no longer needed for point estimate.
    p1hat_w <- sum(wt_MLE * x1 / n1) / sum(wt_MLE)
    p2hat_w <- sum(wt_MLE * x2 / n2) / sum(wt_MLE)
    # as per M&N p218 (big R)
    p1d_w <- sum(wt_MLE * at_MLE$p1d) / sum(wt_MLE)
    if (contrast != "p") {
      p2d_w <- sum(wt_MLE * at_MLE$p2d) / sum(wt_MLE)
    } else {
      p2d_w <- NULL
    }
    if (!is.null(wt)) {
      weighting <- "User-defined"
    }
  } else {
    # if removing empty strata leaves only one stratum treat as unstratified
    p1hat_w <- p1hat
    p2hat_w <- p2hat
    p1d_w <- p1d_MLE
    if (contrast != "p") p2d_w <- p2d_MLE else p2d_w <- NULL
    wt_MLE <- NULL
    ## Sdot <- Q_each <- NULL
  }

  # z- or t- quantile required for specified significance level
  if (stratified == TRUE && random == TRUE) {
    qtnorm <- qt(1 - (1 - level) / 2, nstrat - 1) # for t-distribution method
  } else {
    qtnorm <- qnorm(1 - (1 - level) / 2)
  }

  # Use bisection routine to locate lower and upper confidence limits
  lower <- bisect(
    ftn = function(theta) {
      myfun(theta) - qtnorm
    }, contrast = contrast, distrib = distrib,
    precis = precis + 1, uplow = "low"
  )
  upper <- bisect(
    ftn = function(theta) {
      myfun(theta) + qtnorm
    }, contrast = contrast, distrib = distrib,
    precis = precis + 1, uplow = "up"
  )

  # Produce prediction interval if required
  if (stratified == TRUE && random == TRUE && prediction == TRUE) {
    qtnorm <- qt(1 - (1 - level) / 2, nstrat - 2)
    lowpred <- bisect(
      ftn = function(theta) {
        myfun(theta, predswitch = TRUE) - qtnorm
      }, contrast = contrast, distrib = distrib,
      precis = precis + 1, uplow = "low"
    )
    uppred <- bisect(
      ftn = function(theta) {
        myfun(theta, predswitch = TRUE) + qtnorm
      }, contrast = contrast, distrib = distrib,
      precis = precis + 1, uplow = "up"
    )
    pred <- cbind(round(cbind(Lower = lowpred, Upper = uppred), precis))
  } else {
    pred <- NULL
  }

  # get estimate & CI for each stratum
  if (stratified == TRUE) {
    at_MLE_unstrat <- scoretheta(
      theta = point, x1 = x1, x2 = x2, n1 = n1,
      n2 = n2, bcf = bcf, contrast = contrast,
      distrib = distrib, stratified = FALSE,
      weighting = weighting, MNtol = MNtol, wt = wt,
      random = random, skew = skew, ORbias = ORbias,
      RRtang = RRtang, cc = cc,
      simpleskew = simpleskew, level = level
    )
    point_FE_unstrat <- bisect(
      ftn = function(theta) {
        myfun(theta, randswitch = FALSE, ccswitch = 0,
              stratswitch = FALSE, lev = 0) - 0
      }, contrast = contrast, distrib = distrib,
      precis = precis + 1, uplow = "low"
    )
    if (contrast == "RD") {
      if (distrib == "bin") point_FE_unstrat[x1 == n1 & x2 == n2] <- 0
      point_FE_unstrat[x1 == 0 & x2 == 0] <- 0
    }
    if (contrast == "RR") {
      point_FE_unstrat[x1 == 0 & x2 == 0] <- 1
    }
    lower_unstrat <- bisect(
      ftn = function(theta) {
        myfun(theta, stratswitch = FALSE) - qnorm(1 - (1 - level) / 2)
      }, contrast = contrast, distrib = distrib,
      precis = precis + 1, uplow = "low"
    )
    upper_unstrat <- bisect(
      ftn = function(theta) {
        myfun(theta, stratswitch = FALSE) + qnorm(1 - (1 - level) / 2)
      }, contrast = contrast, distrib = distrib,
      precis = precis + 1, uplow = "up"
    )
    if (contrast %in% c("RR", "OR")) {
      lower_unstrat[x1 == 0] <- 0
      upper_unstrat[x2 == 0] <- Inf
    }
  }

  # fix some extreme cases with zero counts
  # if (contrast == "RR" && skew == FALSE) p2d_w[sum(x2) == 0] <- 0

  if (contrast %in% c("RR", "OR")) {
    if (stratified == FALSE) {
      lower[x1 == 0] <- 0
      upper[x2 == 0] <- Inf
    } else {
      lower[sum(x1) == 0] <- 0
      upper[sum(x2) == 0] <- Inf
      point[sum(x2) == 0 & sum(x1) == 0] <- NA
      upper[upper > 1e+15] <- Inf
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
  } else {
    inputs <- NULL
  }

  estimates <- cbind(
    round(cbind(Lower = lower, MLE = point, Upper = upper), precis),
    level = level, inputs,
    p1hat = p1hat_w, p2hat = p2hat_w, p1mle = p1d_w, p2mle = p2d_w
  )
  if (stratified == FALSE) {
    estimates <- cbind(estimates, V = at_MLE$V)
  }

  # optionally add p-value for a test of null hypothesis: theta <= theta0
  # default value of theta0 depends on contrast
  if (contrast == "RD") {
    theta00 <- 0
  } else if (contrast == "p") {
    theta00 <- 0.5
  } else {
    theta00 <- 1
  }
  if (is.null(theta0)) {
    theta0 <- theta00
  }
  scorezero <- scoretheta(
    theta = theta00, x1 = x1, x2 = x2, n1 = n1, n2 = n2,
    stratified = stratified,
    wt = wt, weighting = weighting, MNtol = MNtol,
    random = random, bcf = bcf, contrast = contrast,
    distrib = distrib, skew = skew, ORbias = ORbias,
    RRtang = RRtang, cc = cc, simpleskew = simpleskew, level = level
  )
  scoreth0 <- scoretheta(
    theta = theta0, x1 = x1, x2 = x2, n1 = n1, n2 = n2,
    stratified = stratified, wt = wt,
    weighting = weighting, MNtol = MNtol, random = random,
    bcf = bcf, contrast = contrast, distrib = distrib,
    skew = skew, ORbias = ORbias, RRtang = RRtang,
    cc = cc, simpleskew = simpleskew, level = level
  )
  pval_left <- scoreth0$pval
  pval_left[scoreth0$dsct < 0] <- NA
  pval_right <- 1 - pval_left
  chisq_zero <- scorezero$score^2
  chisq_zero[scorezero$dsct < 0] <- NA
  pval2sided <- pchisq(chisq_zero, 1, lower.tail = FALSE)
  if (random == TRUE) {
    pval2sided <- pf(chisq_zero, 1, nstrat - 1,
      lower.tail = FALSE
    )
  }
  if (simpleskew == TRUE && skew == TRUE) {
    pval <- NULL # simple version of skewness-corrected score is
                 # not valid for producing a p-value
    if (warn == TRUE) {
      warning(paste0("p-values not calculable with simpleskew == TRUE"),
        call. = FALSE
      )
    }
  } else {
    scorenull <- scoreth0$score
    scorenull[scoreth0$dsct < 0] <- NA
    pval <- cbind(
      chisq = chisq_zero, pval2sided, theta0 = theta0,
      scorenull, pval_left, pval_right
    )
  }

  # Add qualitative interaction test as per equation S4 of Laud 2017
  if (stratified == TRUE && nstrat > 1) {
    # V is evaluated at the fixed effects MLE
    Qc_j <- (scoreth0$Stheta)^2 / at_FE$V
    Qc_m <- sum(Qc_j[scoreth0$Stheta > 0])
    Qc_p <- sum(Qc_j[scoreth0$Stheta < 0])
    Qc <- min(Qc_m, Qc_p)
    Qcprob <- 0
    for (h in 1:(nstrat - 1)) {
      Qcprob <- Qcprob + (1 - pchisq(Qc, h)) *
        dbinom(h, size = nstrat - 1, prob = 0.5)
    }
  }

  # Warn for negative discriminant in stratified skewness corrected score
  if (stratified == TRUE) {
    fullseq <- seq(-1, 1, length.out = 400)
    # generalise for equal spacing on transformed scale for RR/OR
    if (contrast %in% c("RR", "OR")) {
      fullseq <- (round(tan(pi * (fullseq + 1) / 4), 10))
    }
    if (contrast == "p") {
      fullseq <- (fullseq + 1) / 2
    }
    # identify affected range from discriminant function
    dt <- array(sapply(fullseq, function(x) mydsct(x)),
      dim = c(1, length(fullseq))
    )
    anydtflag <- FALSE
    if (min(dt, na.rm = TRUE) < 0) {
      if (skew == TRUE && simpleskew == FALSE && random == FALSE) {
        badrange <- range(fullseq[dt < 0], na.rm = TRUE)
        anydtflag <- TRUE
        if (warn == TRUE) {
          warning(paste0(
            "Negative discriminant (min: ", round(min(dt, na.rm = T), 4),
            ") in quadratic skewness corrected score between (",
            paste(round(badrange, 5), collapse = ","),
            "). Simplified skewness correction used in this range."
          ),
          call. = FALSE
          )
        }
        if (scoreth0$dsct < 0 && warn == TRUE) {
          warning(paste0("1-sided p-value not calculable for theta0 = ",
                         theta0),
                  call. = FALSE
          )
        }
        if (scorezero$dsct < 0 && warn == TRUE) {
          warning(paste0("2-sided p-value not calculable"),
                  call. = FALSE
          )
        }
      }
    }
  }

  # Set plot x-axis limits if not specified
  if (is.null(xlim)) {
    if (contrast == "RD") {
      if (distrib == "bin") {
        xlim <- cbind(
          pmax(-1, (lower - (upper - lower) / 2)),
          pmin(1, (upper + (upper - lower) / 2))
        )
      } else if (distrib == "poi") {
        xlim <- cbind(lower - (upper - lower) / 2, upper + (upper - lower) / 2)
      }
    } else if (contrast == "p" && distrib == "bin") {
      xlim <- cbind(
        pmax(0, (lower - (upper - lower) / 2)),
        pmin(1, (upper + (upper - lower) / 2))
      )
    } else {
      xlim <- cbind(pmax(0, (0.5 * lower)), pmin(plotmax, (1.5 * upper)))
    }
  }
  if (!is.array(xlim)) {
    if (stratified == TRUE) {
      xlim <- array(xlim, dim = c(1, length(xlim)))
    } else {
      xlim <- array(rep(xlim, each = nstrat), dim = c(nstrat, length(xlim)))
    }
  }
  # sequence of x-axis values for plotting
  rangen <- 400
  if (stratified == TRUE) {
    dim1 <- 1
    myseq <- array(seq(xlim[, 1], xlim[, 2], length.out = rangen),
                   dim = c(rangen, dim1))
  } else {
    dim1 <- nstrat
    myseq <- array(sapply(1:nstrat, function(i)
      seq(xlim[i, 1], xlim[i, 2], length.out = rangen)),
      dim = c(rangen, dim1))
  }
  # Create flag to identify where negative discriminants occur within the
  # plotting range, i.e. in the vicinity of the confidence interval
  dtseg <- array(sapply(1:rangen, function(i) mydsct(myseq[i, ])),
                 dim = c(dim1, rangen)
  )
  dtflag <- FALSE
  if (min(dtseg, na.rm = TRUE) < 0 && skew == TRUE && !simpleskew) {
    dtflag <- TRUE
  }
  # Optional plot of the score function.
  # Ideally this would be in a separate function, but it is unlikely to be used
  # much in practice - only included for code development and validation
  # purposes.
  if (plot == TRUE || hetplot == TRUE) {
    if (stratified == TRUE && nstrat > 1) {
      if (hetplot == TRUE) { # Note some problems for OR may need fixing
        if (sum(sqrt(V_FE)) > 0) {
          qqnorm(Stheta_FE / sqrt(V_FE))
          abline(coef = c(0, 1))
          plot(
            x = 1 / sqrt(V_FE), y = Stheta_FE / sqrt(V_FE),
            xlab = expression("1/" * sqrt("V"["j"])),
            ylab = expression("S"["j"] * "(" * theta * ")/" * sqrt("V")),
            xlim = c(0, max(1 / sqrt(V_FE))),
            ylim = range(c(-2.5, 2.5, Stheta_FE / sqrt(V_FE))),
            main = expression("Galbraith plot for S"["j"] * "(" * theta * ")")
          )
          abline(coef = c(0, 0))
          abline(coef = c(1.96, 0), lty = 2)
          abline(coef = c(-1.96, 0), lty = 2)
          xrange <- seq(0.1, max(1 / sqrt(V_FE)), length.out = 30)
          lines(xrange, (1.96 * sqrt(1 - xrange^2 / sum(1 / V_FE))), lty = 3)
          lines(xrange, (-1.96 * sqrt(1 - xrange^2 / sum(1 / V_FE))), lty = 3)
        }
      }
    }
    # score for plotting
    sc <- array(sapply(1:rangen, function(i) myfun(myseq[i, ])),
                  dim = c(dim1, rangen)
    )
    # simpleskew version for displaying in event of negative discriminant
    ssc <- array(sapply(1:rangen,
                        function(i) myfun(myseq[i, ], ssswitch = TRUE)),
                dim = c(dim1, rangen)
    )

    if (stratified == FALSE && nstrat > 0) {
      qnval <- qtnorm
      if (is.null(ylim)) ylim <- c(-2.5, 2.5) * qnval
      for (i in 1:nstrat) {
        plot(myseq[, i], sc[i, ],
          type = "l", xlim = xlim[i, ], ylim = ylim,
          xlab = contrast, yaxs = "i", ylab = "Score", col = "blue",
          main = paste0(
            "Score function for ",
            ifelse(distrib == "bin", "binomial", "Poisson"),
            " ", contrast, "\n", x1[i], "/", n1[i],
            ifelse(contrast == "p", "",
              paste0(" vs ", x2[i], "/", n2[i])
            )
          )
          # log = ifelse(contrast == "RD", "", "x")
        )
        if (dtflag == TRUE) {
          lines(myseq[,i], ssc[i, ], col = "grey", lty = 2)
          # lines(myseq, uc[i, ], col = "green", lty = 3)
        }
        text(
          x = c(lower[i], point[i], upper[i]), y = c(-1.5, -1.75, -2) * qnval,
          labels = formatC(c(lower[i], point[i], upper[i]),
            format = "fg",
            4, flag = "#"
          ),
          pos = 4, offset = -0.2, cex = 0.8, xpd = TRUE
        )
        abline(h = c(-1, 1) * qnval)
        abline(h = 0, lty = 2)
        lines(rep(lower[i], 2), c(ylim[1], -1.5 * qnval - 0.3), lty = 3)
        lines(rep(lower[i], 2), c(-1.5 * qnval + 0.4, qnval), lty = 3)
        lines(rep(upper[i], 2), c(ylim[1], -2 * qnval - 0.3), lty = 3)
        lines(rep(upper[i], 2), c(-2 * qnval + 0.4, -qnval), lty = 3)
        lines(rep(point[i], 2), c(ylim[1], -1.75 * qnval - 0.3), lty = 3)
        lines(rep(point[i], 2), c(-1.75 * qnval + 0.4, 0), lty = 3)
      }
    } else if (stratified == TRUE) {
      qtval <- qtnorm
      if (is.null(ylim)) ylim <- c(-2.5, 2.5) * qtval
      # Remove negative discriminant points from plot
      # though they are used in bisection routine
      if (skew == TRUE) sc[1, dtseg < 0] <- NA
      plot(myseq, sc[1, ],
        type = "l", ylim = ylim, xlab = contrast,
        ylab = "Score", yaxs = "i", col = "blue",
        main = paste(
          "Score function for",
          ifelse(distrib == "bin", "binomial", "Poisson"),
          contrast, "\n", paste(x1, collapse = ","), "/",
          paste(n1, collapse = ","),
          ifelse(contrast == "p", "",
            paste0(
              " vs\n ", paste(x2, collapse = ","), " / ",
              paste(n2, collapse = ",")
            )
          )
        )
        # log = ifelse(contrast == "RD", "", "x")
      )
      if (min(dtseg, na.rm = TRUE) < 0 && skew == TRUE && !simpleskew) {
        lines(myseq, ssc[1, ], col = "grey", lty = 2)
        # lines(myseq, uc[1, ], col = "green", lty = 2)
      }
      abline(h = c(-1, 1) * qtval)
      abline(h = 0, lty = 2)
      lines(rep(lower, 2), c(ylim[1], -1.5 * qtval - 0.3), lty = 3)
      lines(rep(lower, 2), c(-1.5 * qtval + 0.4, qtval), lty = 3)
      lines(rep(upper, 2), c(ylim[1], -2 * qtval - 0.3), lty = 3)
      lines(rep(upper, 2), c(-2 * qtval + 0.4, -qtval), lty = 3)
      lines(rep(point, 2), c(ylim[1], -1.75 * qtval - 0.3), lty = 3)
      lines(rep(point, 2), c(-1.75 * qtval + 0.4, 0), lty = 3)
      text(
        x = c(lower, point, upper), y = c(-1.5, -1.75, -2) * qtval,
        labels = formatC(c(lower, point, upper),
          format = "fg", 4,
          flag = "#"
        ),
        pos = 4, offset = -0.2, cex = 0.8, xpd = TRUE
      )
    }
  }

  outlist <- list(estimates = estimates, pval = pval)
  if (stratified == TRUE && nstrat > 1) {
    Qtest <- c(
      Q = Q_FE, pval_het = pval_het, I2 = I2, tau2 = tau2_FE, Qc = Qc,
      pval_qualhet = Qcprob
    ) # Qc_m, Qc_p,
    # NB Qc_m + Qc_p = Q only when theta0=MLE
    wtpct <- 100 * wt_MLE / sum(wt_MLE)
    wt1pct <- 100 * wt_FE / sum(wt_FE)
    if (random == TRUE && prediction == TRUE) {
      outlist <- append(outlist, list(prediction = pred))
    }
    outlist <- append(
      outlist,
      list(
        Qtest = Qtest, weighting = weighting, # dtflag = dtflag, anydtflag = anydtflag,
        stratdata = cbind(
          x1j = x1, n1j = n1, x2j = x2, n2j = n2,
          p1hatj = p1hat, p2hatj = p2hat, wt_fixed = wt_FE,
          wtpct_fixed = wt1pct, wtpct_rand = wtpct,
          theta_j = point_FE_unstrat, lower_j = lower_unstrat,
          upper_j = upper_unstrat, V_j = V_FE, Stheta_j = Stheta_FE,
          Q_j = Q_each
        )
      )
    )
  }
  outlist <- append(outlist, list(call = c(
    distrib = distrib, contrast = contrast, level = level, skew = skew,
    simpleskew = simpleskew, ORbias = ORbias, RRtang = RRtang,
    bcf = bcf, cc = cc, random = random
  )))
  return(outlist)
}

#' Skewness-corrected asymptotic score ("SCAS") confidence intervals for
#' comparisons of independent binomial or Poisson rates.
#'
#' Wrapper function for the SCAS method. Score-based confidence intervals for
#' the rate (or risk) difference ("RD") or ratio ("RR") for independent binomial
#' or Poisson rates, or for odds ratio ("OR", binomial only), or the single rate
#' ("p"). (This is the "GNbc" method from Laud & Dane, developed from Gart &
#' Nam, and generalised as "SCAS" in Laud 2017) including optional
#' continuity correction.  This function is vectorised in x1, x2, n1, and n2.
#' Vector inputs may also be combined into a single stratified analysis (e.g.
#' meta-analysis). This method assumes the contrast is constant across strata
#' (fixed effects).  For a 'random-effects' method use tdasci (or scoreci with
#' random = TRUE).
#'
#' @param x1,x2 Numeric vectors of numbers of events in group 1 & group 2
#'   respectively.
#' @param n1,n2 Numeric vectors of sample sizes (for binomial rates) or exposure
#'   times (for Poisson rates) in each group.
#' @inheritParams scoreci
#' @param theta0 Number to be used in a one-sided significance test (e.g.
#'   non-inferiority margin). 1-sided p-value will be <0.025 iff 2-sided 95\% CI
#'   excludes theta0. By default, a two-sided test against theta0 = 0 (for RD)
#'   or 1 (for RR/OR) is also output.
#' @return A list containing the following components: \describe{
#'   \item{estimates}{a matrix containing estimates of the rates in each group
#'   and of the requested contrast, with its confidence interval} \item{pval}{a
#'   matrix containing details of the corresponding 2-sided significance test
#'   against the null hypothesis that p_1 = p_2, and one-sided significance
#'   tests against the null hypothesis that theta >= or <= theta0}
#'   \item{call}{details of the function call} } If stratified = TRUE, the
#'   following outputs are added: \describe{ \item{Qtest}{a vector of values
#'   describing and testing heterogeneity} \item{weighting}{a string indicating
#'   the selected weighting method} \item{stratdata}{a matrix containing stratum
#'   estimates and weights}}
#' @export
scasci <- function(x1,
                   n1,
                   x2 = NULL,
                   n2 = NULL,
                   distrib = "bin",
                   contrast = "RD",
                   level = 0.95,
                   cc = FALSE,
                   theta0 = NULL,
                   precis = 6,
                   plot = FALSE,
                   hetplot = FALSE,
                   xlim = NULL,
                   ylim = NULL,
                   plotmax = 100,
                   stratified = FALSE,
                   weighting = NULL,
                   MNtol = 1E-8,
                   wt = NULL,
                   warn = TRUE,
                   ...) {
  scoreci(
    x1 = x1,
    n1 = n1,
    x2 = x2,
    n2 = n2,
    distrib = distrib,
    contrast = contrast,
    level = level,
    skew = TRUE,
    simpleskew = FALSE,
    ORbias = TRUE,
    bcf = TRUE,
    cc = cc,
    theta0 = theta0,
    precis = precis,
    plot = plot,
    hetplot = hetplot,
    plotmax = plotmax,
    xlim = xlim,
    ylim = ylim,
    stratified = stratified,
    weighting = weighting,
    MNtol = MNtol,
    wt = wt,
    random = FALSE,
    prediction = FALSE,
    ...
  )
}

#' t-distribution asymptotic score ("TDAS") confidence intervals for
#' comparisons of independent binomial or Poisson rates.
#'
#' Wrapper function for the TDAS method. Score-based stratified confidence
#' intervals for the rate (or risk) difference ("RD") or ratio ("RR") for
#' independent binomial or Poisson rates, or for odds ratio ("OR", binomial
#' only), or for prevalence or incidence rate ("p"). This function combines
#' vector inputs into a single stratified random effects analysis
#' (e.g. meta-analysis), incorporating any stratum variability into the
#' confidence interval.
#'
#' @param x1,x2 Numeric vectors of numbers of events in group 1 & group 2
#'   respectively.
#' @param n1,n2 Numeric vectors of sample sizes (for binomial rates) or exposure
#'   times (for Poisson rates) in each group.
#' @param skew Logical (default TRUE) indicating whether to apply skewness
#'   correction (for the SCAS method recommended in Laud 2017) or not (for
#'   the Miettinen-Nurminen method) to the per-stratum estimates provided
#'   in the output. Has no effect on the TDAS interval itself.
#' @inheritParams scoreci
#' @param theta0 Number to be used in a one-sided significance test (e.g.
#'   non-inferiority margin). 1-sided p-value will be <0.025 iff 2-sided 95\% CI
#'   excludes theta0. By default, a two-sided test against theta0 = 0 (for RD)
#'   or 1 (for RR/OR) is also output.
#' @return A list containing the following components: \describe{
#'   \item{estimates}{a matrix containing estimates of the rates in each group
#'   and of the requested contrast, with its confidence interval} \item{pval}{a
#'   matrix containing details of the corresponding 2-sided significance test
#'   against the null hypothesis that p_1 = p_2, and one-sided significance
#'   tests against the null hypothesis that theta >= or <= theta0}
#'   \item{Qtest}{a vector of values
#'   describing and testing heterogeneity} \item{weighting}{a string indicating
#'   the selected weighting method} \item{stratdata}{a matrix containing stratum
#'   estimates and weights}
#'   \item{call}{details of the function call}}
#' @export
tdasci <- function(x1,
                   n1,
                   x2 = NULL,
                   n2 = NULL,
                   distrib = "bin",
                   contrast = "RD",
                   level = 0.95,
                   cc = FALSE,
                   theta0 = NULL,
                   precis = 6,
                   plot = FALSE,
                   hetplot = FALSE,
                   plotmax = 100,
                   xlim = NULL,
                   ylim = NULL,
                   weighting = NULL,
                   MNtol = 1E-8,
                   wt = NULL,
                   skew = TRUE, # gives SCAS intervals in stratdata by default
                   prediction = FALSE,
                   warn = TRUE,
                   ...) {
  scoreci(
    x1 = x1,
    n1 = n1,
    x2 = x2,
    n2 = n2,
    distrib = distrib,
    contrast = contrast,
    level = level,
    ORbias = TRUE,
    bcf = TRUE,
    cc = cc,
    theta0 = theta0,
    precis = precis,
    plot = plot,
    hetplot = hetplot,
    plotmax = plotmax,
    xlim = xlim,
    ylim = ylim,
    stratified = TRUE,
    weighting = weighting,
    MNtol = MNtol,
    wt = wt,
    random = TRUE,
    prediction = prediction,
    skew = skew,
    simpleskew = FALSE,
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
bisect <- function(ftn,
                   contrast,
                   distrib,
                   precis,
                   max.iter = 100,
                   uplow = "low") {
  tiny <- (10^-(precis)) / 2
  nstrat <- length(eval(ftn(1)))
  hi <- rep(1, nstrat)
  lo <- rep(-1, nstrat)
  dp <- 2
  niter <- 1
  while (niter <= max.iter && any(dp > tiny | is.na(hi))) {
    dp <- 0.5 * dp
    mid <- pmax(-1, pmin(1, round((hi + lo) / 2, 10)))
    # rounding avoids machine precision problem with, e.g. 7/10-6/10
    if (contrast == "RD" && distrib == "bin") {
      scor <- ftn(mid)
    } else if (contrast == "RD" && distrib == "poi") {
      scor <- ftn(round(tan(pi * mid / 2), 10))
      # avoid machine precision producing values outside [-1, 1]
    } else if (contrast %in% c("RR", "OR") ||
      (contrast == "p" && distrib == "poi")) {
      scor <- ftn(round(tan(pi * (mid + 1) / 4), 10))
      # avoid machine precision producing values outside [-1, 1]
    } else if (contrast == "p" && distrib == "bin") {
      scor <- ftn((mid + 1) / 2)
    }
    check <- (scor <= 0) | is.na(scor)
    # ??scor=NA only happens when |p1-p2|=1 and |theta|=1
    # (in which case hi==lo anyway), or if p1=p2=0
    hi[check] <- mid[check]
    lo[!check] <- mid[!check]
    niter <- niter + 1
  }
  if (uplow == "low") {
    best <- lo
  } else {
    best <- hi
  }
  if (contrast == "RD" && distrib == "bin") {
    return(best)
  } else if ((contrast %in% c("RD") && distrib == "poi")) {
    return(tan(best * pi / 2))
  } else if (contrast %in% c("RR", "OR") ||
    (contrast == "p" && distrib == "poi")) {
    return(tan((best + 1) * pi / 4))
  } else if (contrast == "p" && distrib == "bin") {
    return((best + 1) / 2)
  }
}

# Internal function to evaluate the score at a given value of theta, given the
# observed data using the MLE solution (and notation) given in F&M, extended
# in Laud 2017. This function is vectorised in x1,x2
scoretheta <- function(theta,
                       x1,
                       n1,
                       x2 = NULL,
                       n2 = NULL,
                       distrib = "bin",
                       contrast = "RD",
                       bcf = TRUE,
                       skew = TRUE,
                       simpleskew = FALSE,
                       level = 0.95,
                       ORbias = TRUE,
                       RRtang = TRUE,
                       cc = FALSE,
                       stratified = FALSE,
                       wt = NULL,
                       weighting = "IVS",
                       MNtol = 1E-8,
                       random = FALSE,
                       prediction = FALSE,
                       ...) {
  nstrat <- length(n1)
  lambda <- switch(as.character(bcf),
    "TRUE" = (n1 + n2) / (n1 + n2 - 1),
    "FALSE" = 1
  )
  p1hat <- x1 / n1
  p2hat <- x2 / n2
  x <- x1 + x2
  N <- n1 + n2

  bias <- 0 # Initialise bias modified later for OR only

  # RMLE of p1|theta depends on whether theta=RD or RR or OR, and on whether a
  # binomial or Poisson distribution is assumed
  if (contrast == "RD") {
    Stheta <- (p1hat - p2hat) - theta
    if (distrib == "bin") {
      # NOTE: Farrington & Manning notation (p1448) uses
      #       theta for N2/N1, and s_0 for RD
      a <- N
      b <- (n1 + 2 * n2) * theta - N - x
      c_ <- (n2 * theta - N - 2 * x2) * theta + x
      d <- x2 * theta * (1 - theta)
      v <- (b / a / 3)^3 - b * c_ / (6 * a * a) + d / a / 2
      s <- sqrt(pmax(0, (b / a / 3)^2 - c_ / a / 3))
      # NaNs produced when? - machine precision again?
      u <- ifelse(v > 0, 1, -1) * s
      w <- (pi + acos(pmax(-1, pmin(1, ifelse(
        u == 0 & v == 0, 0, v / u^3
      ))))) / 3
      # avoids machine precision errors passing impossible values
      # (outside (-1,1)) to acos
      p2d <- pmin(1, pmax(0, round(2 * u * cos(w) - b / a / 3, 10)))
      # again, machine precision errors can give values a tiny fraction either
      # side of 0 or 1
      p1d <- pmin(1, pmax(0, p2d + theta))
      V <- pmax(0, (p1d * (1 - p1d) / n1 + p2d * (1 - p2d) / n2) * lambda)
      # V set to zero if machine precision error in cos function produces
      # a negative V
      mu3 <- p1d * (1 - p1d) * (1 - 2 * p1d) / (n1^2) -
        p2d * (1 - p2d) * (1 - 2 * p2d) / (n2^2)
    } else if (distrib == "poi") {
      A <- N
      B <- N * theta - x
      C_ <- -x2 * theta
      num <- (-B + sqrt(B^2 - 4 * A * C_))
      p2d <- ifelse(num == 0, 0, num / (2 * A))
      p1d <- p2d + theta
      V <- (p1d / n1 + p2d / n2)
      mu3 <- (p1d / (n1^2) - p2d / (n2^2))
    }
  } else if (contrast == "RR") {
    Stheta <- p1hat - p2hat * theta
    if (distrib == "bin") {
      A <- N * theta
      B <- (-(n1 * theta + x1 + n2 + x2 * theta))
      C_ <- x
      num <- (-B - sqrt(pmax(0, B^2 - 4 * A * C_)))
      p2d <- ifelse(A == 0, -C_ / B, ifelse(
        (num == 0 | C_ == 0), 0, num / (2 * A)
      ))
      p1d <- p2d * theta
      V <- pmax(0, (p1d * (1 - p1d) / n1 +
        (theta^2) * p2d * (1 - p2d) / n2) * lambda)
      mu3 <- (p1d * (1 - p1d) * (1 - 2 * p1d) / (n1^2) -
        (theta^3) * p2d * (1 - p2d) * (1 - 2 * p2d) / (n2^2))
      if (RRtang == TRUE) {
        Stheta <- (p1hat - p2hat * theta) / p2d
        V <- pmax(0, (p1d * (1 - p1d) / n1 +
          (theta^2) * p2d * (1 - p2d) / n2) * lambda / p2d^2)
        mu3 <- (p1d * (1 - p1d) * (1 - 2 * p1d) / (n1^2) -
          (theta^3) * p2d * (1 - p2d) * (1 - 2 * p2d) / (n2^2)) / p2d^3
        mu3[(x1 == 0 & x2 == 0)] <- 0
        Stheta[(x1 == 0 & x2 == 0)] <- 0
      }
      V[is.na(V)] <- Inf
    } else if (distrib == "poi") {
      # incidence density version from M&N p223
      p2d <- (x1 + x2) / (n1 * theta + n2)
      p1d <- p2d * theta
      V <- pmax(0, (p1d / n1 + (theta^2) * p2d / n2))
      mu3 <- (p1d / (n1^2) - (theta^3) * p2d / (n2^2))
      if (RRtang == TRUE) {
        # Apply Tang score for Poisson parameter
        # EXPERIMENTAL - needs to be evaluated
        Stheta <- (p1hat - p2hat * theta) / p2d
        V <- pmax(0, (p1d / n1 + (theta^2) * p2d / n2) / p2d^2)
        mu3 <- (p1d / (n1^2) - (theta^3) * p2d / (n2^2)) / p2d^3
        mu3[(x1 == 0 & x2 == 0)] <- 0
        Stheta[(x1 == 0 & x2 == 0)] <- 0
      }
    }
  } else if (contrast == "OR") {
    if (distrib == "bin") {
      A <- n2 * (theta - 1)
      B <- n1 * theta + n2 - x * (theta - 1)
      C_ <- -x
      num <- (-B + Re(sqrt(as.complex(B^2 - 4 * A * C_))))
      # If A=0 then we solve a linear equation instead
      p2d <- ifelse(A == 0, -C_ / B, ifelse(num == 0, 0, num / (2 * A)))
      p1d <- p2d * theta / (1 + p2d * (theta - 1))
      p1d[theta == 0] <- 0
      Stheta <- (p1hat - p1d) / (p1d * (1 - p1d)) -
        (p2hat - p2d) / (p2d * (1 - p2d))
      V <- pmax(0, (1 / (n1 * p1d * (1 - p1d)) +
        1 / (n2 * p2d * (1 - p2d))) * lambda)
      # set to zero if machine precision error in cos function
      # produces a negative V
      mu3 <- (1 - 2 * p1d) / ((n1 * p1d * (1 - p1d))^2) -
        (1 - 2 * p2d) / ((n2 * p2d * (1 - p2d))^2)
      Stheta[(x1 == 0 & x2 == 0) | (x1 == n1 & x2 == n2)] <- 0
      mu3[(x1 == 0 & x2 == 0) | (x1 == n1 & x2 == n2)] <- 0

      if (ORbias == TRUE) {
        bias <- (p1d - p2d) / (n1 * p1d * (1 - p1d) + n2 * p2d * (1 - p2d))
        bias[p1d == p2d] <- 0
        Stheta <- Stheta - bias
      }
      # note V = V_G^-1 = (npq1 + npq2)/(npq1npq2) and this version puts bias
      # within (Stheta-bias)/sqrt(V) instead of Stheta/sqrt(V) - bias
    } else if (distrib == "poi") {
      print("Odds ratio not applicable to Poisson rates")
    }
  } else if (contrast == "p") {
    Stheta <- p1hat - theta
    Stheta[n1 == 0] <- 0
    if (distrib == "bin") {
      V <- (pmax(0, (theta * (1 - theta) / n1)))
      mu3 <- (theta * (1 - theta) * (1 - 2 * theta) / (n1^2))
    } else if (distrib == "poi") {
      V <- theta / n1
      mu3 <- theta / (n1^2)
    }
    mu3[n1 == 0] <- 0 # quick fix
    p1d <- theta
    p2d <- NA
  }

  # continuity corrections
  corr <- 0
  if (cc > 0) {
    if (contrast == "OR") {
      corr <- cc * (1 / (n1 * p1d * (1 - p1d)) + 1 / (n2 * p2d * (1 - p2d)))
      # cc=0.5 gives Cornfield correction. Try cc=0.25 instead
    } else if (contrast == "RR") {
      corr <- cc * (1 / (n1) + theta / (n2)) # try 0.125 or 0.25
      if (RRtang == TRUE) corr <- corr / p2d
    } else if (contrast == "RD") {
      corr <- cc * (1 / pmin(n1, n2)) # cc=0.5 gives Hauck-Anderson. Try 0.25
      if (stratified == TRUE) {
        corr <- (3 / 16) * (sum(n1 * n2 / (n1 + n2)))^(-1)
      }
      # from Mehrotra & Railkar, also Zhao et al.
    } else if (contrast == "p") {
      corr <- cc / n1
    }
  }

  if (stratified == TRUE) {
    pval <- NA
    if (is.null(wt)) {
      if (weighting == "MH") {
        wt <- n1 * n2 / (n1 + n2)
        # MH weights for RD, applied across other comparative parameters too
        # (without theoretical justification for OR)
        if (contrast == "p") wt <- n1
      } else if (weighting == "IVS") {
        # IVS: inverse variance weights updated wih V_tilde
        if (all(V == 0) || all(V == Inf | is.na(V))) {
          wt <- rep(1, nstrat)
        } else {
          wt <- 1 / V
        }
      } else if (weighting == "INV") {
        # INV: inverse variance weights updated wih V_tilde
        # Removing bcf from the weights ensures alignment with CMH test
        if (all(V == 0) || all(V == Inf | is.na(V))) {
          wt <- rep(1, nstrat)
        } else {
          wt <- lambda / V
        }
      } else if (weighting == "MN") {
        if (contrast == "RR" && distrib == "poi") {
          wt <- (1 / n1 + theta / n2)^(-1)
        } else if (contrast == "RR" && distrib == "bin") {
          # M&Ns iterative weights - quite similar to MH
          wtdiff <- 1
          wtx <- (1 / n1 + theta / n2)^(-1)
          while (wtdiff > MNtol) {
            p2ds <- sum(wtx * p2d) / sum(wtx)
            p1ds <- sum(wtx * p1d) / sum(wtx)
            wt <- ((1 - p1ds) / (n1 * (1 - p2ds)) + theta / n2)^(-1)
            wt[p2ds == 1] <- 0
            wtdiff <- max(abs(wtx - wt))
            wtx <- wt
          }
        } else if (contrast == "RD") {
          # M&Ns iterative weights - quite similar to MH wtx <- n1*n2/(n1+n2)
          wt <- wtx <- (1 / n1 + 1 / n2)^(-1)
          wtdiff <- 1
          while (wtdiff > MNtol) {
            p2ds <- sum(wtx * p2d) / sum(wtx)
            p1ds <- sum(wtx * p1d) / sum(wtx)
            if (distrib == "bin") {
              wt <- (
                (p1ds * (1 - p1ds) / (p2ds * (1 - p2ds))) / n1 + 1 / n2
              )^(-1)
            } else if (distrib == "poi") {
              wt <- ((p1ds / p2ds) / n1 + 1 / n2)^(-1)
            }
            wt[p2ds == 0 || p2ds == 1] <- 0
            if (sum(wt) == 0) wt <- wt + 1 # Fix for when all weights are zero
            wtdiff <- max(abs(wtx - wt))
            wtx <- wt
          }
        } else if (contrast == "OR") {
          # M&Ns weights are very similar in structure to IVS
          wt <- n1 * n2 * ((1 - p1d) * p2d)^2 /
            (n1 * p1d * (1 - p1d) + n2 * p2d * (1 - p2d))
        }
      }
    } else {
      weighting <- "User-defined"
    }

    Sdot <- sum(wt * Stheta) / sum(wt)
    # NB the skewness correction is omitted for the heterogeneity test statistics.
    Q_j <- ((Stheta - Sdot)^2) / V
    # NB it is necessary to include Sdot here for TDAS method to work.
    # - for the heterogeneity test evaluated at MLE,
    #   Sdot will equal 0 if skew = FALSE
    # NB it is necessary to use equation S2 here for TDAS method to work
    Q <- sum(Q_j)
    W <- sum(wt)

    if (weighting %in% c("IVS", "INV")) { # Check if this works for INV as well
      tau2 <- max(0, (Q - (nstrat - 1))) / (W - (sum(wt^2) / W))
      # published formula for IVS weights
    } else {
      tau2 <- max(0, (Q - (nstrat - 1 * (sum(V * wt^2) / W)))) /
        (sum(1 / V) - 1 * (sum(wt^2) / W))
      # only needed if want to output tau2.
    }
    if (
      random == TRUE && weighting %in% c("IVS", "INV") &&
        !(all(V == 0) || all(V == Inf | is.na(V)))
    ) {
      wt <- 1 / (V + tau2)
    }

    Sdot <- sum(wt * Stheta) / sum(wt)
    Vdot <- sum(((wt / sum(wt))^2) * V)
    # Alternative form for IVS weights avoids the need to exclude
    # non-informative strata for OR. For use in possible future update
    if (weighting %in% c("IVS")) Vdot <- sum(wt / (sum(wt))^2)
    if (weighting %in% c("INV")) Vdot <- sum(lambda * wt / (sum(wt))^2)

    if (contrast == "OR" && cc > 0) {
      # corr <- cc * Vdot # This gives essentially the same correction as Gart 1985
      corr <- cc * sum(((wt / sum(wt))^2) * (1 / (n1 * p1d * (1 - p1d)) +
        1 / (n2 * p2d * (1 - p2d))))
    }
    if (contrast == "RR" && cc > 0) { # EXPERIMENTAL cc for stratified RR
      # Tentative
      corr <- cc * sum(((wt / sum(wt))^2) * (1 / (n1) + theta / (n2)))
      # corr <- cc * Vdot #more tentative - I think wrong
    }
    corr <- corr * sign(Sdot)
    score1 <- sum((wt / sum(wt)) * (Stheta - corr)) / pmax(0, sqrt(Vdot))
    score1[sum((wt / sum(wt)) * (Stheta - corr)) == 0] <- 0 # Avoids NA scores
    scterm <- sum(((wt / sum(wt))^3) * mu3) / (6 * Vdot^(3 / 2))
    scterm[sum(((wt / sum(wt))^3) * mu3) == 0] <- 0 # Avoids NA scores
    A <- scterm
    B <- 1
    C_ <- -(score1 + scterm)
    num <- (-B + sqrt(pmax(0, B^2 - 4 * A * C_)))
    dsct <- B^2 - 4 * A * C_
    score <- ifelse(
      (skew == FALSE | scterm == 0),
      score1, num / (2 * A)
    )
    if (skew == TRUE) {
      qtnorm <- qnorm(1 - (1 - level) / 2)
      scoresimp <- score1 - (qtnorm^2 - 1) * scterm
      # Bartlett version, not recommended:
      # scorebart <- score1 - (score1^2 - 1) * scterm
      if (simpleskew == TRUE) {
        score <- scoresimp
      }
      if (!is.na(dsct) & dsct < 0) {
        score <- scoresimp
      }
    }
    score[abs(Sdot) < abs(sum((wt / sum(wt)) * corr))] <- 0
    pval <- pnorm(score)
    VS <- sum(wt / sum(wt) * (Stheta - Sdot)^2) / (nstrat - 1)
    if (random == TRUE) {
      tnum <- sum((wt / sum(wt)) * (Stheta - corr))
      tdenom <- sqrt(VS)
      if (prediction == TRUE) {
        # Prediction interval from Higgins et al. 2009,
        # using Hartung-Knapp variance estimate
        tdenom <- sqrt(VS + tau2)
      }
      score1 <- tnum / tdenom
      ##  Aborted (& unnecessary) attempt to add skewness correction
      ## to TDAS (for contrast="p")
      # score1 <- (tnum + scterm)/tdenom
      ## Or
      # A <- scterm
      # B <- 1
      # C_ <- -(score1 + scterm)
      # num <- (-B + sqrt(max(0, B^2 - 4 * A * C_)))
      # score <- ifelse((skew == FALSE | scterm == 0), score1, num/(2 * A))
      # score[abs(Sdot) < abs(sum( (wt/sum(wt)) * corr))] <- 0
      score <- score1
      pval <- pt(score, nstrat - 1)
    }
    p2ds <- sum(wt * p2d / sum(wt))
    p1ds <- sum(wt * p1d / sum(wt))
  } else if (stratified == FALSE) {
    p1ds <- p1d
    p2ds <- p2d

    corr <- corr * sign(Stheta)
    # Calculation of score & p-value involves solving for z_p:
    # z_p = Stheta / sqrt(V) - (z_p^2) * mu3 / (6 * V^(3 / 2)) +
    #       mu3 / (6 * V^(3 / 2))
    # Note that in the special case of mu3=0, this reduces to
    # the skew = FALSE case.
    # i.e. z_p = Stheta / sqrt(V)
    scterm <- mu3 / (6 * V^(3 / 2))
    scterm[mu3 == 0] <- 0
    score1 <- (Stheta - corr) / sqrt(V)
    score1[Stheta == 0] <- 0
    A <- scterm
    B <- 1
    C_ <- -(score1 + scterm)
    num <- (-B + sqrt(pmax(0, B^2 - 4 * A * C_)))
    dsct <- B^2 - 4 * A * C_
    score <- ifelse((skew == FALSE | scterm == 0),
      score1, num / (2 * A)
    )
    if (skew == TRUE) {
      qtnorm <- qnorm(1 - (1 - level) / 2)
      scoresimp <- score1 - (qtnorm^2 - 1) * scterm
      # Bartlett version, not recommended
      # scorebart <- score1 - (score1^2 - 1) * scterm
      if (simpleskew == TRUE) {
        score <- scoresimp
      }
      # below is unnecessary as this never happens for unstratified datasets
      score[!is.na(dsct) & dsct < 0] <- scoresimp[!is.na(dsct) & dsct < 0]
    }
    score[abs(Stheta) < abs(corr)] <- 0

    pval <- pnorm(score)
  }

  outlist <- list(
    score = score, p1d = p1d, Stheta = Stheta, num = num, V = V,
    p2d = p2d, mu3 = mu3, pval = pval, dsct = dsct
  )
  if (stratified) {
    outlist <- append(outlist, list(
      Sdot = Sdot, Vdot = Vdot, tau2 = tau2,
      VS = VS, Q_j = Q_j, Q = Q, wt = wt, p1ds = p1ds,
      p2ds = p2ds
    ))
  }
  return(outlist)
}
