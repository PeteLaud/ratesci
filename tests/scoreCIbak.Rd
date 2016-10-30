\name{scoreCI}
\alias{scoreCI}
%%\alias{score}
%%\alias{bisect}
\title{Confidence Intervals for Comparison of Rates}
\description{
Score-based confidence intervals for the rate difference (RD) or ratio (RR), for independent binomial or Poisson rates, or for odds ratio (OR) (binomial only)
%%including options for variance bias correction (MN), skewness correction (GN) and continuity correction 
}
\usage{
scoreCI(x1, x2=NULL, n1, n2=NULL, dist = "bin", contrast = "RD", level = 0.95,
skew = TRUE, bcf = TRUE, cc = 0, delta = NULL, precis = 6, plot=FALSE,
plotmax=1000, stratified=FALSE, weights="IV", wt=NULL, tdas=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{vector of numerators (number(s) of events) for group 1}
  \item{x2}{vector of numerators for group 2}
  \item{n1}{vector of denominators for group 1 (numbers of subjects for binomial, or exposure durations for Poisson rates)}
  \item{n2}{vector of denominators for group 2}
  \item{dist}{indication of whether data represent binomial ("bin") or Poisson ("poi") rates}
  \item{contrast}{indication of which contrast is required ("RD", "RR" or "OR")}
  \item{level}{confidence level}
  \item{skew}{option to include skewness corrections from Laud 2016, for improved equal-tailed performance}
  \item{bcf}{option allow omission of variance bias correction factor from Miettinen & Nurminen, for comparison with legacy methods (e.g. Gart & Nam)}
  \item{cc}{option to include "continuity correction" for more conservative coverage}
  \item{delta}{optional non-inferiority margin for one-sided significance test}
%  \item{leftside}{option to decide direction of one-sided test. TRUE (default) tests H0: theta<=delta, FALSE tests H0: theta>=delta}
  \item{precis}{required level of precision for output of confidence intervals (number of dps)}
  \item{plot}{option to produce plot of the score function, for checking the outcome of the bisection routine}
  \item{plotmax}{upper limit of x-axis for plot}
  \item{stratified}{stratified option}
  \item{weights}{weights ["MH" (recommended for tdas=FALSE), "IV" (recommended for tdas=TRUE) or "MN")}
  \item{wt}{user-defined weights}
  \item{tdas}{t-distribution asymptotic score method}
  \item{...}{other args}
%  \item{omitlower}{utility option to omit calculation of lower confidence limit to reduce computation time}
%  \item{omitupper}{utility option to omit calculation of upper confidence limit to reduce computation time}
%%  \item{stratified}{option to perform stratified calculations}
%%  \item{random}{option to use random effects model}
}
\details{
[See article submitted to Pharmaceutical Statistics journal.]
}
\value{
A list containing the following components
\item{estimates}{a matrix containing estimates of the rates in each group and of the requested contrast, with its confidence interval}
\item{pval}{a matrix containing details of the corresponding 2-sided and one-sided significance tests}
\item{nstrat}{numeric indicating the number of strata included in the analysis}
%%\item{Qtest}{a vector containing details of the heterogeneity test, if requested}
%%\item{wt}{vector of stratum weights, if stratified==TRUE}
}
\references{
[Laud PJ. Equal-tailed confidence intervals for comparison of rates: Submitted to Pharmaceutical Statistics for peer review]

Miettinen OS, Nurminen M. Comparative analysis of two rates. Statistics in Medicine 1985; 4:213-226.

Gart JJ, Nam JM. Approximate interval estimation of the ratio of binomial parameters: A review and corrections for skewness. Biometrics 1988; 44(2):323-338.

Gart JJ, Nam JM. Approximate interval estimation of the difference in binomial parameters: correction for skewness and extension to multiple tabes. Biometrics 1990; 46(3):637-643.

}
\author{Peter J Laud}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Intervals from Hartung & Knapp
#SCAS:
#ratesCI(x1=12, x2=1, n1=16, n2=16, dist="bin", contrast="RD", skew=TRUE) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
