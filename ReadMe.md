ratesci
=====

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/ratesci)](http://cran.r-project.org/package=ratesci)
[![Travis-CI Build Status](https://travis-ci.org/PeteLaud/ratesci.svg?branch=master)](https://travis-ci.org/PeteLaud/ratesci)

ratesci is an [R](http://www.r-project.org) package to compute confidence intervals for rate (or risk) difference ('RD'), rate ratio ('RR', also known as relative risk), or odds ratio ('OR'). All three contrasts apply for binomial proportions, and the first two may also be used for the comparison of Poisson 'exposure-adjusted' incidence rates. The scoreci function incorporates 'skewness-corrected' asymptotic score ('SCAS') methods, which ensure equal-tailed coverage (or central location), in other words for a nominal 95% confidence interval, the one-sided non-coverage probability is (on average) close to 2.5% on each side. Stratified calculations are also catered for (e.g. meta-analysis, including random effects), as well as confidence intervals for the single binomial or Poisson rate, and for matched pairs (with the pairbinci function).  Corresponding hypothesis tests against any specified null parameter value are provided in each case.  Omission of the skewness correction is also allowed, resulting in the often-recommended 'Miettinen-Nurminen' asymptotic score methods, which can have inferior one-sided coverage, especially for RR.

For large (single-stratum) sample sizes, the 'MOVER-B' methods (moverci function) improve on traditional approximate methods with respect to one-sided and two-sided coverage, particularly in the case of RR. Being based on Bayesian methods, these also allow the option to incorporate prior beliefs about the rates in each group - by default, the 'non-informative' Jeffreys priors are used. These methods are adapted from the Newcombe 'square-and-add' method, which is also included for reference.

For those wishing to achieve strictly conservative coverage, so-called 'continuity corrections' are provided as approximations to 'exact' methods. The performance of these adjustments has not been extensively evaluated, but they appear to be more successful for SCAS than for MOVER, in terms of achieving conservative coverage. 

An online calculator based on this package is available [here](http://ssu.shef.ac.uk/ratesci/calc.php)


#### Installation

You can install:

- the latest released version from CRAN with

    ``` r
    install.packages("ratesci")
    ```

- the latest development version from the [GitHub repository](https://github.com/PeteLaud/ratesci) with

    ``` r
    if (packageVersion("devtools") < 1.6) {
      install.packages("devtools")
    }
    devtools::install_github("PeteLaud/ratesci")
    ```

#### Example use

```r
library(ratesci)
scoreci(x1 = 5, n1 = 56, x2 = 0, n2 = 29)
```


