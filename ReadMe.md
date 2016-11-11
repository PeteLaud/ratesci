### R/ratesci: Confidence intervals for comparisons of binomial or Poisson rates
ratesci
=====

[Pete Laud] #(http://kbroman.org)

[R/ratesci] is an [R](http://www.r-project.org) package to compute confidence intervals for rate (or risk) difference (RD), ratio (RR) or odds ratio, for binomial or Poisson rates. It incorporates novel "skewness-corrected" asymptotic score methods, which ensure equal-tailed coverage, in other words for a nominal 95\% confidence interval, the one-sided non-coverage probability is (on average) close to 2.5% on each side. For large sample sizes, similar coverage is achieved with the more approximate 'MOVER-B' methods. Being based on Bayesian methods, these also allow the option to incorporate prior beliefs about the rates in each group - otherwise, the 'non-informative' Jeffreys priors are used.

#### Installation

You can install:

#-   the latest released version from CRAN with

#Install R/qtlcharts from CRAN using

#```r
#install.packages("ratesci")
#```

-   the latest development version from the [GitHub repository](https://github.com/PeteLaud/ratesci) with

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


#### Licenses

The R/ratesci package is distributed under
[GPL-3 (GNU General Public License version 3)[http://www.gnu.org/licenses/gpl-3.0.en.html].

