# ratesci 0.2-0

## New features
* Added pairbinci function for all comparisons of paired binomial rates.
* Added option to suppress warnings in scoreci.
* Added Galbraith plot (for assessing stratum heterogeneity) to scoreci.
* Added qualitative interaction test to scoreci.
* Added stratum estimates & CIs to scoreci output when stratified = TRUE.

## Bug fixes
* Fixed bug for contrast = "p" in moverci.
* Fixed bug in tdasci wrapper function.
* Fixed bug for stratified OR.
* Altered adjustment options for boundary cases in moverci.
* Changed point estimate used in moverci to posterior median for type = "jeff",
  to ensure consistent calculations with informative priors.
