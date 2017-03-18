# ratesci 0.1-0.9000

## New features
* Added pairbinci function for all comparisons of paired binomial rates
* Added option to suppress warnings in scoreci
* Added Galbraith plot to scoreci

## Bug fixes
* Fixed bug for contrast = "p" in moverci
* Fixed bug in tdasci wrapper function
* Fixed bug for stratified OR
* Altered adjustment options for boundary cases in moverci 
* Changed point estimate used in moverci to posterior median for type = "jeff",
  to ensure consistent calculations with informative priors
