# ratesci 0.2-0.9000

* Clarified documentation regarding continuity corrections
* Fixed bug in pairbinci for contrast="OR"
* Fixed bug in moverci for contrast="p" and type="wilson"
* Added score methods for paired binomial RD and RR (Tango & Tang)
* Added non-iterative SCAS methods for single binomial or Poisson rate
* Added transformed mid-p method for paired OR for comparison with transformed SCAS
* Added bias correction to SCAS method for OR (based on Gart 1985)
* Corrected error in cc for stratified SCAS method for OR

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
