# ratesci 0.3-0

## New features
* Added bias correction in `scoreci()` for OR SCAS method (derived from Gart 1985).
* Added score methods (Tango & Tang) as default for paired binomial RD and RR in `pairbinci()`.
* Added transformed mid-p method for paired OR for comparison with transformed SCAS.
* Added `scaspci()` for non-iterative SCAS methods for single binomial or Poisson rate.
* Added `rateci()` for selected methods for single binomial or Poisson rate.

## Bug fixes
* Fixed bug in `pairbinci()` for contrast="OR".
* Fixed bug in `moverci()` for contrast="p" and type="wilson".
* Corrected error in cc for stratified SCAS method for OR.
* Clarified documentation regarding continuity corrections.
* Set Stheta to 0 if |Stheta|<cc in `scoreci()`
* Fixed stratified calulations for contrast = "p" in `scoreci()`.

# ratesci 0.2-0

## New features
* Added `pairbinci()` for all comparisons of paired binomial rates.
* Added option to suppress warnings in scoreci.
* Added Galbraith plot (for assessing stratum heterogeneity) to `scoreci()`.
* Added qualitative interaction test to `scoreci()`.
* Added stratum estimates & CIs to `scoreci()` output when stratified = TRUE.

## Bug fixes
* Fixed bug for contrast = "p" in `moverci()`.
* Fixed bug in `tdasci()` wrapper function.
* Fixed bug for stratified OR.
* Altered adjustment options for boundary cases in `moverci()`.
* Changed point estimate used in `moverci()` to posterior median for type = "jeff",
  to ensure consistent calculations with informative priors.
