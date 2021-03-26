# ratesci 0.3-0.9000

## To-do
* [Add checks on valid determinants when solving quadratics 
  - produce warning & highlight on plot of score function.]
* [Check input vector lengths match]

## New features
* MN weighting in `scoreci()` now iterates to convergence (@jonjvallejo, #20).
* `tdasci()` default uses skew=TRUE for stratum CIs.
* Added optional prediction interval for random effects method.
* Added xlim argument to control plot output.
* Added sda argument to control sparse data adjustment (needs work to improve versatility).
* Added INV option for weights that omit the bias correction.
* Added RRtang argument to apply Tang's alternative score for stratified RR with INV/IVS weights.
* Added simplified skewness correction option.
* Changed ORbias default to TRUE.

## Bug fixes
* MN weighting in `scoreci()` corrected for distrib="poi".
* Fixed bug in `scoreci()` for calculation of stratum CIs with random=TRUE.
* Fixed bug in `scoreci()` for distrib="poi" and contrast="p" (#7).
* Fixed finite precision bug in `scaspci()`.
* Fixed bug in `rateci()` for closed-form calculation of continuity-corrected SCAS.

## Other
* Renamed tdas argument to 'random'.
* Removed redundant t2 variable.

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
