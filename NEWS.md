# ratesci 0.4-0 (2021-12-04)

## New features
### In `scoreci()`:
* MN weighting now iterates to convergence (@jonjvallejo, #20).
* Added optional prediction interval for random effects method (also in `tdasci()`).
* Added xlim and ylim arguments to control plot output.
* Added sda & fda arguments for optional sparse/full data adjustment 
    when x1 + x2 = 0 or x1 + x2 = n1 + n2 in a stratum.
* Added INV option for weights that omit the variance bias correction.
* Added RRtang argument to apply Tang's alternative score for RR (recommended 
for stratified analysis with INV/IVS weights. Experimental for Poisson RR).
    `Stheta = (p1hat - p2hat * theta) / p2d`  (see Tang 2020)
* Added simplified skewness correction option (causes p-values to be omitted, see Tang 2021 & Laud 2021).
* Introduced warning and plot features for very rare occasions when quadratic 
  skewness correction cannot be calculated due to a negative discriminant.
* p-value suppressed where affected by negative discriminants.
* Changed ORbias default to TRUE (see Laud 2018).
* Changed weighting default to MH for RD & RR, INV for OR (for consistency with CMH test).
* Added hetplot argument to separate heterogeneity plots from score function plot.
* Uninformative strata are now retained in the analysis except if: 
  * contrast = OR with MH weighting;
  * contrast = RR with IVS/INV weighting if RRtang = FALSE;
  * random = TRUE (needs further evaluation);
  * excluded using new option dropzeros = TRUE.
### In `tdasci()`:
* Default uses skew = TRUE for stratum CIs.

## Bug fixes
* MN weighting in `scoreci()` corrected for distrib="poi".
* Fixed bug in `scoreci()` for calculation of stratum CIs with random=TRUE.
* Fixed bug in `scoreci()` for distrib = "poi" and contrast = "p" (#7).
* Fixed finite precision bug in `scaspci()`.
* Fixed bug in `rateci()` for closed-form calculation of continuity-corrected SCAS.
* Fixed bug in `scoreci()` for stratified zero scores calculated as NA, resulting in UL = 0. (Thanks to Lidia Mukina for reporting the bug.)
* Fixed variable plot ranges for vectorised inputs.

## Other
* Renamed tdas argument to 'random'.
* Removed redundant t2 variable.

# ratesci 0.3-0 (2018-02-15)

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

# ratesci 0.2-0 (2017-04-21)

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
