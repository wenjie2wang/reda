# reda v0.2.1.9000

## Major changes

* Borrowed the power from R package **splines2** for piecewise constat and
  splines based baseline rate function, and thus improved the performance of
  **reda** in model fitting.


# reda v0.2.1

## Bug fixes

* Updated S4 method `plotMcf,sampleMcf`: Replaced `show_guide` with
  `show.legend` in function `geom_text` to incorporate updates in package
  `ggplot2` v1.0.1.

* Minor updates that clear checking note from CRAN.


# reda v0.2.0

## New features

* Renamed main function name from `heart` to `rateReg` and added new argument.

* Implementation of spline baseline rate function.

* Updated object class of fitted model.

* Added function `AIC` and `BIC`.

* Replaced sample simulated dataset for demonstration.


# reda v0.1.0

## New features

* First version of reda mainly providing function to fit gamma frailty model
  with piecewise constant baseline rate function for recurrent event data.


