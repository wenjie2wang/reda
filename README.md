# reda

The R package **reda** mainly provides function to fit gamma frailty model with
either a piecewise constant or a spline as the baseline rate function
for recurrent event data. What's more, some handy functions are designed,
such as computing and plotting sample nonparametric mean cumulative function,
or so-called Nelson-Aalen estimator. Most functions in this package
are S4 methods that produce S4 class objects.


## Development

The latest version of package is under development in branch 'dev'.


## Installation of Stable Version

You can install the stable version in branch 'master' or on
[CRAN](http://cran.rstudio.com/package=reda):

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/reda)](http://cran.r-project.org/package=reda)

```r
install.packages("reda", dependencies = TRUE)
```


## Basic Usage

```r
help(pacakge = "reda")
library(reda)
## help on main function for model fitting
?rateReg
```

See [package help manual](https://cran.rstudio.com/web/packages/reda/reda.pdf)
and [vignette](https://cran.rstudio.com/web/packages/reda/vignettes/reda-intro.html)
for details and demonstration.


## License

The R package reda is free software: You can redistribute it and/or
modify it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or
any later version (at your option).
See the [GNU General Public License](http://www.gnu.org/licenses/) for details.

The R package reda is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
