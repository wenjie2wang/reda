# reda

[![CRAN_Status_Badge][r-pkg-badge]][cran-url]
[![Build Status][gha-icon]][gha-url]
[![codecov][codecov-main]][codecov]

## Overview

The R package **reda** provides functions for

- simulating survival, recurrent event, and multiple event data from stochastic
  process point of view;
- exploring and modeling recurrent event data through the mean cumulative
  function (MCF) by the Nelson-Aalen estimator of the cumulative hazard rate
  function, and gamma frailty model with spline rate function;
- comparing two-sample recurrent event responses with the pseudo-score tests.


## Installation of CRAN Version

You can install the released version from [CRAN][cran-url].

```R
install.packages("reda")
```

## Development

The latest version of the package is under development at [GitHub][github-url].
If it is able to pass the automated package checks, one may install it by

```R
if (! require(remotes)) install.packages("remotes")
remotes::install_github("wenjie2wang/reda", upgrade = "never")
```

## Getting Started

[Online documentation][homepage] provides function documentations and includes
package vignettes for

- [exploring and modeling recurrent event data][reda-intro].
- [introduction to formula response function Recur()][reda-Recur]
- [simulating survival and recurrent event data][reda-simulate].


## License

[GNU General Public License][gpl] (≥ 3)


[r-pkg-badge]: https://www.r-pkg.org/badges/version/reda
[cranlog-badge]: https://cranlogs.r-pkg.org/badges/reda
[cran-url]: https://CRAN.R-project.org/package=reda
[github-url]: https://github.com/wenjie2wang/reda
[gha-icon]: https://github.com/wenjie2wang/reda/workflows/R-CMD-check/badge.svg
[gha-url]: https://github.com/wenjie2wang/reda/actions
[codecov]: https://codecov.io/gh/wenjie2wang/reda
[codecov-main]: https://codecov.io/gh/wenjie2wang/reda/branch/main/graph/badge.svg
[homepage]: https://wwenjie.org/reda/
[reda-intro]: https://wwenjie.org/reda/articles/reda-intro.html
[reda-Recur]: https://wwenjie.org/reda/articles/reda-Recur.html
[reda-simulate]: https://wwenjie.org/reda/articles/reda-simulate.html
[gpl]: https://www.gnu.org/licenses/
