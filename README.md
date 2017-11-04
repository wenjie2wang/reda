# reda

The R package **reda** provides functions that fit gamma frailty model with
spline baseline rate function for recurrent event data, compute and plot
parametric mean cumulative function (MCF) from a fitted model as well as
nonparametric sample MCF (Nelson-Aalen estimator) are provided.  Most functions
are S4 methods that produce S4 class objects.


## Installation of CRAN Version

[![CRAN_Status_Badge][r-pkg-badge]][cran-url]
[![Build Status][travis-master]][travis]
[![codecov][codecov-master]][codecov]


You can install the released version from [CRAN][cran-url].

```R
install.packages("reda")
```


## Development

[![Build Status][travis-dev]][travis]
[![codecov][codecov-dev]][codecov]

The latest version of package is under development at [GitHub][github-url] in
branch 'dev'.  If it is able to pass the building check by Travis CI, you may
consider installing it with the help of **remotes** by

```R
if (! require(remotes)) install.packages("remotes")
remotes::install_github("wenjie2wang/reda", ref = "dev")
```


## Get Started

- [Package vignette][vignette] provides a quick demonstration for the basic
  usage of main functions.

- [Package help manual][pdf-manual] is also available for more technical
  details.


## License

The R package reda is free software: You can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or any later version (at
your option).  See the [GNU General Public License][gpl-url] for details.

The R package reda is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.


[r-pkg-badge]: http://www.r-pkg.org/badges/version/reda
[cranlog-badge]: https://cranlogs.r-pkg.org/badges/splines2
[cran-url]: https://CRAN.R-project.org/package=reda
[travis]: https://travis-ci.org/wenjie2wang/reda
[travis-master]: https://travis-ci.org/wenjie2wang/reda.svg?branch=master
[travis-dev]: https://travis-ci.org/wenjie2wang/reda.svg?branch=dev
[github-url]: https://github.com/wenjie2wang/reda
[vignette]: https://wenjie-stat.me/reda/
[pdf-manual]: https://wenjie-stat.me/reda/reda.pdf
[gpl-url]: http://www.gnu.org/licenses/
[codecov]: https://codecov.io/gh/wenjie2wang/reda
[codecov-master]: https://codecov.io/gh/wenjie2wang/reda/branch/master/graph/badge.svg
[codecov-dev]: https://codecov.io/gh/wenjie2wang/reda/branch/dev/graph/badge.svg
