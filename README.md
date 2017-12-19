# reda

## Overview

The R pacakge **reda** provides functions for

- simulating survival, recurrent event, and multiple event data from stochastic
  process point of view;
- exploring and modeling recurrent event data through the mean cumulative
  function (MCF) or also called the Nelson-Aalen estimator of the cumulative
  hazard rate function, and gamma frailty model with spline rate function;
- comparing two-sample recurrent event responses with the pseudo-score tests.


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


## Getting Started

- [Online documentation][homepage]
    - [Package vignette][reda-intro] on exploring and modeling recurrent event
      data.
    - [Package vignette][reda-simulate] on simulating survival and recurrent
      event data.


## License

The R package reda is free software: You can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or any later version (at
your option).  See the [GNU General Public License][gpl-url] for details.

The R package reda is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.


[r-pkg-badge]: https://www.r-pkg.org/badges/version/reda
[cranlog-badge]: https://cranlogs.r-pkg.org/badges/splines2
[cran-url]: https://CRAN.R-project.org/package=reda
[travis]: https://travis-ci.org/wenjie2wang/reda
[travis-master]: https://travis-ci.org/wenjie2wang/reda.svg?branch=master
[travis-dev]: https://travis-ci.org/wenjie2wang/reda.svg?branch=dev
[github-url]: https://github.com/wenjie2wang/reda
[homepage]: https://wenjie-stat.me/reda/
[reda-intro]: https://wenjie-stat.me/reda/articles/reda-intro.html
[reda-simulate]: https://wenjie-stat.me/reda/articles/reda-simulate.html
[gpl-url]: https://www.gnu.org/licenses/
[codecov]: https://codecov.io/gh/wenjie2wang/reda
[codecov-master]: https://codecov.io/gh/wenjie2wang/reda/branch/master/graph/badge.svg
[codecov-dev]: https://codecov.io/gh/wenjie2wang/reda/branch/dev/graph/badge.svg
