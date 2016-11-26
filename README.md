# reda

The R package **reda** provides function fitting gamma frailty model with
baseline rate function modeled by regression splines for recurrent event data,
as well as functions computing and ploting the mean cumulative funciton from
sample (Nelson-Aalen estimator) and the fitted frailty model.


## Installation of CRAN Version

[![CRAN_Status_Badge][1]][2]
[![Build Status][3]][4]

You can install the released version from [CRAN][2].

```R
install.packages("reda")
```


## Development

[![Build Status][5]][4]

The latest version of package is under development at [GitHub][6] in branch
'dev'.  If it is able to pass the building check by Travis CI, you may consider
installing it with the help of **devtools** by

```R
devtools::install_github("wenjie2wang/reda", ref = "dev")
```

or cloning this reposotory to local and install by makefile as follows:

```
git clone https://github.com/wenjie2wang/reda.git
cd reda
make install
```


## Get Started

- [Package vignettes][7] provides a quick demonstration for the basic usage of
  main functions.

- [Package help manual][8] is also available for more technical details.


## License

The R package reda is free software: You can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or any later version (at
your option).  See the [GNU General Public License][9] for details.

The R package reda is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.


[1]: http://www.r-pkg.org/badges/version/reda
[2]: https://CRAN.R-project.org/package=reda
[3]: https://travis-ci.org/wenjie2wang/reda.svg?branch=master
[4]: https://travis-ci.org/wenjie2wang/reda
[5]: https://travis-ci.org/wenjie2wang/reda.svg?branch=dev
[6]: https://github.com/wenjie2wang/reda
[7]: http://wenjie-stat.me/reda/
[8]: http://wenjie-stat.me/reda/reda.pdf
[9]: http://www.gnu.org/licenses/
