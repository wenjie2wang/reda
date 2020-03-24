##
## R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
## Copyright (C) 2015-2020
##
## This file is part of the R package reda.
##
## The R package reda is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package reda is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##

## set default options
reda_default_options <- list(
    reda.Recur.maxPrint = 3L
)

## set options for reda
.onLoad <- function(libname, pkgname) {
  op <- options()

  toset <- ! names(reda_default_options) %in% names(op)
  if (any(toset))
      options(reda_default_options[toset])

  invisible(NULL)
}
