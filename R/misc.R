################################################################################
##
##   R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
##   Copyright (C) 2015-2017
##
##   This file is part of the R package reda.
##
##   The R package reda is free software: You can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   any later version (at your option). See the GNU General Public License
##   at <http://www.gnu.org/licenses/> for details.
##
##   The R package reda is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##
################################################################################


### some trivial internal functions ============================================
## wrap messages and keep proper line length
wrapMessages <- function(..., strwrap.args = list()) {
    x <- paste(...)
    wrap_x <- do.call(strwrap, c(list(x = x), strwrap.args))
    paste(wrap_x, collapse = "\n")
}

## is x a numeric vector
isNumVector <- function(x) {
    is.numeric(x) && is.vector(x)
}

## is x a numeric value
isNumOne <- function(x) {
    isNumVector(x) && identical(length(x), 1L)
}

## is x a character vector
isCharVector <- function(x) {
    is.character(x) && is.vector(x)
}

## is x a character value
isCharOne <- function(x) {
    isCharVector(x) && identical(length(x), 1L)
}

## is x a Survr object
is.Survr <- function(x) {
    inherits(x, "Survr")
}

## is x a sampleMcf object
is.sampleMcf <- function(x) {
    inherits(x, "sampleMcf")
}

## is x a rateReg object
is.rateReg <- function(x) {
    inherits(x, "rateReg")
}
