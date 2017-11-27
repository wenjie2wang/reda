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


##' Recurrent Event Data Analysis
##'
##' The R package \pkg{reda} provides functions for simulating, exploring and
##' modeling recurrent event data.
##'
##' The main functions are summarized as follows:
##' \itemize{
##'
##' \item \code{simEventData}: Simulating survival, recurrent event, and
##' multiple event data from stochastic process point of view.
##'
##' \item \code{mcf}: Estimating the mean cumulative function (MCF) from a
##' fitted gamma frailty model, or from a sample recurrent event data by using
##' the nonparametic MCF estimator (the Nelson-Aelen estimator of the cumulative
##' hazard function).
##'
##' \item \code{mcfDiff}: Comparing two-sample MCFs by the pseudo-score tests
##' and estimating their difference over time.
##'
##' \item \code{rateReg}: Fitting Gamma fraitly model with spline baseline rate
##' function.
##' }
##'
##' See the package vignettes for more introduction and demonstration.
##'
##' @importFrom methods is new setClass setGeneric setMethod validObject
##' @docType package
##' @name reda-package
NULL
