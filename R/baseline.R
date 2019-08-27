##
## R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
## Copyright (C) 2015-2019
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


## collation after class.R
##' @include class.R
NULL


##' Estimated Baseline Rate Function
##'
##' An S4 class generic function that returns the estimated baseline rate
##' function.
##'
##' @aliases baseRate
##'
##' @param object An object used to dispatch a method.
##' @param ... Other arguments for future usage.
##'
##' @return A \code{baseRate} object.
##'
##' @examples
##' ## See examples given in function rateReg.
##' @seealso
##' \code{\link{rateReg}} for model fitting;
##' \code{\link{summary,rateReg-method}} for summary of a fitted model;
##' \code{\link{plot,baseRate.rateReg-method}} for ploting method.
##' @export
setGeneric(
    name = "baseRate",
    def = function(object, ...) {
        standardGeneric("baseRate")
    }
)


##' @describeIn baseRate Estiamted baseline rate from a fitted model.
##'
##' @param level An optional numeric value
##' indicating the confidence level required. The default value is 0.95.
##' @param control An optional list to specify the time grid
##' where the baseline rate function is estimated.
##' The availble elements of the control list include
##' \code{grid}, \code{length.out}, \code{from} and \code{to}.
##' The time grid can be directly specified via element \code{grid}.
##' A dense time grid is suggested.
##' Element \code{length.out} represents the length of grid points.
##' The dafault value is 1,000.
##' Element \code{from} means the starting point of grid with default 0.
##' Element \code{to} represnts the endpoint of grid
##' with the right boundary knot as default.
##' When \code{grid} is missing, the grid will be generated
##' by \code{seq} (from package \pkg{base})
##' with arguments \code{from}, \code{to} and \code{length.out}.
##'
##' @aliases baseRate,rateReg-method
##'
##' @importFrom stats qnorm
##'
##' @export
setMethod(
    f = "baseRate",
    signature = "rateReg",
    definition = function(object, level = 0.95, control = list(), ...) {

        ## baseline rate coefficients
        alpha <- object@estimates$alpha[, 1L, drop = FALSE]

        ## time grid
        splinesList <- object@spline
        knots <- splinesList$knots
        degree <- splinesList$degree
        Boundary.knots <- splinesList$Boundary.knots
        controlist <- c(control, list(Boundary.knots_ = Boundary.knots))
        control <- do.call("rateReg_mcf_control", controlist)
        gridTime <- control$grid

        ## remove the grid on the right boundary for piece-wise constant bases?
        ## gridTime <- gridTime[gridTime < Boundary.knots[2L]]

        ## reconstruct spline basis matrix
        bList <- list(x = gridTime, knots = knots, degree = degree,
                      intercept = TRUE, Boundary.knots = Boundary.knots)
        splineName <- splinesList$spline
        bMat <- if (splineName == "bSplines") {
                    do.call(splines2::bSpline, bList)
                } else if (splineName == "mSplines") {
                    do.call(splines2::mSpline, bList)
                } else {
                    stop("Unknown splines type.")
                }

        ## point estimates
        estVec <- as.numeric(bMat %*% alpha)

        ## variance-covariance matrix
        nBeta <- nrow(object@estimates$beta)
        ind <- seq_len(nBeta + 1L)
        covMat <- solve(object@fisher)[- ind, - ind, drop = FALSE]
        varVec <- apply(bMat, 1L, function (bVec, covMat) {
            crossprod(bVec, covMat) %*% bVec
        }, covMat = covMat)

        seVec <- tryCatch(sqrt(varVec), warning = function(w) w)
        if ("warning" %in% class(seVec)) {
            stop(wrapMessages(
                "The variance-covariance matrix is not positive definite.",
                "Please check possible error",
                "(or adjust spline bases and perhaps",
                "try different set of starting values)."
            ))
        }

        ## confidence interval for the given level
        confBand <- stats::qnorm((1 + level) / 2) * seVec
        lower <- estVec - confBand
        upper <- estVec + confBand

        ## prepare for output
        outDat <- data.frame(time = gridTime, baseRate = estVec,
                             se = seVec, lower = lower, upper = upper)
        methods::new("baseRate.rateReg",
                     baseRate = outDat,
                     level = level)
    })
