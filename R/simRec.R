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


##' Generate Simulated Survival Data or Recurrent Events for One Process
##'
##' The Cox proportional hazard model (Cox) Andersen
##' The thinning method (Lewis and Shedler, 1979) is applied for generating
##' simulated recurrent event times from Gamma frailty model for one process or
##' usually one subject or machine in applications. Conditional on predictors
##' (or covariates) and the unobserved Gamma frailty effect, the process is a
##' Poisson process.
##'
##' For argument \code{z}, \code{zCoef}, and \code{rho}, a function of time can
##' be specified for time-varying effect.  The (first) argument of the input
##' function has to be the time (not need to be named as "time" though). Other
##' arguments of the function can be specified through a named list in
##' \code{arguments}.
##'
##' @param z Time-invariant or time-varying covariates. The default value is
##'     \code{0} for no covariate effect.  This argument should be a numeric
##'     vector for time-invariant covariates or a function of time that returns
##'     a numeric vector for time-varying covariates.
##' @param zCoef Time-invariant or time-varying coefficients of covariates. The
##'     default value is \code{1}. Similar to the argument \code{x}, this
##'     argument should be a numeric vector for time-invariant coefficients or a
##'     function of time that returns a numeric vector for time-varying
##'     coefficients. The length of the numeric vector specified or returned
##'     from \code{x} and \code{x.coef} has to be always the same.
##' @param rho Baseline rate (or intensity) function for the Poisson process.
##'     The default is \code{1} for a homogenous Poisson process of unit
##'     intensity. This argument can be either a non-negative numeric value for
##'     a homogenous Poisson process or a function of time for a non-homogenous
##'     process.
##' @param rhoCoef Coefficients of baseline rate function. The default value is
##'     \code{1}. It can be useful when \code{rho} is a function generating
##'     spline bases.
##' @param origin The time origin set to be \code{0} by default.
##' @param endTime The end of follow-up time set to be \code{5} by default.
##' @param frailty Frailty effect.
##'
##' @param arguments Other arguments that can be specified through a named list
##'     for those time-varying functions.
##' @param ... Other arguemtns for future usage.
##'
##'
##' @references
##'
##' Cox, D. R. (1972). Regression models and life-tables.
##' \emph{Journal of the Royal Statistical Society. Series B
##' (Methodological)}, 34(2), 187--220.
##'
##' Lewis, P. A., & G. S. Shedler. (1979). Simulation of
##' Nonhomogeneous Poisson Processes by Thinning.
##' \emph{Naval Research Logistics Quarterly},
##' 26(3), Wiley Online Library: 403--13.
##'
##' @examples
##' library(reda)
##' set.seed(1216)
##'
##' ### example of customized time-varying baseline rate function
##' ## the baseline rate function has to be non-negative
##' simRec(rho = function(timeVec) { sin(timeVec) + 1 })
##' ## other arguments can be specified through a named list in 'arguments'
##' simRec(rho = function(a, b) { cos(a + b) + 1 },
##'        arguments = list(rho = list(b = 1)))
##' ## use quadratic M-splines without one internal knot
##' simRec(rho = "mSpline", rhoCoef = c(0.2, 0.5, 0.3, 0.4),
##'        arguments = list(rho = list(degree = 2, df = 4, intercept = TRUE,
##'                         Boundary.knots = c(0, 5))))
##'
##' ### example of time-invariant covariates and coefficients
##' simRec(z = 1, zCoef = 0.5)
##' simRec(z = c(0.5, 1), zCoef = c(1, 0))
##'
##' ### example of time-varying covariates and time-varying coefficients
##' zFun <- function(timeVec, intercept) {
##'    c(timeVec / 10 + intercept, as.numeric(timeVec > 3))
##' }
##' zCoefFun <- function(timeVec, shift) {
##'   c(sqrt(timeVec + shift), 1)
##' }
##' simRec(z = zFun, zCoef = zCoefFun,
##'        arguments = list(z = list(intercept = 0.1),
##'                         zCoef = list(shift = 0.1)))
##'
##' ### example of frailty effect
##' ## The default distribution is Gamma distribution
##' simRec(z = c(0.5, 1), zCoef = c(1, 0), frailty = TRUE,
##'        arguments = list(frailty = list(shape = 2, scale = 0.5)))
##' ## equivalent to the following function call
##' simRec(z = c(0.5, 1), zCoef = c(1, 0), frailty = "rgamma",
##'        arguments = list(frailty = list(shape = 2, scale = 0.5)))
##' ## lognormal with mean zero
##' logNorm <- function(a) exp(rnorm(n = 1, mean = 0, sd = a))
##' simRec(z = c(0.5, 1), zCoef = c(1, 0), frailty = logNorm,
##'        arguments = list(frailty = list(a = 1)))
##'
##' @importFrom stats optimize rexp qexp
simRec <- function(z = 0, zCoef = 1, rho = 1, rhoCoef = 1,
                   origin = 0, endTime = 5, frailty = FALSE,
                   recurrent = TRUE,
                   arguments = list(x = list(),
                                    xCoef = list(),
                                    rho = list(),
                                    frailty = list()),
                   ...)
{
    ## simple internal functions
    isNumVector <- function(x) is.numeric(x) & is.vector(x)
    isCharOne <- function(x) {
        is.character(x) & is.vector(x) & identical(length(x), 1L)
    }

    ## quick checks
    zVecIdx <- isNumVector(z)
    if (! (zVecIdx | is.function(z) | isCharOne(z)))
        stop(wrapMessages(
            "'z' has to be a numeric vector, a function or a function name"
        ))
    zCoefVecIdx <- isNumVector(zCoef)
    if (! (zCoefVecIdx | is.function(zCoef) | isCharOne(zCoef)))
        stop(wrapMessages(
            "'zCoef' has to be a numeric vector, a function or a function name"
        ))
    rhoVecIdx <- isNumVector(rho) & identical(length(rho), 1L)
    if (! (is.function(rho) | rhoVecIdx | isCharOne(rho)))
        stop("'rho' has to be a numeric vector or a function")
    if (rhoVecIdx & rho < 0)
        stop("The baseline rate function 'rho' has to be non-negative.")

    ## gpet arguments
    z_args <- arguments$z
    zCoef_args <- arguments$zCoef
    rho_args <- arguments$rho
    frailty_args <- arguments$frailty

    ## covariate: time-varying or time-invariant
    zFun <- ifelse(zVecIdx, function(zVec) z, z)
    zCoefFun <- ifelse(zCoefVecIdx, function(zVec) zCoef, zCoef)
    rhoFun <- ifelse(rhoVecIdx, function(zVec) rho, rho)

    ## prepare frailty effect
    noFrailty <- FALSE
    if (is.logical(frailty)) {
        if (frailty) {
            ## default as Gamma distribution
            frailty <- "rgamma"
        } else {
            noFrailty <- TRUE
            frailtyEffect <- 1
        }
    }
    if (is.function(frailty) |
        (is.character(frailty) & identical(length(frailty), 1L))) {
        argsNames <- names(as.list(args(frailty)))
        ## add "n = 1" for common distribution from stats library
        if ("n" %in% argsNames) {
            frailty_args <- c(list(n = 1),
                              frailty_args[names(frailty_args) != "n"])
        }
        frailtyEffect <- do.call(frailty, frailty_args)
    } else if (! noFrailty) {
        stop(wrapMessages(
            "The argument 'frailty' has to be a logical value (TRUE or FALSE),",
            "a function (or its name)."
        ))
    }

    ## rate function
    rateFun <- function(timeNum) {
        zVec <- do.call(zFun, c(list(timeNum), z_args))
        zCoefVec <- do.call(zCoefFun, c(list(timeNum), zCoef_args))
        covEffect <- as.numeric(exp(zVec %*% zCoefVec))
        rhoMat <- do.call(rhoFun, c(list(timeNum), rho_args))
        rhoVec <- as.numeric(rhoMat %*% rhoCoef)
        out <- frailtyEffect * rhoVec * covEffect
        out
    }

    ## step 1: calculate the supremum value of rate function
    rhoMaxObj <- stats::optimize(rateFun, interval = c(origin, endTime),
                                 maximum = TRUE)
    rho_max <- rhoMaxObj$objective

    ## step 2: generate W_i in batch for possible better performance
    ## estimate the number of W_i before censoring
    if (recurrent) {
        batchNum <- ceiling((endTime - origin) / stats::qexp(0.10, rho_max))
        eventTime <- NULL
        lastEventTime <- origin
        while (lastEventTime < endTime) {
            W <- stats::rexp(n = batchNum, rate = rho_max)
            ## step 3: update evnet times
            eventTime <- c(eventTime, lastEventTime + cumsum(W))
            lastEventTime <- eventTime[length(eventTime)]
        }
    } else {
        eventTime <- origin + stats::rexp(n = 1, rate = rho_max)
    }
    ## only keep event time before end time
    eventTime <- eventTime[eventTime <= endTime]
    len_eventTime <- length(eventTime)
    if (! len_eventTime) {
        ## no any event
        xOut <- numeric(0)
    } else {
        ## step 4: thinning
        rho_t <- sapply(eventTime, rateFun)
        U <- runif(n = len_eventTime)
        ind <- U <= rho_t / rho_max
        xOut <- eventTime[ind]
    }
    ## return
    out <- list(x = xOut,
                z = z,
                zCoef = zCoef,
                rho = rho,
                rhoCoef = rhoCoef,
                origin = origin,
                endTime = endTime,
                frailty = frailty,
                frailtyEffect = frailtyEffect,
                recurrent = recurrent,
                arguments = arguments)
    class(out) <- "simRec"
    out
}


##' Generate dataset for a given number of subjects
simRecData <- function(nProcess,
                       x = 0,
                       x.coef = 1,
                       origin = 0, endTime = 5,
                       rho.fun = function(x) 1,
                       rho.args = list(),
                       rho.coef = 1,
                       frailty.fun = rgamma,
                       frailty.args = list(shape = 0.5, scale = 2),
                       ...)
{
    xMat <- as.matrix(x)
    ## determine number of process from x
    if (missing(nProcess)) {
        nProcess <- nrow(xMat)
    } else {
        ## prepare x, x.coef, origin, and endTime
        xMat <- if (nrow(xMat) < nProcess) {
                    apply(xMat, 2L, function(oneCol) {
                        rep(oneCol, length.out = nProcess)
                    })
                }
    }
    x.coef <- rep(x.coef, length.out = ncol(xMat))
    origin <- rep(origin, length.out = nProcess)
    endTime <- rep(endTime, length.out = nProcess)
    ## function convert results from simRec to data.frame
    simRec2data <- function(ID, obj) {
        timeVec <- as.numeric(obj)
        attrList <- attr(obj, "settings")
        xVec <- attrList$x
        ## if no event
        if (! length(obj))
            return(data.frame(ID = ID,
                              origin = attrList$origin,
                              time = attrList$endTime,
                              event = 0,
                              as.list(xVec)))
        ## else for any event
        data.frame(ID = ID,
                   origin = attrList$origin,
                   time = c(timeVec, attrList$endTime),
                   event = c(rep(1, length(timeVec)), 0),
                   as.list(xVec))
    }
    ## generate simulated data for each process
    resList <- lapply(seq_len(nProcess), function(i) {
        res <- simRec(x = xMat[i, ],
                      x.coef = x.coef,
                      origin = origin[i],
                      endTime = endTime[i],
                      rho.fun = rho.fun,
                      rho.args = rho.args,
                      rho.coef = rho.coef,
                      frailty.fun = frailty.fun,
                      frailty.args = frailty.args,
                      ...
                      )
        simRec2data(ID = i, res)
    })
    ## return
    do.call(rbind, resList)
}

foo <- function(a) {exp(rnorm(n = 1, mean = 0, sd = a))}
simRec(z = c(0.5, 1), zCoef = c(1, 0),
       frailty = "foo",
       arguments = list(frailty = list(a = 1)))
