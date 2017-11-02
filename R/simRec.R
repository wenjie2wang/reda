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
##' ### example of time-invariant covariates and coefficients
##' simRec(z = c(0.5, 1), zCoef = c(1, 0))
##' simRec(z = 1, zCoef = 0.5, recurrent = FALSE)
##'
##' ### example of customized time-varying baseline rate function
##' ## the baseline rate function has to be non-negative
##' simRec(rho = function(timeVec) { sin(timeVec) + 1 })
##' ## other arguments can be specified through a named list in 'arguments'
##' simRec(rho = function(a, b) { cos(a + b) + 1 },
##'        arguments = list(rho = list(b = 1)))
##' ## quadratic I-splines with one internal knot
##' ## (using function 'iSpline' from splines2 package)
##' simRec(rho = "iSpline", rhoCoef = c(0.2, 0.5, 0.3, 0.4),
##'        arguments = list(rho = list(degree = 2, df = 4, intercept = TRUE,
##'                         Boundary.knots = c(0, 5))))
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
                   origin = 0, endTime = 5,
                   frailty = FALSE,
                   recurrent = TRUE,
                   arguments = list(x = list(),
                                    xCoef = list(),
                                    rho = list(),
                                    frailty = list()),
                   method = c("thinning", "inverse.cdf"),
                   ...)
{
    ## record function call
    Call <- match.call()
    ## match method
    method <- match.arg(method)

    ## some simple internal functions
    isNumVector <- function(x) is.numeric(x) && is.vector(x)
    isNumOne <- function(x) isNumVector(x) && identical(length(x), 1L)
    isCharOne <- function(x) {
        is.character(x) && is.vector(x) && identical(length(x), 1L)
    }

    ## check covariate z
    zVecIdx <- isNumVector(z)
    if (! (zVecIdx || is.function(z) || isCharOne(z)))
        stop(wrapMessages(
            "'z' has to be a numeric vector, a function or a function name"
        ))
    ## check coefficients zCoef
    zCoefVecIdx <- isNumVector(zCoef)
    if (! (zCoefVecIdx || is.function(zCoef) || isCharOne(zCoef)))
        stop(wrapMessages(
            "'zCoef' has to be a numeric vector, a function or a function name"
        ))
    ## check baseline rate function rho
    rhoVecIdx <- isNumOne(rho)
    if (! (rhoVecIdx || is.function(rho) || isCharOne(rho)))
        stop("'rho' has to be a numeric vector or a function")
    if (rhoVecIdx && rho < 0)
        stop("The baseline rate function 'rho' has to be non-negative.")
    ## check origin and endTime
    if (! (isNumOne(origin) && isNumOne(endTime) && origin < endTime))
        stop(wrapMessages(
            "The 'origin' and 'endTime'",
            "has to be two numerical values s.t. 'origin' < 'endTime'."
        ))

    ## get arguments
    z_args <- arguments$z
    zCoef_args <- arguments$zCoef
    rho_args <- arguments$rho
    frailty_args <- arguments$frailty

    ## covariate: time-varying or time-invariant
    zFun <- ifelse(zVecIdx, function(zVec) z, z)
    zCoefFun <- ifelse(zCoefVecIdx, function(zVec) zCoef, zCoef)
    rhoFun <- ifelse(rhoVecIdx, function(zVec) rho, rho)

    ## prepare frailty effect
    if (is.logical(frailty)) {
        ## default as Gamma distribution
        frailty <- if(frailty) "rgamma" else NULL
    }
    if (is.function(frailty) || isCharOne(frailty)) {
        argsNames <- names(as.list(args(frailty)))
        ## add "n = 1" for common distribution from stats library
        if ("n" %in% argsNames) {
            frailty_args <- c(list(n = 1),
                              frailty_args[names(frailty_args) != "n"])
        }
        frailtyEffect <- do.call(frailty, frailty_args)
        ## check frailty effect of length one numeric
        if (! isNumOne(frailtyEffect)) {
            stop(wrapMessages(
                "The Frailty effect (returned from function)",
                "must be of length one."
            ))
        }
    } else if (is.null(frailty)) {
        frailtyEffect <- 1
    } else {
        stop(wrapMessages(
            "The argument 'frailty' has to be a logical value (TRUE or FALSE),",
            "a function (or its name)."
        ))
    }

    ## rate function
    rateFun <- function(timeNum, forOptimize = TRUE) {
        zVec <- do.call(zFun, c(list(timeNum), z_args))
        zCoefVec <- do.call(zCoefFun, c(list(timeNum), zCoef_args))
        covEffect <- as.numeric(exp(zVec %*% zCoefVec))
        rhoMat <- do.call(rhoFun, c(list(timeNum), rho_args))
        rhoVec <- as.numeric(rhoMat %*% rhoCoef)
        rho_t <- frailtyEffect * rhoVec * covEffect
        if (forOptimize)
            return(rho_t)
        ## else return
        list(rho_t = rho_t,
             rhoMat = matrix(rhoMat, nrow = 1),
             zVec = matrix(zVec, nrow = 1),
             zCoefVec = matrix(zCoefVec, nrow = 1))
    }

    ## the baseline rate function has to be non-negative
    rhoMinObj <- stats::optimize(rateFun, interval = c(origin, endTime),
                                 maximum = FALSE)
    if (rhoMinObj$objective < 0)
        stop(wrapMessages(
            "The rate function has to be non-negative",
            "from 'origin' to 'endTime'."
        ))

    ## step 1: calculate the supremum value of rate function
    rhoMaxObj <- stats::optimize(rateFun, interval = c(origin, endTime),
                                 maximum = TRUE)
    rho_max <- rhoMaxObj$objective
    ## if the supremum is finite, use thinning method
    if (is.infinite(rho_max) && identical(method, "thinning")) {
        method <- "inverse.cdf"
        warning(wrapMessages(
            "The rate function may go to infinite.",
            "The Inverse CDF method was used."
        ))
    }

    ## thinning method
    if (identical(method, "thinning")) {
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
            resList <- lapply(eventTime, rateFun, forOptimize = FALSE)
            rho_t <- sapply(resList, function(a) a$rho_t)
            U <- runif(n = len_eventTime)
            ind <- U <= rho_t / rho_max
            xOut <- eventTime[ind]
        }
        ## update zMat, zCoefMat, and rhoMat
        if (length(xOut)) {
            zMat <- do.call(rbind, lapply(resList[ind], function(a) {
                a$zVec
            }))
            zCoefMat <- do.call(rbind, lapply(resList[ind], function(a) {
                a$zCoefVec
            }))
            rhoMat <- do.call(rbind, lapply(resList[ind], function(a) {
                a$rhoMat
            }))
        } else {
            ## only return values on end time
            resList <- rateFun(endTime, forOptimize = FALSE)
            zMat <- resList$zVec
            zCoefMat <- resList$zCoefVec
            rhoMat <- resList$rhoMat
        }

    } else {
        ## (naive) method based on inverse CDF
        vecRateFun <- Vectorize(rateFun)
        intRate <- tryCatch(
            stats::integrate(vecRateFun, lower = origin,
                             upper = endTime)$value,
            error = function(e) e)
        if ("error" %in% class(intRate))
            stop(wrapMessages(
                "The integral of rate function",
                "is probably divergent"
            ))
        numEvent <- stats::rpois(n = 1, lambda = intRate)
        if (! recurrent)
            numEvent <- min(numEvent, 1)
        if (! length(numEvent)) {
            xOut <- numeric(0)
            xVec <- endTime
        } else {
            U <- sort(runif(n = numEvent))
            ## denFun may still go to infinite
            ## hard to apply rejection sampling
            ## use inverse CDF method numerically / the hard way
            invFun <- function(prob) {
                foo <- function(timeNum) {
                    stats::integrate(vecRateFun, lower = origin,
                                         upper = timeNum)$value / intRate - prob
                }
                stats::uniroot(foo, interval = c(origin + .Machine$double.eps,
                                                 endTime))$root
            }
            xVec <- xOut <- sapply(U, invFun)
        }
        ## compute zMat, zCoefMat, and rhoMat
        resList <- lapply(xVec, rateFun, forOptimize = FALSE)
        zMat <- do.call(rbind, lapply(resList, function(a) a$zVec))
        zCoefMat <- do.call(rbind, lapply(resList, function(a) a$zCoefVec))
        rhoMat <- do.call(rbind, lapply(resList, function(a) a$rhoMat))
    }

    ## prepare outputs
    ## for covariates
    if (zVecIdx) {
        zFun <- zArgs <- NULL
    } else {
        zFun <- z
        zArgs <- if (length(z_args)) z_args else list()
    }
    z <- zMat
    ## for covariate coefficients
    if (zCoefVecIdx) {
        zCoefFun <- zArgs <- NULL
    } else {
        zCoefFun <- zCoef
        zArgs <- if (length(zCoef_args)) zCoef_args else list()
    }
    zCoef <- zCoefMat
    ## for baseline rate function
    if (rhoVecIdx) {
        rhoFun <- rhoArgs <- NULL
    } else {
        rhoFun <- rho
        rhoArgs <- if (length(rho_args)) rho_args else list()
    }
    rho <- rhoMat
    ## for frailty
    frailtyArgs <- if (is.null(frailty)) {
                       NULL
                   } else {
                       if (length(frailty_args)) frailty_args else list()
                   }

    ## return
    methods::new("simRec", xOut,
                 call = Call,
                 z = list(
                     z = z,
                     zFun = zFun,
                     zArgs = zArgs,
                     zTimeVarying = ! zVecIdx
                 ),
                 zCoef = list(
                     zCoef = zCoef,
                     zCoefFun = zCoefFun,
                     zArgs = zArgs,
                     zCoefTimeVarying = ! zCoefVecIdx
                 ),
                 rho = list(
                     rho = rho,
                     rhoFun = rhoFun,
                     rhoArgs = rhoArgs,
                     rhoTimeVarying = ! rhoVecIdx
                 ),
                 rhoCoef = rhoCoef,
                 origin = origin,
                 endTime = endTime,
                 frailty = list(
                     frailtyEffect = frailtyEffect,
                     frailtyFun = frailty,
                     frailtyArgs = frailtyArgs
                 ),
                 recurrent = recurrent,
                 method = method
                 )

}
