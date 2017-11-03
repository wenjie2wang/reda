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


## collation after class.R
##' @include class.R
NULL


##' Simulated Survival times or Recurrent Events
##'
##' The function \code{simRec} generates simulated recurrent events or survival
##' time (the first event time) from one stochastic process. The function
##' \code{simRecData} generates simulated recurrent event data or survival data
##' with the help of the function \code{simRec}.
##'
##' For each process, a time-invariant or time-varying baseline hazard rate
##' (intensity) function of failure can be specified.  Covariates and their
##' coefficients can be specified and are incorporated based on the Cox
##' proportional hazard model (Cox, 1972) for survival data or Andersen-Gill
##' model (Andersen and Gill, 1982) for recurrent events. In addition, a frailty
##' effect can be considered.  Conditional on predictors (or covariates) and the
##' unobserved frailty effect, the process is a Poisson process.
##'
##' The thinning method (Lewis and Shedler, 1979) is applied for bounded hazard
##' rate function by default. The method based on inverse cumulative
##' distribution function (CDF) is also available for possibly unbounded but
##' integrable rate function over the given time period.
##'
##' For covariates \code{z}, covariate coefficients \code{zCoef}, and baseline
##' hazard rate function \code{rho}, a function of time can be specified for
##' time-varying effect.  The (first) argument of the input function has to be
##' the time (not need to be named as "time" though). Other arguments of the
##' function can be specified through a named list in \code{arguments}.
##'
##' @aliases simRec
##'
##' @usage
##' simRec(z = 0, zCoef = 1, rho = 1, rhoCoef = 1, origin = 0,
##'        endTime = 3, frailty = FALSE, recurrent = TRUE,
##'        method = c("thinning", "inverse.cdf"),
##'        arguments = list(x = list(), xCoef = list(),
##'                         rho = list(), frailty = list()), ...)
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
##'     from \code{z} and \code{zCoef} has to be always the same.
##' @param rho Baseline rate (or intensity) function for the Poisson process.
##'     The default is \code{1} for a homogenous Poisson process of unit
##'     intensity. This argument can be either a non-negative numeric value for
##'     a homogenous Poisson process or a function of time for a non-homogenous
##'     process.
##' @param rhoCoef Coefficients of baseline rate function. The default value is
##'     \code{1}. It can be useful when \code{rho} is a function generating
##'     spline bases.
##' @param origin The time origin set to be \code{0} by default.
##' @param endTime The end of follow-up time set to be \code{3} by default.
##' @param frailty Frailty effect. An optional logical value indicating whether
##'     to consider a frailty model or a function that produces the frailty
##'     effect.  The default value is \code{FALSE} for no frailty effect. If
##'     \code{TRUE}, a frailty factor from Gamma distribution will be used and
##'     the shape and scale parameter has to be specified through a list named
##'     \code{frailty} in \code{arguments}. Similar to \code{z}, \code{zCoef},
##'     and \code{rho}, a function or its name can be specified for other
##'     distribution of the frailty effect. The specified function should
##'     randomly return a positive numeric value. For example, the functions
##'     that generate random numbers following a certain distribution from
##'     \code{stats} package can directly used. Again, all the arguments of the
##'     function can be specified through a list named \code{frailty} in
##'     \code{arguments}. Note that \code{n = 1} will be implicitly specified if
##'     the function has an argument named \code{n}, which is designed for those
##'     common functions generating random numbers from \code{stats} package.
##' @param recurrent A logical value with default value \code{TRUE} indicating
##'     whether to generate recurrent event data or survival data (i.e. the
##'     first event only).
##' @param method A character string specifying the method for generating
##'     simulated recurrent or survival data. The default method is thinning
##'     method (Lewis and Shedler, 1979). Another available option is the method
##'     based on inverse cumulative distribution function (CDF). When the rate
##'     function may go to infinite, the inverse CDF method is used and a
##'     warning will be thrown out if the thinning method is initially
##'     specified.
##' @param arguments Other arguments that can be specified through a named list
##'     for those time-varying functions. The input arguments will be evaluated
##'     within function \code{simRec} before feed into the function, which can
##'     be useful for randomly setting function parameters for each process in
##'     function \code{simRecData}.
##' @param ... Other arguemtns for future usage.
##'
##' @return The function \code{simRec} returns a \code{simRec} S4 class object
##'     and the function \code{simRecData} returns a \code{data.frame}.
##'
##' @references
##'
##' Andersen, P. K., & Gill, R. D. (1982). Cox's regression model for counting
##' processes: A large sample study. \emph{The annals of statistics}, 10(4),
##' 1100--1120.
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
##' ### time-invariant covariates and coefficients
##' ## one process
##' simRec(z = c(0.5, 1), zCoef = c(1, 0))
##' simRec(z = 1, zCoef = 0.5, recurrent = FALSE)
##'
##' ## simulated data
##' simRecData(1, z = c(0.5, 1), zCoef = c(1, 0), endTime = 2)
##' simRecData(3, z = cbind(rnorm(3), 1), zCoef = c(1, 0))
##' simRecData(5, z = matrix(rnorm(5)), zCoef = 0.5, recurrent = FALSE)
##'
##' ### time-varying covariates and time-varying coefficients
##' zFun <- function(timeVec, intercept) {
##'    c(timeVec / 10 + intercept, as.numeric(timeVec > 1))
##' }
##' zCoefFun <- function(timeVec, shift) {
##'   c(sqrt(timeVec + shift), 1)
##' }
##' simRec(z = zFun, zCoef = zCoefFun,
##'        arguments = list(z = list(intercept = 0.1),
##'                         zCoef = list(shift = 0.1)))
##'
##' ## same function of time for all processes
##' simRecData(3, z = zFun, zCoef = zCoefFun,
##'            arguments = list(z = list(intercept = 0.1),
##'                             zCoef = list(shift = 0.1)))
##'
##' ## same function within one process but different between processes
##' ## use quote function in the arguments
##' simDat <- simRecData(3, z = zFun, zCoef = zCoefFun,
##'                      arguments = list(
##'                          z = list(intercept = quote(rnorm(1) / 10)),
##'                          zCoef = list(shift = 0.1)
##'                     ))
##' ## check the intercept randomly generated,
##' ## which should be the same within each ID but different between IDs.
##' unique(with(simDat, cbind(ID, intercept = round(X.1 - time / 10, 3))))
##'
##' ### non-negative time-varying baseline hazard rate function
##' simRec(rho = function(timeVec) { sin(timeVec) + 1 })
##' simRecData(3, origin = rnorm(3), endTime = rnorm(3, 5),
##'            rho = function(timeVec) { sin(timeVec) + 1 })
##' ## specify other arguments
##' simRec(rho = function(a, b) { cos(a + b) + 1 },
##'        arguments = list(rho = list(b = 1)))
##' simRecData(z = cbind(rnorm(3), rbinom(3, 1, 0.5)),
##'            rho = function(a, b) { cos(a + b) + 1 },
##'            arguments = list(rho = list(b = 1)))
##'
##' ## quadratic I-splines with one internal knot at "time = 1"
##' ## (using function 'iSpline' from splines2 package)
##' simRec(rho = "iSpline", rhoCoef = c(0.2, 0.5, 0.3, 0.4),
##'        arguments = list(rho = list(degree = 2, knots = 1, intercept = TRUE,
##'                                    Boundary.knots = c(0, 3))))
##'
##' ### frailty effect
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
##' @importFrom stats integrate optimize qexp rexp runif rgamma rpois uniroot
##' @export
simRec <- function(z = 0, zCoef = 1, rho = 1, rhoCoef = 1,
                   origin = 0, endTime = 3,
                   frailty = FALSE, recurrent = TRUE,
                   method = c("thinning", "inverse.cdf"),
                   arguments = list(x = list(),
                                    xCoef = list(),
                                    rho = list(),
                                    frailty = list()),
                   ...)
{
    ## record function call
    Call <- match.call()
    ## match method
    method <- match.arg(method)

    ## check covariate z
    zVecIdx <- isNumVector(z)
    if (! (zVecIdx || is.function(z) || isCharOne(z)))
        stop(wrapMessages(
            "The covariates 'z' has to be a numeric vector / matrix,",
            "a function or its name."
        ))
    ## check coefficients zCoef
    zCoefVecIdx <- isNumVector(zCoef)
    if (! (zCoefVecIdx || is.function(zCoef) || isCharOne(zCoef)))
        stop(wrapMessages(
            "The covariate coefficients 'zCoef' has to be a numeric vector,",
            "a function or its name."
        ))
    ## check baseline rate function rho
    rhoVecIdx <- isNumOne(rho)
    if (! (rhoVecIdx || is.function(rho) || isCharOne(rho)))
        stop("'rho' has to be a numeric vector or a function")
    if (rhoVecIdx && rho < 0)
        stop("The baseline hazard rate function 'rho' has to be non-negative.")
    ## check origin and endTime
    if (! (isNumOne(origin) && isNumOne(endTime) &&
           origin < endTime) && is.finite(endTime))
        stop(wrapMessages(
            "The 'origin' and 'endTime'",
            "has to be two numerical values s.t. 'origin' < 'endTime < Inf'."
        ))

    ## get arguments
    z_args <- lapply(arguments[["z"]], eval)
    zCoef_args <- lapply(arguments[["zCoef"]], eval)
    rho_args <- lapply(arguments[["rho"]], eval)
    frailty_args <- lapply(arguments[["frailty"]], eval)

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

    ## values at end time (censoring time)
    cenList <- rateFun(endTime, forOptimize = FALSE)
    zMat_cen <- cenList$zVec
    zCoefMat_cen <- cenList$zCoefVec
    rhoMat_cen <- cenList$rhoMat

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
            zMat <- zMat_cen
            zCoefMat <- zCoefMat_cen
            rhoMat <- rhoMat_cen
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
        if (identical(numEvent, 0)) {
            xOut <- numeric(0)
            ## only return values on end time
            zMat <- zMat_cen
            zCoefMat <- zCoefMat_cen
            rhoMat <- rhoMat_cen
        } else {
            U <- sort(runif(n = numEvent))
            if (! recurrent)
                U <- U[1L]
            ## density function may go to infinite
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
            xOut <- sapply(U, invFun)
            ## compute zMat, zCoefMat, and rhoMat
            resList <- lapply(xOut, rateFun, forOptimize = FALSE)
            zMat <- do.call(rbind, lapply(resList, function(a) a$zVec))
            zCoefMat <- do.call(rbind, lapply(resList, function(a) a$zCoefVec))
            rhoMat <- do.call(rbind, lapply(resList, function(a) a$rhoMat))
        }
    }

    ## prepare outputs
    ## for covariates
    if (zVecIdx) {
        zFun <- zArgs <- NULL
    } else {
        zFun <- z
        zArgs <- if (length(z_args)) z_args else list()
    }
    ## for covariate coefficients
    if (zCoefVecIdx) {
        zCoefFun <- zArgs <- NULL
    } else {
        zCoefFun <- zCoef
        zArgs <- if (length(zCoef_args)) zCoef_args else list()
    }
    ## for baseline rate function
    if (rhoVecIdx) {
        rhoFun <- rhoArgs <- NULL
    } else {
        rhoFun <- rho
        rhoArgs <- if (length(rho_args)) rho_args else list()
    }
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
                     z = zMat,
                     zFun = zFun,
                     zArgs = zArgs,
                     zTimeVarying = ! zVecIdx
                 ),
                 zCoef = list(
                     zCoef = zCoefMat,
                     zCoefFun = zCoefFun,
                     zArgs = zArgs,
                     zCoefTimeVarying = ! zCoefVecIdx
                 ),
                 rho = list(
                     rho = rhoMat,
                     rhoFun = rhoFun,
                     rhoArgs = rhoArgs,
                     rhoTimeVarying = ! rhoVecIdx
                 ),
                 rhoCoef = rhoCoef,
                 frailty = list(
                     frailtyEffect = frailtyEffect,
                     frailtyFun = frailty,
                     frailtyArgs = frailtyArgs
                 ),
                 origin = origin,
                 endTime = endTime,
                 censoring = list(
                     z = zMat_cen,
                     zCoef = zCoefMat_cen,
                     rho = rhoMat_cen
                 ),
                 recurrent = recurrent,
                 method = method
                 )

}


##' @rdname simRec
##' @aliases simRecData
##' @usage
##' simRecData(nProcess = 1, z = 0, zCoef = 1, rho = 1, rhoCoef = 1,
##'            origin = 0, endTime = 3, frailty = FALSE, recurrent = TRUE,
##'            method = c("thinning", "inverse.cdf"),
##'            arguments = list(z = list(), zCoef = list(),
##'                             rho = list(), frailty = list()), ...)
##'
##' @param nProcess Number of stochastic processes. A positive number should be
##'     speicified. The default value is \code{1}.
##'
##' @export
simRecData <- function(nProcess = 1, z = 0, zCoef = 1,
                       rho = 1, rhoCoef = 1,
                       origin = 0, endTime = 3,
                       frailty = FALSE, recurrent = TRUE,
                       method = c("thinning", "inverse.cdf"),
                       arguments = list(z = list(),
                                        zCoef = list(),
                                        rho = list(),
                                        frailty = list()),
                       ...)
{
    ## record function call
    Call <- match.call()
    ## match method
    method <- match.arg(method)

    ## check covariate z
    if (isNumVector(z)) z <- matrix(z, nrow = 1)  # convert vector z to matrix
    if (! (is.matrix(z) || is.function(z) || isCharOne(z)))
        stop(wrapMessages(
            "The covariates 'z' has to be a numeric vector / matrix,",
            "a function or its name."
        ))

    ## if covariates are given as a matrix
    if (isZmatIdx <- is.matrix(z)) {
        if (missing(nProcess)) {
            ## determine number of process from z
            nProcess <- nrow(z)
        } else {
            ## recycle z if necessary
            if (nrow(z) < nProcess) {
                z <- apply(z, 2L, function(oneCol) {
                    rep(oneCol, length.out = nProcess)
                })
            }
        }
        zCoef <- rep(zCoef, length.out = ncol(z))
    }

    origin <- rep(origin, length.out = nProcess)
    endTime <- rep(endTime, length.out = nProcess)

    ## generate simulated data for each process
    resList <- if (isZmatIdx) {
                   lapply(seq_len(nProcess), function(i) {
                       res <- simRec(z = z[i, ],
                                     zCoef = zCoef,
                                     rho = rho,
                                     rhoCoef = rhoCoef,
                                     origin = origin[i],
                                     endTime = endTime[i],
                                     frailty = frailty,
                                     recurrent = recurrent,
                                     method = method,
                                     arguments = arguments,
                                     ...)
                       simRec2data(ID = i, res)
                   })
               } else {
                   lapply(seq_len(nProcess), function(i) {
                       res <- simRec(z = z,
                                     zCoef = zCoef,
                                     rho = rho,
                                     rhoCoef = rhoCoef,
                                     origin = origin[i],
                                     endTime = endTime[i],
                                     frailty = frailty,
                                     recurrent = recurrent,
                                     method = method,
                                     arguments = arguments,
                                     ...)
                       simRec2data(ID = i, res)
                   })
               }
    ## prepare for output
    out <- do.call(rbind, resList)
    if (! recurrent) {
        uniIdx <- duplicated(out$ID, fromLast = TRUE)
        out <- out[uniIdx, ]
    }
    ## return
    out
}


### internal functions =========================================================
## function convert results from simRec to data.frame
simRec2data <- function(ID, obj) {
    timeVec <- obj@.Data
    ## if no event
    if (! length(timeVec))
        return(data.frame(ID = ID,
                          time = obj@endTime,
                          event = 0,
                          origin = obj@origin,
                          X = obj@censoring$z))
    ## else for any event
    data.frame(ID = ID,
               time = c(timeVec, obj@endTime),
               event = c(rep(1, length(timeVec)), 0),
               origin = obj@origin,
               X = rbind(obj@z$z, obj@censoring$z))
}
