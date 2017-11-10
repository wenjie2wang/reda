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
##' The function \code{simEve} generates simulated recurrent events or survival
##' time (the first event time) from one stochastic process. The function
##' \code{simEveData} provides a simple wrapper that calls \code{simEve}
##' internally and collects the generated survival data or recurrent events into
##' a data frame.
##'
##' For each process, a time-invariant or time-varying baseline hazard rate
##' (intensity) function of failure can be specified.  Covariates and their
##' coefficients can be specified and are incorporated based on the Cox
##' proportional hazard model (Cox, 1972) for survival data or Andersen-Gill
##' model (Andersen and Gill, 1982) for recurrent events. In addition, a frailty
##' effect can be considered.  Conditional on predictors (or covariates) and the
##' unobserved frailty effect, the process is by default a Poisson process,
##' where the interarrival times between two successive arrivals/events follow
##' exponential distribution. A general renewal process can be specified through
##' \code{interarrival} for other distributions of the interarrival times.
##'
##' The thinning method (Lewis and Shedler, 1979) is applied for bounded hazard
##' rate function by default. The method based on inverse cumulative
##' distribution function (CDF) is also available for possibly unbounded but
##' integrable rate function over the given time period. The inverse CDF method
##' will be used when the rate function may go to infinite and a warning will be
##' thrown out if the thinning method is specified originally.
##'
##' For the covariates \code{z}, the covariate coefficients \code{zCoef}, and
##' the baseline hazard rate function \code{rho}, a function of time can be
##' specified for time-varying effect.  The first argument of the input function
##' has to be the time variable (not need to be named as "time" though). Other
##' arguments of the function can be specified through a named list in
##' \code{arguments}, while the first argument should not be specified.
##'
##' For the frailty effect \code{frailty}, the starting point \code{origin}, and
##' the end point of the process \code{endTime}, functions that generate random
##' numbers can be specified. An argument \code{n = 1} will be implicitly
##' specified if the function has an argument named \code{n}, which is designed
##' for those common functions generating random numbers from \code{stats}
##' package. Similar to \code{z}, \code{zCoef}, and \code{rho}, other arguments
##' of the function can be specified through a named list in \code{arguments}.
##'
##' For time-varying covariates, the function \code{simEveData} assumes
##' covariates can be observed only at event times and censoring times. Thus,
##' covariate values are returned only at these time points. If we want more
##' observed covariate values to be recorded, we may write a simple wrapper
##' function of \code{simEve} similar to \code{simEveData}.
##'
##' @aliases simEve
##'
##' @usage
##' simEve(z = 0, zCoef = 1, rho = 1, rhoCoef = 1, origin = 0, endTime = 3,
##'        frailty = FALSE, recurrent = TRUE, interarrival = "rexp",
##'        method = c("thinning", "inverse.cdf"), arguments = list(), ...)
##'
##' @param z Time-invariant or time-varying covariates. The default value is
##'     \code{0} for no covariate effect.  This argument should be a numeric
##'     vector for time-invariant covariates or a function of time that returns
##'     a numeric vector for time-varying covariates.
##' @param zCoef Time-invariant or time-varying coefficients of covariates. The
##'     default value is \code{1}. Similar to the argument \code{z}, this
##'     argument should be a numeric vector for time-invariant coefficients or a
##'     function of time that returns a numeric vector for time-varying
##'     coefficients. The length of the numeric vector specified or returned
##'     from \code{z} and \code{zCoef} has to be always the same.
##' @param rho Baseline rate (or intensity) function for the Poisson process.
##'     The default is \code{1} for a homogenous process of unit intensity. This
##'     argument can be either a non-negative numeric value for a homogenous
##'     process or a function of time for a non-homogenous process.
##' @param rhoCoef Coefficients of baseline rate function. The default value is
##'     \code{1}. It can be useful when \code{rho} is a function generating
##'     spline bases.
##' @param origin The time origin set to be \code{0} by default. It should be
##'     either a numeric value less than \code{endTime} or a function that
##'     returns such a numeric value.
##' @param endTime The end of follow-up time set to be \code{3} by default.
##'     Similar to \code{origin}, \code{endTime} should be either a numeric
##'     value greater than \code{origin} or a function that returns such a
##'     numeric value.
##' @param frailty Frailty effect. An optional logical value indicating whether
##'     to consider a frailty model or a function that produces the frailty
##'     effect.  The default value is \code{FALSE} for no frailty effect. If
##'     \code{TRUE}, a frailty factor from Gamma distribution will be used and
##'     the shape and scale parameter has to be specified through a list named
##'     \code{frailty} in \code{arguments}. Similar to \code{z}, \code{zCoef},
##'     and \code{rho}, a function or a function name can be specified for other
##'     distribution of the frailty effect. The specified function should
##'     randomly return a positive numeric value. For example, the functions
##'     that generate random numbers following a certain distribution from
##'     \code{stats} package can directly used. All the arguments of the
##'     function can be specified through a list named \code{frailty} in
##'     \code{arguments}.
##' @param recurrent A logical value with default value \code{TRUE} indicating
##'     whether to generate recurrent event data or survival data (i.e. the
##'     first event only).
##' @param interarrival A function object (or a function name) for randomly
##'     generating (positive) interarrival time between two successive
##'     arrivals/events.  The default value is \code{"rexp"} for generating
##'     interarrival times following exponential distribution, which leads to a
##'     Poisson process. If the assumption of exponential interarrival times
##'     cannot be justified, we may consider a renewal process, (a
##'     generalization of Poisson process), in which interarrival times between
##'     events independently follows an identical distribution. A customized
##'     function can be specified in this case. It must have at least one
##'     argument named \code{rate} for the expected number of arrivals/events in
##'     unit time and returns one positive numerical value. If the function
##'     contains an argument named \code{n}, it is assumed that the function
##'     returns \code{n} interarrival times in one function call to possibly
##'     speed up the random number generation procedure.  Other arguments can be
##'     specified through a named list inside \code{arguments}.
##' @param method A character string specifying the method for generating
##'     simulated recurrent or survival data. The default method is thinning
##'     method (Lewis and Shedler, 1979). Another available option is the method
##'     based on inverse cumulative distribution function (CDF). When the rate
##'     function may go to infinite, the inverse CDF method is used and a
##'     warning will be thrown out if the thinning method is initially
##'     specified.
##' @param arguments A list that consists of named lists for specifying other
##'     arguments in the corresponding functions. For example, if a function of
##'     time named \code{foo} with two arguments, \code{x} (for time) and
##'     \code{y} is specified for the time-varying covariates, the value of its
##'     second argument \code{y} can be specified by letting \code{arguments =
##'     list(z = list(y = 1)}.  A partial matching on names is not allowed to
##'     avoid possible misspecification. The input arguments will be evaluated
##'     within function \code{simEve}, which can be useful for randomly setting
##'     function parameters for each process in function \code{simEveData}. See
##'     examples and vignettes for details.
##' @param ... Other arguemtns for future usage.
##'
##' @return The function \code{simEve} returns a \code{simEve} S4 class object
##'     and the function \code{simEveData} returns a \code{data.frame}.
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
##'
##' ### time-invariant covariates and coefficients
##' ## one process
##' simEve(z = c(0.5, 1), zCoef = c(1, 0))
##' simEve(z = 1, zCoef = 0.5, recurrent = FALSE)
##'
##' ## simulated data
##' simEveData(1, z = c(0.5, 1), zCoef = c(1, 0), endTime = 2)
##' simEveData(3, z = cbind(rnorm(3), 1), zCoef = c(1, 0))
##' simEveData(5, z = matrix(rnorm(5)), zCoef = 0.5, recurrent = FALSE)
##'
##'
##' ### time-varying covariates and time-varying coefficients
##' zFun <- function(time, intercept) {
##'    c(time / 10 + intercept, as.numeric(time > 1))
##' }
##' zCoefFun <- function(x, shift) {
##'   c(sqrt(x + shift), 1)
##' }
##' simEve(z = zFun, zCoef = zCoefFun,
##'        arguments = list(z = list(intercept = 0.1),
##'                         zCoef = list(shift = 0.1)))
##'
##' ## same function of time for all processes
##' simEveData(3, z = zFun, zCoef = zCoefFun,
##'            arguments = list(z = list(intercept = 0.1),
##'                             zCoef = list(shift = 0.1)))
##'
##' ## same function within one process but different between processes
##' ## use quote function in the arguments
##' simDat <- simEveData(3, z = zFun, zCoef = zCoefFun,
##'                      arguments = list(
##'                          z = list(intercept = quote(rnorm(1) / 10)),
##'                          zCoef = list(shift = 0.1)
##'                     ))
##' ## check the intercept randomly generated,
##' ## which should be the same within each ID but different between IDs.
##' unique(with(simDat, cbind(ID, intercept = round(X.1 - time / 10, 3))))
##'
##'
##' ### non-negative time-varying baseline hazard rate function
##' simEve(rho = function(timeVec) { sin(timeVec) + 1 })
##' simEveData(3, origin = rnorm(3), endTime = rnorm(3, 5),
##'            rho = function(timeVec) { sin(timeVec) + 1 })
##' ## specify other arguments
##' simEve(rho = function(a, b) { cos(a + b) + 1 },
##'        arguments = list(rho = list(b = 1)))
##' simEveData(z = cbind(rnorm(3), rbinom(3, 1, 0.5)),
##'            rho = function(a, b) { cos(a + b) + 1 },
##'            arguments = list(rho = list(b = 1)))
##'
##' ## quadratic B-splines with one internal knot at "time = 1"
##' ## (using function 'bSpline' from splines2 package)
##' simEve(rho = "bSpline", rhoCoef = c(0.4, 0.5, 0.3, 0.6),
##'        arguments = list(rho = list(degree = 2, knots = 1, intercept = TRUE,
##'                                    Boundary.knots = c(0, 3))))
##'
##'
##' ### frailty effect
##' ## The default distribution is Gamma distribution
##' set.seed(123)
##' simEve(z = c(0.5, 1), zCoef = c(1, 0), frailty = TRUE,
##'        arguments = list(frailty = list(shape = 2, scale = 0.5)))
##' ## equivalent to the following function call
##' set.seed(123)
##' simEve(z = c(0.5, 1), zCoef = c(1, 0), frailty = "rgamma",
##'        arguments = list(frailty = list(shape = 2, scale = 0.5)))
##'
##' ## lognormal with mean zero (on the log scale)
##' set.seed(123)
##' logNorm <- function(a) exp(rnorm(n = 1, mean = 0, sd = a))
##' simEve(z = c(0.5, 1), zCoef = c(1, 0), frailty = logNorm,
##'        arguments = list(frailty = list(a = 1)))
##' ## equivalent to the following function call
##' set.seed(123)
##' simEve(z = c(0.5, 1), zCoef = c(1, 0), frailty = "rlnorm",
##'        arguments = list(frailty = list(sdlog = 1)))
##'
##'
##' ### renewal process
##' ## interarrival times following uniform distribution
##' rUnif <- function(n, rate, min) runif(n, min, max = 2 / rate - min)
##' simEve(interarrival = rUnif,
##'        arguments = list(interarrival = list(min = 0.1)))
##'
##' ## interarrival times following Gamma distribution with scale one
##' set.seed(123)
##' simEve(interarrival = function(n, rate) rgamma(n, shape = 1 / rate))
##' ## or equivalently
##' set.seed(123)
##' simEve(interarrival = function(rate) rgamma(n = 1, shape = 1 / rate))
##'
##' @importFrom stats integrate optimize qexp rexp runif rgamma rpois uniroot
##' @export
simEve <- function(z = 0, zCoef = 1,
                   rho = 1, rhoCoef = 1,
                   origin = 0, endTime = 3,
                   frailty = FALSE,
                   recurrent = TRUE,
                   interarrival = "rexp",
                   method = c("thinning", "inverse.cdf"),
                   arguments = list(), ...)
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
            "a function or a function name."
        ))
    ## check coefficients zCoef
    zCoefVecIdx <- isNumVector(zCoef)
    if (! (zCoefVecIdx || is.function(zCoef) || isCharOne(zCoef)))
        stop(wrapMessages(
            "The covariate coefficients 'zCoef' has to be a numeric vector,",
            "a function or a function name."
        ))
    ## check baseline rate function rho
    rhoVecIdx <- isNumOne(rho)
    if (! (rhoVecIdx || is.function(rho) || isCharOne(rho)))
        stop("The baseline hazard rate function",
             "'rho' has to be a numeric vector, a function or a function name.")
    if (rhoVecIdx && rho < 0)
        stop("The baseline hazard rate function 'rho' has to be non-negative.")
    ## check function for interarrival time
    if (! (is.function(interarrival) || isCharOne(interarrival)))
        stop("The 'interarrival' has to be a function or a function name.")
    interarrivalFun <- if (isCharOne(interarrival)) {
                           eval(parse(text = interarrival))
                       } else {
                           interarrival
                       }
    defaultIntArvIdx <- missing(interarrival) ||
        identical(interarrivalFun, stats::rexp)

    ## get arguments
    z_args <- lapply(arguments[["z"]], eval)
    zCoef_args <- lapply(arguments[["zCoef"]], eval)
    rho_args <- lapply(arguments[["rho"]], eval)

    origin_args <- lapply(arguments[["origin"]], eval)
    origin_args <- origin_args[names(origin_args) != "n"]

    endTime_args <- lapply(arguments[["endTime"]], eval)
    endTime_args <- endTime_args[names(endTime_args) != "n"]

    frailty_args <- lapply(arguments[["frailty"]], eval)
    frailty_args <- frailty_args[names(frailty_args) != "n"]

    interarrival_args <- lapply(arguments[["interarrival"]], eval)
    interarrival_args <-
        interarrival_args[! names(interarrival_args) %in% c("n", "rate")]
    intArvArgs <- names(as.list(args(interarrival)))
    if (! "rate" %in% intArvArgs)
        stop(wrapMessages(
            "The function for interarrival times must have",
            "one arguments named 'rate' for the expected number of",
            "events/arrivals in unit time."
        ))

    ## check origin and endTime
    if (is.function(origin) || isCharOne(origin)) {
        ## add "n = 1" for common distribution from stats library
        if ("n" %in% names(as.list(args(origin))))
            origin_args <- c(list(n = 1), origin_args)
        originFun <- origin
        origin_args <- if (! length(origin_args)) list()
        origin <- do.call(originFun, origin_args)
    } else {
        originFun <- origin_args <- NULL
    }
    if (is.function(endTime)|| isCharOne(endTime)) {
        ## add "n = 1" for common distribution from stats library
        if ("n" %in% names(as.list(args(endTime))))
            endTime_args <- c(list(n = 1), endTime_args)
        endTimeFun <- endTime
        endTime_args <- if (! length(endTime_args)) list()
        endTime <- do.call(endTimeFun, endTime_args)
    } else {
        endTimeFun <- endTime_args <- NULL
    }
    if (! (isNumOne(origin) && isNumOne(endTime) &&
           origin < endTime && is.finite(endTime))) {
        stop(wrapMessages(
            "The 'origin' and 'endTime'",
            "has to be two numerical values s.t. 'origin' < 'endTime < Inf'."
        ))
    }

    ## covariate: time-varying or time-invariant
    zFun <- if (zVecIdx) function(x) { z } else z
    zCoefFun <- if (zCoefVecIdx) function(x) { zCoef } else zCoef
    rhoFun <- if (rhoVecIdx) function(x) { rho } else rho

    ## prepare frailty effect
    if (is.logical(frailty)) {
        ## default as Gamma distribution
        frailty <- if(frailty) "rgamma" else NULL
    }
    if (is.function(frailty) || isCharOne(frailty)) {
        ## add "n = 1" for common distribution from stats library
        if ("n" %in% names(as.list(args(frailty))))
            frailty_args <- c(list(n = 1), frailty_args)
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
            "a function (or a function name)."
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
    if (is.infinite(rho_max)) {
        if (identical(method, "thinning")) {
            method <- "inverse.cdf"
            warning(wrapMessages(
                "The rate function may go to infinite.",
                "The Inverse CDF method was used."
            ))
        }
    }

    ## values at end time (censoring time)
    cenList <- rateFun(endTime, forOptimize = FALSE)
    zMat_cen <- cenList$zVec
    zCoefMat_cen <- cenList$zCoefVec
    rhoMat_cen <- cenList$rhoMat

    ## thinning method
    if (identical(method, "thinning")) {
        ## step 2: generate W_i in batch for possible better performance
        ## take care of possible interarrival arguments
        interarrivalArgs <- c(list(rate = rho_max), interarrival_args)
        if (recurrent) {
            ## estimate the number of W_i before censoring
            batchNum <- ceiling((endTime - origin) / stats::qexp(0.20, rho_max))
            if ("n" %in% intArvArgs)
                interarrivalArgs <- c(list(n = batchNum), interarrivalArgs)
            eventTime <- NULL
            lastEventTime <- origin
            while (lastEventTime < endTime) {
                W <- do.call(interarrival, interarrivalArgs)
                if (! isNumVector(W) || any(W <= 0))
                    stop("The interarrival times must be positive!")
                ## step 3: update evnet times
                eventTime <- c(eventTime, lastEventTime + cumsum(W))
                lastEventTime <- eventTime[length(eventTime)]
            }
        } else {
            if ("n" %in% intArvArgs)
                interarrivalArgs <- c(list(n = 1), interarrivalArgs)
            W <- do.call(interarrival, interarrivalArgs)
            if (! isNumOne(W) || any(W <= 0))
                stop("The interarrival time must be a positive nuumber!")
            eventTime <- origin + W
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
        ## error if rate function is not integrable
        if ("error" %in% class(intRate))
            stop(wrapMessages(
                "The integral of rate function",
                "is probably divergent."
            ))
        ## determine number of events, numEvent
        if (defaultIntArvIdx) {
            numEvent <- stats::rpois(n = 1, lambda = intRate)
        } else {
            ## take care of possible interarrival arguments
            interarrivalArgs <- c(list(rate = intRate), interarrival_args)
            ## take a larger step when argument 'n' is available
            if ("n" %in% intArvArgs)
                interarrivalArgs <- c(list(n = 10), interarrivalArgs)
            ## initialize for the while loop
            numEvent <- 0L
            lastEventTime <- origin
            while (lastEventTime < endTime) {
                W <- do.call(interarrival, interarrivalArgs)
                if (any(W <= 0))
                    stop("The interarrival times must be positive!")
                ## step 3: update evnet times
                eventTime <- lastEventTime + cumsum(W)
                numEvent <- numEvent + sum(eventTime < endTime)
                lastEventTime <- eventTime[length(eventTime)]
            }
        }
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
    methods::new("simEve", xOut,
                 call = Call,
                 z = list(
                     z = zMat,
                     fun = zFun,
                     args = zArgs,
                     timevarying = ! zVecIdx
                 ),
                 zCoef = list(
                     zCoef = zCoefMat,
                     fun = zCoefFun,
                     args = zArgs,
                     timevarying = ! zCoefVecIdx
                 ),
                 rho = list(
                     rho = rhoMat,
                     fun = rhoFun,
                     args = rhoArgs,
                     timevarying = ! rhoVecIdx
                 ),
                 rhoCoef = rhoCoef,
                 frailty = list(
                     frailty = frailtyEffect,
                     fun = frailty,
                     args = frailtyArgs
                 ),
                 origin = list(
                     origin = origin,
                     fun = originFun,
                     args = origin_args
                 ),
                 endTime = list(
                     endTime = endTime,
                     fun = endTimeFun,
                     args = endTime_args
                 ),
                 censoring = list(
                     z = zMat_cen,
                     zCoef = zCoefMat_cen,
                     rho = rhoMat_cen
                 ),
                 recurrent = recurrent,
                 interarrival = list(
                     fun = interarrival,
                     args = interarrival_args
                 ),
                 method = method
                 )

}


##' @rdname simEve
##' @aliases simEveData
##' @usage
##' simEveData(nProcess = 1, z = 0, zCoef = 1, rho = 1, rhoCoef = 1,
##'            origin = 0, endTime = 3, frailty = FALSE, recurrent = TRUE,
##'            interarrival = "rexp", method = c("thinning", "inverse.cdf"),
##'            arguments = list(), ...)
##'
##' @param nProcess Number of stochastic processes. A positive number should be
##'     speicified. The default value is \code{1}.
##'
##' @export
simEveData <- function(nProcess = 1,
                       z = 0, zCoef = 1,
                       rho = 1, rhoCoef = 1,
                       origin = 0, endTime = 3,
                       frailty = FALSE,
                       recurrent = TRUE,
                       interarrival = "rexp",
                       method = c("thinning", "inverse.cdf"),
                       arguments = list(),
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
            "a function or a function name."
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

    ## take care of origin and endTime before simEve
    originFunIdx <- is.function(origin) || isCharOne(origin)
    if (! originFunIdx) {
        origin <- rep(origin, length.out = nProcess)
    }
    endTimeFunIdx <- is.function(endTime) || isCharOne(endTime)
    if (! endTimeFunIdx) {
        endTime <- rep(endTime, length.out = nProcess)
    }

    ## generate simulated data for each process
    resList <- lapply(seq_len(nProcess), function(i) {
        res <- simEve(z = if (isZmatIdx) z[i, ] else z,
                      zCoef = zCoef,
                      rho = rho,
                      rhoCoef = rhoCoef,
                      origin = if (originFunIdx) origin else origin[i],
                      endTime = if (endTimeFunIdx) endTime else endTime[i],
                      frailty = frailty,
                      recurrent = recurrent,
                      interarrival = interarrival,
                      method = method,
                      arguments = arguments,
                      ...)
        simEve2data(ID = i, res)
    })

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
## function convert results from simEve to data.frame
simEve2data <- function(ID, obj) {
    timeVec <- obj@.Data
    ## if no event
    if (! length(timeVec))
        return(data.frame(ID = ID,
                          time = obj@endTime$endTime,
                          event = 0,
                          origin = obj@origin$origin,
                          X = obj@censoring$z))
    ## else for any event
    data.frame(ID = ID,
               time = c(timeVec, obj@endTime$endTime),
               event = c(rep(1, length(timeVec)), 0),
               origin = obj@origin$origin,
               X = rbind(obj@z$z, obj@censoring$z))
}
