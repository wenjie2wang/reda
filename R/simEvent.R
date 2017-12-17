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
##' The function \code{simEvent} generates simulated recurrent events or
##' survival time (the first event time) from one stochastic process. The
##' function \code{simEventData} provides a simple wrapper that calls
##' \code{simEvent} internally and collects the generated survival data or
##' recurrent events into a data frame. More examples are available in one of
##' the package vignettes in addition to the function documentation.
##'
##' For each process, a time-invariant or time-varying baseline hazard rate
##' (intensity) function of failure can be specified.  Covariates and their
##' coefficients can be specified and incorporated by the specified relative
##' risk functions. The default is the exponential relative risk function, which
##' corresponds to the Cox proportional hazard model (Cox 1972) for survival
##' data or Andersen-Gill model (Andersen and Gill 1982) for recurrent
##' events. Other relative risk function can be specified through the argument
##' \code{relativeRisk}. In addition, a frailty effect can be considered.
##' Conditional on predictors (or covariates) and the unobserved frailty effect,
##' the process is by default a Poisson process, where the interarrival times
##' between two successive arrivals/events follow exponential distribution. A
##' general renewal process can be specified through \code{interarrival} for
##' other distributions of the interarrival times in addition to the exponential
##' distribution.
##'
##' The thinning method (Lewis and Shedler 1979) is applied for bounded hazard
##' rate function by default. The inversion method (Cinlar 1975) is also
##' available for possibly unbounded but integrable rate function over the given
##' time period. The inversion method will be used when the rate function may go
##' to infinite and a warning will be thrown out if the thinning method is
##' specified originally.
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
##' for those common functions generating random numbers from \pkg{stats}
##' package. Similar to \code{z}, \code{zCoef}, and \code{rho}, other arguments
##' of the function can be specified through a named list in \code{arguments}.
##'
##' For time-varying covariates, the function \code{simEventData} assumes
##' covariates can be observed only at event times and censoring times. Thus,
##' covariate values are returned only at these time points. If we want other
##' observed covariate values to be recorded, we may write a simple wrapper
##' function for \code{simEvent} similar to \code{simEventData}.
##'
##' @aliases simEvent
##'
##' @usage
##' simEvent(z = 0, zCoef = 1, rho = 1, rhoCoef = 1, origin = 0, endTime = 3,
##'          frailty = 1, recurrent = TRUE, interarrival = "rexp",
##'          relativeRisk = c("exponential", "linear", "excess", "none"),
##'          method = c("thinning", "inversion"), arguments = list(), ...)
##'
##' @param z Time-invariant or time-varying covariates. The default value is
##'     \code{0} for no covariate effect.  This argument should be a numeric
##'     vector for time-invariant covariates or a function of times that returns
##'     a numeric matrix for time-varying covariates, where each row represents
##'     the covariate vector at one perticular time point.
##' @param zCoef Time-invariant or time-varying coefficients of covariates. The
##'     default value is \code{1}. Similar to the argument \code{z}, this
##'     argument should be a numeric vector for time-invariant coefficients or a
##'     function of times that returns a numeric matrix for time-varying
##'     coefficients, where each row represents the coefficient vector at one
##'     perticular time point. The dimension of the \code{z} and \code{zCoef}
##'     (either specified or generated) has to match with each other.
##' @param rho Baseline rate (or intensity) function for the Poisson process.
##'     The default is \code{1} for a homogenous process of unit intensity. This
##'     argument can be either a non-negative numeric value for a homogenous
##'     process or a function of times for a non-homogenous process. In the
##'     latter case, the function should be able to take a vector of time points
##'     and return a numerical matrix (or vector) with each row representing the
##'     baseline hazard rate vector (or scalar value) at each time point.
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
##' @param frailty A positive number or a function for frailty effect. The
##'     default value is \code{1} for no frailty effect.  Other positive value
##'     can be specified directly for a shared frailty effect within a cluster.
##'     Similar to \code{z}, \code{zCoef}, and \code{rho}, a function can be
##'     specified for other distribution of the frailty effect. The specified
##'     function should randomly return a positive numeric value. The functions
##'     that generate random numbers following a certain distribution from
##'     \code{stats} package can be directly used. The arguments of the function
##'     can be specified through a list named \code{frailty} in
##'     \code{arguments}. For example, if we consider Gamma distribution with
##'     mean one as the distribution of frailty effect, we may specify
##'     \code{frailty = "rgamma"}. The shape and scale parameter needs to be
##'     specified through a list named \code{frailty} in \code{arguments}, such
##'     as \code{arguments = list(frailty = list(shape = 2, scale = 0.5))}.
##' @param recurrent A logical value with default value \code{TRUE} indicating
##'     whether to generate recurrent event data or survival data.
##' @param interarrival A function object for randomly generating (positive)
##'     interarrival time between two successive arrivals/events.  The default
##'     value is \code{"rexp"} (i.e., function \code{\link[stats]{rexp}}) for
##'     generating interarrival times following exponential distribution, which
##'     leads to a Poisson process. If the assumption of exponential
##'     interarrival times cannot be justified, we may consider a renewal
##'     process, (a generalization of Poisson process), in which interarrival
##'     times between events independently follows an identical distribution. A
##'     customized function can be specified in this case. It must have at least
##'     one argument named \code{rate} for the expected number of
##'     arrivals/events in unit time and returns one positive numerical
##'     value. If the function contains an argument named \code{n}, it is
##'     assumed that the function returns \code{n} interarrival times in one
##'     function call to possibly speed up the random number generation
##'     procedure.  Other arguments can be specified through a named list inside
##'     \code{arguments}.
##' @param relativeRisk Relateive risk function for incorporating the covariates
##'     and the covariate coefficients into the intensity function. The
##'     applicable choices include \code{exponential} (the default) for the
##'     regular Cox model or Andersen-Gill model, \code{linear} for linear model
##'     (including an intercept term), \code{excess} for excess model, and
##'     \code{none} for not incorporating the covariates through a relative risk
##'     function. A customized function can be specified. The specified function
##'     must have at least one argument named \code{z} for the covariate vector
##'     and another argument named {zCoef} for covariate coefficient vector.
##'     The function should return a numeric value for given \code{z} vector and
##'     \code{zCoef} vector.  Other arguments can be specified through a named
##'     list inside \code{arguments}.
##' @param method A character string specifying the method for generating
##'     simulated recurrent or survival data. The default method is thinning
##'     method (Lewis and Shedler 1979). Another available option is the
##'     inversion method (Cinlar 1975). When the rate function may go to
##'     infinite, the inversion method is used and a warning will be thrown out
##'     if the thinning method is initially specified.
##' @param arguments A list that consists of named lists for specifying other
##'     arguments in the corresponding functions. For example, if a function of
##'     time named \code{foo} with two arguments, \code{x} (for time) and
##'     \code{y}, is specified for the time-varying covariates, the value of its
##'     second argument, \code{y}, can be specified by \code{arguments = list(z
##'     = list(y = 1)}.  A partial matching on names is not allowed to avoid
##'     possible misspecification. The input arguments will be evaluated within
##'     function \code{simEvent}, which can be useful for randomly setting
##'     function parameters for each process in function
##'     \code{simEventData}. See examples and vignettes for details.
##' @param ... Additional arguements passed from function \code{simEventData} to
##'     fucntion \code{simEvent}. For function \code{simEvent}, \code{...} is
##'     not used.
##'
##' @return
##'
##' The function \code{simEvent} returns a \code{simEvent} S4 class object and
##' the function \code{simEventData} returns a \code{data.frame}.
##'
##' @references
##'
##' Andersen, P. K., & Gill, R. D. (1982). Cox's regression model for counting
##' processes: A large sample study. \emph{The annals of statistics}, 10(4),
##' 1100--1120.
##'
##' Cinlar, Erhan (1975). \emph{Introduction to stochastic processes}. Englewood
##' Cliffs, NJ: Printice-Hall.
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
##' set.seed(123)
##' ### time-invariant covariates and coefficients
##' ## one process
##' simEvent(z = c(0.5, 1), zCoef = c(1, 0))
##' simEvent(z = 1, zCoef = 0.5, recurrent = FALSE)
##'
##' ## simulated data
##' simEventData(z = c(0.5, 1), zCoef = c(1, 0), endTime = 2)
##' simEventData(z = cbind(rnorm(3), 1), zCoef = c(1, 0))
##' simEventData(z = matrix(rnorm(5)), zCoef = 0.5, recurrent = FALSE)
##'
##'
##' ### time-varying covariates and time-varying coefficients
##' zFun <- function(time, intercept) {
##'    cbind(time / 10 + intercept, as.numeric(time > 1))
##' }
##' zCoefFun <- function(x, shift) {
##'   cbind(sqrt(x + shift), 1)
##' }
##' simEvent(z = zFun, zCoef = zCoefFun,
##'          arguments = list(z = list(intercept = 0.1),
##'                           zCoef = list(shift = 0.1)))
##'
##' ## same function of time for all processes
##' simEventData(3, z = zFun, zCoef = zCoefFun,
##'              arguments = list(z = list(intercept = 0.1),
##'                               zCoef = list(shift = 0.1)))
##'
##' ## same function within one process but different between processes
##' ## use quote function in the arguments
##' simDat <- simEventData(3, z = zFun, zCoef = zCoefFun,
##'                        arguments = list(
##'                            z = list(intercept = quote(rnorm(1) / 10)),
##'                            zCoef = list(shift = 0.1)
##'                        ))
##' ## check the intercept randomly generated,
##' ## which should be the same within each ID but different between IDs.
##' unique(with(simDat, cbind(ID, intercept = round(X.1 - time / 10, 6))))
##'
##'
##' ### non-negative time-varying baseline hazard rate function
##' simEvent(rho = function(timeVec) { sin(timeVec) + 1 })
##' simEventData(3, origin = rnorm(3), endTime = rnorm(3, 5),
##'              rho = function(timeVec) { sin(timeVec) + 1 })
##' ## specify other arguments
##' simEvent(z = c(rnorm(1), rbinom(1, 1, 0.5)) / 10,
##'          rho = function(a, b) { sin(a + b) + 1 },
##'          arguments = list(rho = list(b = 0.5)))
##' simEventData(z = cbind(rnorm(3), rbinom(3, 1, 0.5)) / 10,
##'              rho = function(a, b) { sin(a + b) + 1 },
##'              arguments = list(rho = list(b = 0.5)))
##'
##' ## quadratic B-splines with one internal knot at "time = 1"
##' ## (using function 'bSpline' from splines2 package)
##' simEvent(rho = splines2::bSpline, rhoCoef = c(0.8, 0.5, 1, 0.6),
##'          arguments = list(rho = list(degree = 2, knots = 1,
##'                                      intercept = TRUE,
##'                                      Boundary.knots = c(0, 3))))
##'
##'
##' ### frailty effect
##' ## Gamma distribution with mean one
##' simEvent(z = c(0.5, 1), zCoef = c(1, 0), frailty = rgamma,
##'          arguments = list(frailty = list(shape = 2, scale = 0.5)))
##'
##' ## lognormal with mean zero (on the log scale)
##' set.seed(123)
##' simEvent(z = c(0.5, 1), zCoef = c(1, 0), frailty = "rlnorm",
##'          arguments = list(frailty = list(sdlog = 1)))
##' ## or equivalently
##' set.seed(123)
##' logNorm <- function(a) exp(rnorm(n = 1, mean = 0, sd = a))
##' simEvent(z = c(0.5, 1), zCoef = c(1, 0), frailty = logNorm,
##'          arguments = list(frailty = list(a = 1)))
##'
##' ### renewal process
##' ## interarrival times following uniform distribution
##' rUnif <- function(n, rate, min) runif(n, min, max = 2 / rate)
##' simEvent(interarrival = rUnif,
##'          arguments = list(interarrival = list(min = 0)))
##'
##' ## interarrival times following Gamma distribution with scale one
##' set.seed(123)
##' simEvent(interarrival = function(n, rate) rgamma(n, shape = 1 / rate))
##' ## or equivalently
##' set.seed(123)
##' simEvent(interarrival = function(rate) rgamma(n = 1, shape = 1 / rate))
##'
##' ### relative risk functioin
##' set.seed(123)
##' simEvent(relativeRisk = "linear")
##' ## or equivalently
##' rriskFun <- function(z, zCoef, intercept) {
##'     as.numeric(z %*% zCoef) + intercept
##' }
##' set.seed(123)
##' simEvent(relativeRisk = rriskFun,
##'          arguments = list(relativeRisk = list(intercept = 1)))
##'
##' @importFrom Rcpp sourceCpp
##' @useDynLib reda
##'
##' @importFrom stats integrate optimize qexp rexp runif rgamma rpois uniroot
##' @importFrom splines2 bSpline
##' @export
simEvent <- function(z = 0, zCoef = 1,
                     rho = 1, rhoCoef = 1,
                     origin = 0, endTime = 3,
                     frailty = 1,
                     recurrent = TRUE,
                     interarrival = "rexp",
                     relativeRisk = c("exponential", "linear",
                                      "excess", "none"),
                     method = c("thinning", "inversion"),
                     arguments = list(), ...)
{
    ## record function call
    Call <- match.call()
    ## match method
    method <- match.arg(method, c("thinning", "inversion"))

    ## check covariate z
    zVecIdx <- isNumVector(z)
    if (isCharOne(z)) z <- match.fun(z)
    if (! (zVecIdx || is.function(z)))
        stop(wrapMessages(
            "The covariates `z` has to be a numeric vector,",
            "or a function."
        ), call. = FALSE)
    ## check coefficients zCoef
    zCoefVecIdx <- isNumVector(zCoef)
    if (isCharOne(zCoef)) zCoef <- match.fun(zCoef)
    if (! (zCoefVecIdx || is.function(zCoef)))
        stop(wrapMessages(
            "The covariate coefficients `zCoef` has to be a numeric vector,",
            "a function."
        ), call. = FALSE)
    if (zVecIdx && zCoefVecIdx && length(zCoef) < (tmp <- length(z)))
        zCoef <- rep(zCoef, length.out = tmp)
    ## check baseline rate function rho
    rhoVecIdx <- isNumOne(rho)
    if (isCharOne(rho)) rho <- match.fun(rho)
    if (! (rhoVecIdx || is.function(rho)))
        stop(wrapMessages(
            "The baseline hazard rate function",
            "`rho` has to be a numeric number, a function."
        ), call. = FALSE)
    if (rhoVecIdx && rho < 0)
        stop("The baseline hazard rate function has to be non-negative.",
             call. = FALSE)
    ## check the baseline rate function coefficients, rhoCoef
    if (! isNumVector(rhoCoef))
        stop("The `rhoCoef` has to be a numeric vector", call. = FALSE)
    n_rhoCoef <- length(rhoCoef)
    ## check function for interarrival time
    if (isCharOne(interarrival)) interarrival <- match.fun(interarrival)
    if (! is.function(interarrival))
        stop("The `interarrival` has to be a function.", call. = FALSE)
    defaultIntArvIdx <- missing(interarrival) ||
        identical(interarrival, stats::rexp)
    ## match relative risk function
    rriskNames <- c("exponential", "linear", "excess", "none")
    rriskFun <- if (rriskFunIdx <- is.function(relativeRisk)) {
                    .vectorize_rrisk(relativeRisk)
                } else if (isCharVector(relativeRisk)) {
                    rriskInd <- pmatch(relativeRisk, rriskNames)[1L]
                    if (is.na(rriskInd))
                        stop(wrapMessages(
                            "The specified relative risk function",
                            "is not available."
                        ), call. = FALSE)
                    relativeRisk <- rriskNames[rriskInd]
                    ## strange function names that user may not create
                    switch(relativeRisk,
                           "exponential" = rrisk_exponential,
                           "linear" = rrisk_linear,
                           "excess" = rrisk_excess,
                           "none" = rrisk_none)
                } else {
                    stop(wrapMessages(
                        "The specified relative risk function `relativeRisk`",
                        "should be a function."
                    ), call. = FALSE)
                }

    ## get arguments
    z_args <- lapply(arguments[["z"]], eval)
    zCoef_args <- lapply(arguments[["zCoef"]], eval)
    rho_args <- lapply(arguments[["rho"]], eval)
    rho_args <- rho_args[! names(rho_args) %in% c("z", "zCoef")]

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
            "at least one argument named `rate` for the expected number of",
            "events/arrivals in unit time."
        ), call. = FALSE)

    rrisk_args <- lapply(arguments[["relativeRisk"]], eval)
    rrisk_args <- rrisk_args[! names(rrisk_args) %in% c("z", "zCoef")]
    rriskFunArgs <- names(as.list(args(rriskFun)))
    if (any(! c("z", "zCoef") %in% rriskFunArgs))
        stop(wrapMessages(
            "The relative risk function must have",
            "at least one argument named `z` for covariates and",
            "another argument named `zCoef` for covariate coefficients."
        ), call. = FALSE)

    ## check origin
    if (isCharOne(origin)) origin <- match.fun(origin)
    if (originFunIdx <- is.function(origin)) {
        ## add "n = 1" for common distribution from stats library
        if ("n" %in% names(as.list(args(origin))))
            origin_args <- c(list(n = 1), origin_args)
        originFun <- origin
        if (length(origin_args) == 0L) origin_args <- list()
        origin <- do.call(originFun, origin_args)
    } else {
        originFun <- origin_args <- NULL
    }
    ## check endTime similarly to origin
    if (isCharOne(endTime)) endTime <- match.fun(endTime)
    if (endTimeFunIdx <- is.function(endTime)) {
        ## add "n = 1" for common distribution from stats library
        if ("n" %in% names(as.list(args(endTime))))
            endTime_args <- c(list(n = 1), endTime_args)
        endTimeFun <- endTime
        if (length(endTime_args) == 0L) endTime_args <- list()
        endTime <- do.call(endTimeFun, endTime_args)
    } else {
        endTimeFun <- endTime_args <- NULL
    }
    ## check origin and endTime
    if (! (isNumOne(origin) && isNumOne(endTime) &&
           origin < endTime && is.finite(endTime))) {
        stop(wrapMessages(
            "The `origin` and `endTime`",
            "has to be two numerical values s.t. `origin` < `endTime` < `Inf`."
        ), call. = FALSE)
    }
    ## prepare frailty effect
    if (isCharOne(frailty)) frailty <- match.fun(frailty)
    if (isNumOne(frailty)) {
        frailtyEffect <- frailty
        frailtyFun <- NULL
    } else if (frailtyFunIdx <- is.function(frailty)) {
        ## add "n = 1" for common distribution from stats library
        if ("n" %in% names(as.list(args(frailty))))
            frailty_args <- c(list(n = 1), frailty_args)
        frailtyEffect <- do.call(frailty, frailty_args)
        frailtyFun <- frailty
    } else {
        frailtyEffect <- - 1            # leads to errors
    }
    if (! isNumOne(frailtyEffect) || frailtyEffect <= 0)
        stop(wrapMessages(
            "The argument `frailty` has to be a positive number or",
            "a function that generates a positive number."
        ), call. = FALSE)

    ## rate function
    ## takes a number of time points (including one point)
    ## returns rate function at those time points and zMat, zCoefMat, rhoMat
    rateFun <- function(timeVec, forOptimize = TRUE) {
        nTime <- length(timeVec)
        seqTime <- seq_len(nTime)
        repTime <- rep(1L, nTime)
        zMat <- if (zVecIdx) {
                    matrix(z, nrow = 1)[repTime, , drop = FALSE]
                } else {
                    matrix(do.call(z, c(list(timeVec), z_args)),
                           nrow = nTime)[seqTime, , drop = FALSE]
                }
        zCoefMat <- if (zCoefVecIdx) {
                        matrix(zCoef, nrow = 1)[repTime, , drop = FALSE]
                    } else {
                        matrix(do.call(zCoef, c(list(timeVec), zCoef_args)),
                               nrow = nTime)
                    }
        rhoMat <- if (rhoVecIdx) {
                      matrix(rho, nrow = 1)[repTime, , drop = FALSE]
                  } else {
                      tmpIdx <- c("z", "zCoef") %in% names(as.list(args(rho)))
                      matrix(do.call(
                          rho, c(list(timeVec),
                                 list(z = zMat, zCoef = zCoefMat)[tmpIdx],
                                 rho_args)
                      ), nrow = nTime)
                  }
        ## possibly improve performance for time-invariant z and zCoef
        covEffect <-
            if (zVecIdx && zCoefVecIdx) {
                rep(do.call(rriskFun,
                            c(list(z = zMat[1L, , drop = FALSE],
                                   zCoef = zCoefMat[1L, , drop = FALSE]),
                              rrisk_args)), nTime)
            } else {
                do.call(rriskFun,
                        c(list(z = zMat, zCoef = zCoefMat),
                          rrisk_args))
            }
        rhoVec <- as.numeric(rhoMat %*% rhoCoef)
        rho_t <- frailtyEffect * rhoVec * covEffect
        if (forOptimize)
            return(rho_t)
        ## check if it is not for optimize
        if (! isNumVector(covEffect) || any(covEffect <= 0))
            stop(wrapMessages(
                "The relative risk function should return positive values."
            ), call. = FALSE)
        ## the rate function has to be non-negative
        if (any(is.na(rho_t) || rho_t < 0))
            stop(wrapMessages(
                "The rate function has to be non-negative",
                "from `origin` to `endTime`."
            ), call. = FALSE)
        ## else return
        list(rho_t = rho_t,
             rhoMat = rhoMat,
             zMat = zMat,
             zCoefMat = zCoefMat)
    }

    if (method == "thinning") {
        ## step 1: calculate the supremum value of rate function
        rhoMaxObj <- tryCatch(
            stats::optim((origin + endTime) / 2, rateFun,
                         lower = origin, upper = endTime,
                         method = "L-BFGS-B",
                         control = list(fnscale = - 1)),
            error = function(e) e
        )
        ## if the supremum is finite, use thinning method
        if ("error" %in% class(rhoMaxObj)) {
            method <- "inversion"
            warning(wrapMessages(
                "The rate function may go to infinite.",
                "The Inversion method was used."
            ), call. = FALSE)
        } else {
            rho_max <- rhoMaxObj$value
        }
    }
    ## values at end time (censoring time)
    cenList <- rateFun(endTime, forOptimize = FALSE)
    zMat_cen <- cenList$zMat
    zCoefMat_cen <- cenList$zCoefMat
    rhoMat_cen <- cenList$rhoMat

    ## thinning method
    if (method == "thinning") {
        ## step 2: generate W_i in batch for possible better performance
        ## take care of possible interarrival arguments
        interarrivalArgs <- c(list(rate = rho_max), interarrival_args)
        eventTime <- NULL
        lastEventTime <- origin
        ## if we want all the recurrent event times
        if (recurrent) {
            ## estimate the number of W_i before censoring
            batchNum <- ceiling((endTime - origin) / stats::qexp(0.8, rho_max))
            if ("n" %in% intArvArgs)
                interarrivalArgs <- c(list(n = batchNum), interarrivalArgs)
            while (lastEventTime < endTime) {
                W <- do.call(interarrival, interarrivalArgs)
                if (! isNumVector(W) || any(W < 0))
                    stop("The interarrival times must be nonnegative!",
                         call. = FALSE)
                ## step 3: update evnet times
                eventTime <- c(eventTime, lastEventTime + cumsum(W))
                lastEventTime <- eventTime[length(eventTime)]
            }
            ## only keep event time before end time
            eventTime <- eventTime[eventTime <= endTime]
            len_eventTime <- length(eventTime)
            if (len_eventTime == 0L) {
                ## no any event
                xOut <- numeric(0)
            } else {
                ## step 4: thinning
                resList <- rateFun(eventTime, forOptimize = FALSE)
                rho_t <- resList$rho_t
                U <- runif(n = len_eventTime)
                ind <- U <= rho_t / rho_max
                xOut <- eventTime[ind]
            }
        } else {
            ## if only the first event is of interest
            ## we may break the loop once we get the first event
            batchNum <- 10
            if ("n" %in% intArvArgs)
                interarrivalArgs <- c(list(n = batchNum), interarrivalArgs)
            xOut <- numeric(0)
            while (lastEventTime < endTime) {
                W <- do.call(interarrival, interarrivalArgs)
                if (! isNumVector(W) || any(W < 0))
                    stop("The interarrival times must be nonnegative!",
                         call. = FALSE)
                ## step 3: update evnet times
                eventTime <- lastEventTime + cumsum(W)
                lastEventTime <- eventTime[batchNum]
                ## only keep event time before end time
                eventTime <- eventTime[eventTime <= endTime]
                len_eventTime <- length(eventTime)
                if (len_eventTime == 0L) {
                    ## no any event
                    break;
                } else {
                    ## step 4: thinning
                    resList <- rateFun(eventTime, forOptimize = FALSE)
                    rho_t <- resList$rho_t
                    U <- runif(n = len_eventTime)
                    ind <- U <= rho_t / rho_max
                    if (any(ind)) {
                        ind <- which(ind)[1L]
                        xOut <- eventTime[ind]
                        break;
                    }
                }
            }
        }

        ## update zMat, zCoefMat, and rhoMat
        if (length(xOut) > 0L) {
            zMat <- resList$zMat[ind, , drop = FALSE]
            zCoefMat <- resList$zCoefMat[ind, , drop = FALSE]
            rhoMat <- resList$rhoMat[ind, , drop = FALSE]
        } else {
            ## only return values on end time
            zMat <- zMat_cen
            zCoefMat <- zCoefMat_cen
            rhoMat <- rhoMat_cen
        }

    } else {
        ## the inversion method
        intRate <- tryCatch(
            stats::integrate(rateFun, lower = origin,
                             upper = endTime)$value,
            error = function(e) e)
        ## error if rate function is not integrable
        if ("error" %in% class(intRate))
            stop(wrapMessages(
                "The integral of rate function",
                "is probably divergent."
            ), call. = FALSE)
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
                if (any(W < 0))
                    stop("The interarrival times must be nonnegative!",
                         call. = FALSE)
                ## step 3: update evnet times
                eventTime <- lastEventTime + cumsum(W)
                numEvent <- numEvent + sum(eventTime < endTime)
                lastEventTime <- eventTime[length(eventTime)]
            }
        }
        if (numEvent == 0L) {
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
            ## use inversion method numerically / the hard way
            invFun <- function(prob) {
                foo <- function(timeNum) {
                    stats::integrate(rateFun, lower = origin,
                                     upper = timeNum)$value / intRate - prob
                }
                stats::uniroot(foo, interval = c(origin + .Machine$double.eps,
                                                 endTime))$root
            }
            xOut <- sapply(U, invFun)
            ## compute zMat, zCoefMat, and rhoMat
            resList <- rateFun(xOut, forOptimize = FALSE)
            zMat <- resList$zMat
            zCoefMat <- resList$zCoefMat
            rhoMat <- resList$rhoMat
        }
    }

    ## prepare outputs
    ## for covariates
    if (zVecIdx) {
        zFun <- zArgs <- NULL
    } else {
        zFun <- z
        zArgs <- if (length(z_args) > 0L) z_args else list()
    }
    ## for covariate coefficients
    if (zCoefVecIdx) {
        zCoefFun <- zArgs <- NULL
    } else {
        zCoefFun <- zCoef
        zArgs <- if (length(zCoef_args) > 0L) zCoef_args else list()
    }
    ## for baseline rate function
    if (rhoVecIdx) {
        rhoFun <- rhoArgs <- NULL
    } else {
        rhoFun <- rho
        rhoArgs <- if (length(rho_args) > 0L) rho_args else list()
    }
    ## for frailty
    frailtyArgs <- if (is.null(frailtyFun)) {
                       NULL
                   } else {
                       if (length(frailty_args) > 0) frailty_args else list()
                   }

    ## return
    methods::new("simEvent",
                 .Data = xOut,
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
                     fun = frailtyFun,
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
                 relativeRisk = list(
                     fun = relativeRisk,
                     args = rrisk_args
                 ),
                 method = method
                 )

}


##' @rdname simEvent
##' @aliases simEventData
##' @usage
##' simEventData(nProcess, z = 0, origin = 0, endTime = 3, frailty = 1, ...)
##'
##' @param nProcess Number of stochastic processes. If missing, the value will
##'     be the number of row of the specified matrix \code{z}. Otherwise, a
##'     positive number should be speicified.
##'
##' @export
simEventData <- function(nProcess = 1,
                         z = 0,
                         origin = 0,
                         endTime = 3,
                         frailty = 1,
                         ...)
{
    ## record function call
    Call <- match.call()

    ## check covariate z
    ## convert vector z to matrix
    if (isNumVector(z, error_na = TRUE)) z <- matrix(z, nrow = 1)
    if (! (is.matrix(z) || is.function(z)))
        stop(wrapMessages(
            "The covariates 'z' has to be a numeric vector / matrix,",
            "or a function."
        ))
    ## if covariates are given as a matrix
    if (zMatIdx <- isNumMatrix(z, error_na = TRUE)) {
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
    }

    ## take care of origin, endTime, and frailty before simEvent
    originFunIdx <- is.function(origin)
    if (! originFunIdx) {
        if (! isNumVector(origin, error_na = TRUE))
            stop("The time origins, `origin` has to be a numeric vector.")
        origin <- rep(origin, length.out = nProcess)
    }
    endTimeFunIdx <- is.function(endTime)
    if (! endTimeFunIdx) {
        if (! isNumVector(endTime, error_na = TRUE))
            stop("The ends of time, `endTime` has to be a numeric vector.")
        endTime <- rep(endTime, length.out = nProcess)
    }
    frailtyFunIdx <- is.function(frailty)
    if (! frailtyFunIdx) {
        if (! isNumVector(frailty, error_na = TRUE))
            stop("The frailty effects has to be a numeric vector.")
        frailty <- rep(frailty, length.out = nProcess)
    }

    ## split cases outside the loop for better performance
    z_i <- if (zMatIdx) {function(i) z[i, ]} else {function(i) z}
    origin_i <- if (originFunIdx) {
                    function(i) origin
                } else {
                    function(i) origin[i]
                }
    endTime_i <- if (endTimeFunIdx) {
                     function(i) endTime
                 } else {
                     function(i) endTime[i]
                 }
    frailty_i <- if (frailtyFunIdx) {
                     function(i) frailty
                 } else {
                     function(i) frailty[i]
                 }
    ## generate simulated data for each process
    resList <- lapply(seq_len(nProcess), function(i) {
        simEvent2data(ID = i,
                      simEvent(z = z_i(i),
                               origin = origin_i(i),
                               endTime = endTime_i(i),
                               frailty = frailty_i(i),
                               ...)
                      )
    })

    ## prepare for output
    out <- do.call(rbind, resList)
    if (! attr(out, "recurrent")) {
        uniIdx <- ! duplicated(out$ID)
        out <- out[uniIdx, ]
    }
    ## add original function call to attribute
    attr(out, "call") <- Call
    ## reset row names
    row.names(out) <- NULL
    ## return
    out
}


##' Parametrizations of Covariates and Covariate Coefficients
##'
##' This function helps the parametrizations of covariates and covariate
##' coeffcients when users specify a general hazard rate function in function
##' \code{simEvent} and \code{simEventData}. It applies the specified function
##' (or the built-in option) \code{FUN} to the \eqn{i_{th}} row of the covariate
##' matrix \code{z} and the \eqn{i_{th}} row of the coefficient matrix,
##' iteratively, for \eqn{i} from one to the number of rows of the covariate
##' matrix \code{z}.
##'
##' @usage
##' parametrize(z, zCoef, FUN = c("exponential", "linear", "excess"), ...)
##'
##' @param z A numeric matrix, each row of which represents the covariate vector
##'     at one perticular time point.
##' @param zCoef A numeric matrix, each row of which represents the covariate
##'     coeffcient vector at one perticular time point.
##' @param FUN The parametrization of the model parameter(s) with covariates and
##'     covariate coefficients. The built-in options include
##'     \code{"exponential"}, \code{"linear"}, \code{"excess"} for
##'     parametrization in the exponential, linear, excess relative risk model
##'     form, respectively. It can also be a function that at least has argument
##'     \code{z} and \code{zCoef} for incorporating the covariates and covariate
##'     coefficients into the model. The user-specified function should expect
##'     that both the input \code{z} and \code{zCoef} are numeric vectors and
##'     return a numeric value (or can be convected to a numeric value by
##'     \code{as.numeric}).
##' @param ... Other arguments that can be passed to the function \code{FUN}.
##'
##' @return A numeric vector.
##' @examples
##' ## time points
##' timeVec <- c(0.5, 2)
##' ## time-variant covariates
##' zMat <- cbind(0.5, ifelse(timeVec > 1, 1, 0))
##' ## time-varying coefficients
##' zCoefMat <- cbind(sin(timeVec), timeVec)
##'
##' ## the following three ways are equivalent for the exponential form,
##' ## where the first one (using the built-in option) has the best performance
##' parametrize(zMat, zCoefMat, FUN = "exponential")
##' parametrize(zMat, zCoefMat, function(z, zCoef) exp(z %*% zCoef))
##' sapply(1 : 2, function(i) as.numeric(exp(zMat[i, ] %*% zCoefMat[i, ])))
##'
##' @seealso \code{simEvent}
##' @export
parametrize <- function(z, zCoef,
                        FUN = c("exponential", "linear", "excess"),
                        ...)
{
    if (isCharVector(FUN)) {
        funNames <- c("exponential", "linear", "excess", "none")
        idx <- pmatch(FUN, funNames)[1L]
        if (is.na(idx)) {
            FUN <- match.fun(FUN)
        } else {
            FUN <- switch(funNames[idx],
                          "exponential" = rrisk_exponential,
                          "linear" = rrisk_linear,
                          "excess" = rrisk_excess,
                          "none" = rrisk_none)
            return(FUN(z = z, zCoef))
        }
    }
    if (is.function(FUN)) {
        do.call(.vectorize_rrisk(FUN),
                c(list(z = z, zCoef = zCoef),
                  list(...)["..." %in% names(as.list(args(FUN)))]
                  ))
    } else {
        stop(wrapMessages(
            "Cannot found the function `FUN`."
        ), call. = FALSE)
    }
}


### internal functions =========================================================
## function convert results from simEvent to data.frame
simEvent2data <- function(ID, obj) {
    ## using obj@.Data is much slower
    timeVec <- unclass(obj)
    nTime <- length(timeVec)
    out <- if (nTime > 0L) {
               ## if any event
               data.frame(
                   ID = ID,
                   time = c(timeVec, obj@endTime$endTime),
                   event = c(rep(1, nTime), 0),
                   origin = obj@origin$origin,
                   X = rbind(obj@z$z, obj@censoring$z)
               )
           } else {
               ## if no event
               data.frame(
                   ID = ID,
                   time = obj@endTime$endTime,
                   event = 0,
                   origin = obj@origin$origin,
                   X = obj@censoring$z
               )
           }
    attr(out, "recurrent") <- obj@recurrent
    out
}

## vectorize relative risk function based on input FUN
.vectorize_rrisk <- function(FUN) {
    function(z, zCoef, ...) {
        sapply(seq_len(nrow(z)), function(i) {
            as.numeric(FUN(z = z[i, ], zCoef = zCoef[i, ], ...))
        })
    }
}

## not incorporate z and zCoef via the relative risk function
rrisk_none <- function(z, zCoef, ...) {
    rep(1, nrow(z))
}
