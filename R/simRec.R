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


##' Generate Simulated Recurrent Events for One Process
##'
##' The thinning method (Lewis and Shedler, 1979) is applied for generating
##' simulated recurrent event times from Gamma frailty model for one process or
##' usually one subject or machine in applications. Conditional on predictors
##' (or covariates) and the unobserved Gamma frailty effect, the process is a
##' Poisson process.
##'
##' x and beta can be vector of length more than one.
##' alpha should be vector of length more than one if df > 1
##' if knots and degree are set, df will be ignore
##' tau can be different for different processes
##'
##' @references
##' Lewis, P. A., & G. S. Shedler. (1979). Simulation of
##' Nonhomogeneous Poisson Processes by Thinning.
##' \emph{Naval Research Logistics Quarterly},
##' 26(3), Wiley Online Library: 403--13.
##'
##' @importFrom stats optimize rexp qexp
simRec <- function(x = 0, x.coef = 1,
                   origin = 0, endTime = 5,
                   rho.fun = function(x) 1,
                   rho.args = list(),
                   rho.coef = 1,
                   frailty.fun = rgamma,
                   frailty.args = list(shape = 0.5, scale = 2),
                   ...)
{
    ## TODO: add checks on inputs
    ## if (any(bsRateFun < 0))
    ##         stop("Baseline rate function must be nonnegative")
    ## }

    ## covariate effect
    expXbeta <- as.numeric(exp(x %*% x.coef))

    ## prepare frailty effect
    frailty.args <- c(list(n = 1), frailty.args[names(frailty.args) != "n"])
    frailty <- do.call(frailty.fun, frailty.args)

    ## rate function
    rateFun <- function(xVec) {
        rhoMat <- do.call(rho.fun, c(list(xVec), rho.args))
        rhoVec <- as.numeric(rhoMat %*% rho.coef)
        frailty * rhoVec * expXbeta
    }

    ## step 1: calculate the supremum value of rate function
    rhoMaxObj <- stats::optimize(rateFun,
                                 lower = origin,
                                 upper = endTime,
                                 maximum = TRUE)
    rho_max <- rhoMaxObj$objective

    ## step 2: generate W_i in batch for possible better performance
    ## estimate the number of W_i before censoring
    batchNum <- ceiling((endTime - origin) / stats::qexp(0.10, rho_max))
    eventTime <- NULL
    lastEventTime <- origin
    while (lastEventTime < endTime) {
        W <- stats::rexp(n = batchNum, rate = rho_max)
        ## step 3: update evnet times
        eventTime <- c(eventTime, lastEventTime + cumsum(W))
        lastEventTime <- eventTime[length(eventTime)]
    }
    eventTime <- eventTime[eventTime <= endTime]

    ## prepare output in a list
    outList <- list(x = x, x.coef = x.coef,
                    origin = origin, endTime = endTime)
    if (! length(eventTime)) {
        ## no any event
        out <- numeric(0)
    } else {
        ## step 4: thinning
        rho_t <- rateFun(eventTime)
        U <- runif(n = length(eventTime))
        ind <- U <= rho_t / rho_max
        out <- eventTime[ind]
    }
    ## return
    attr(out, "settings") <- outList
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
