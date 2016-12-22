################################################################################
##
##   R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
##   Copyright (C) 2015-2016
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


##' Generate Simulated Event Times for Each Process (Each Subject)
##' x and beta can be vector of length more than one
##'
##' alpha should be vector of length more than one if df > 1
##' if knots and degree are set, df will be ignore
##' tau can be different for different processes
##'
##' @importFrom splines2 bSpline mSpline
simuData <- function(ID = 1, x = 0, beta = 1, theta = 0.5, alpha = 0.05,
                     df = NULL, knots = NULL, degree = 0L,
                     Boundary.knots = c(0, 168), tau = 168,
                     rho0 = NULL, control = list(), ...) {

    ## quick check:
    ## length(theta) == 1
    ## length(alpha) == length(knots) | length(alpha) == df - degree
    ## length(beta) == length(x)
    ## Boundary.knots == sort(Boundary.knots)
    ## Boundary.knots[2] >= tau

    ## step 1: Calculate the supremum value of rate function for each process
    if (! is.null(rho0) & is.function(rho0)) {
        ## if the baseline rate function, rho0 is specified
        xTime <- seq(from = Boundary.knots[1L], to = Boundary.knots[2L],
                     length.out = 1e3)
        bsRateFun <- rho0(xTime)
        if (any(bsRateFun < 0))
            stop("Baseline rate function must be nonnegative")
    } else {
        xTime <- seq(from = Boundary.knots[1L], to = Boundary.knots[2L],
                     length.out = 1e3)
        bsMat <- splines2::bSpline(x = xTime, df = df, knots = knots,
                                   degree = degree, intercept = TRUE,
                                   Boundary.knots = Boundary.knots)
        bsRateFun <- bsMat %*% alpha
    }
    ## rate functions for a single process
    r <- rgamma(n = 1, shape = theta, scale = 1 / theta)
    tempComp <- as.numeric(r * exp(crossprod(beta, x)))
    rho <- tempComp * bsRateFun
    rho_m <- max(rho)

    ## step 2: Simulate W_i every time
    ## estimate the number of W_i before censoring
    ## iterLim <- ceiling(max(1, tau / qexp(p = 0.05, rate = rho_m)))
    iterLim <- 1e3
    eventTime <- NULL
    lastEventTime <- 0
    for (i in seq_len(iterLim)) {
        W <- rexp(n = iterLim, rate = rho_m)
        ## step 3
        eventTime <- c(eventTime, lastEventTime + cumsum(W))
        lastEventTime <- eventTime[iterLim * i]
        if (lastEventTime > tau)
            break
        if (i == iterLim)
            warning("Fix me by more iterations?")
    }
    eventTime <- eventTime[eventTime < tau]

    ## step 4
    tempn <- length(eventTime)
    U <- runif(n = tempn, min = 0, max = 1)

    if (length(eventTime) == 0L) {
        ## no any event
        rho_t <- 0
    } else if (! is.null(rho0) & is.function(rho0)) {
        rho_t <- tempComp * rho0(eventTime)
    } else {
        bsMat <- splines2::bSpline(x = eventTime, knots = knots,
                                   degree = degree, intercept = TRUE,
                                   Boundary.knots = Boundary.knots)
        rho_t <- tempComp * as.numeric(bsMat %*% alpha)
    }

    ind <- U <= rho_t / rho_m
    timeout <- c(eventTime[ind], tau)
    eventout <- c(rep(1L, length(timeout) - 1L), 0L)
    nRecord <- length(timeout)
    nBeta <- length(beta)
    xNames <- paste0("x", seq_len(nBeta))
    x <- matrix(t(x), ncol = nBeta, nrow = nRecord, byrow = TRUE)
    resMat <- cbind(rep(ID, nRecord), timeout, eventout, x)
    colnames(resMat) <- c("ID", "time", "event", xNames)
    resMat
}


## generate covariates
xFun <- function() c(rbinom(1, 1, 0.5), round(rnorm(1, 0, 1), 2))


## generate dataset for a given number of subjects
## package foreach is used
simuDataset <- function(nSubject, ..., Boundary.knots, x0 = xFun, rho0 = NULL,
                        tau0 = rep(Boundary.knots[2L], nSubject)) {
    tmpList <- vector(mode = "list", length = nSubject)
    for (i in seq_len(nSubject))
        tmpList[[i]] <- simuData(ID = i, ..., Boundary.knots = Boundary.knots,
                                 x = x0(), tau = tau0[i], rho0 = rho0)
    data.frame(do.call(rbind, tmpList))
}
