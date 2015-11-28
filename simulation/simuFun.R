################################################################################
## simulation study on performance of spline rate function
## 200 patients with follow-up period: 24 * 7 = 168 days.
## covariate: treatment group, factor with level: treatment and control
## set five internal knots due to 6 visits
################################################################################

### attach packages needed =====================================================
require(splines)
if (! require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if (! require(tidyr)) {install.packages("tidyr"); library(tidyr)}
if (! require(foreach)) {install.packages("foreach"); library(foreach)}
if (! require(doParallel)) {install.packages("doParallel"); library(doParallel)}
if (! require(doRNG)) {install.packages("doRNG"); library(doRNG)}
if (! require(snow)) {install.packages("snow"); library(snow)}
if (! require(rlecuyer)) {install.packages("rlecuyer"); library(rlecuyer)}
if (! require(plyr)) {install.packages("plyr"); library(plyr)}


### function part ==============================================================
## generate event times for each process (each subject)
## x and beta can be vector of length more than one
## alpha should be vector of length more than one if df > 1
## if knots and degree are set, df will be ignore
## tau can be different for different processes
simuData <- function (ID = 1, beta = 0.3, theta = 0.5, alpha = 0.06,
                     df = NULL, knots = NULL, degree = 0L, 
                     boundaryKnots = c(0, 168), x = 0.5, tau = 168,
                     rho0 = NULL, ...) {

    ## quick check:
    ## length(theta) == 1
    ## length(alpha) == length(knots) | length(alpha) == df - degree
    ## length(beta) == length(x)
    ## boundaryKnots == sort(boundaryKnots)
    ## boundaryKnots[2] >= tau
    
    ## internal handy function
    whereT <- function (tt, bKnots) {
        min(which(tt <= bKnots))
    }
    ## set df, knots and degree
    ind <- (is.null(df) + 1) * is.null(knots) + 1
    ## ind == 1: knots is not NULL; df <- length(knots) + 1
    ## ind == 2: df is not NULL, knots is NULL; number of piece <- df
    ## ind == 3: df, knots are both NULL; one-piece constant, df <- 1
    df <- switch(ind, length(knots) + 1, df, 1)
    if (ind > 1) {
        tknots <- df + 1
        knots <- seq.int(from = boundaryKnots[1], to = boundaryKnots[2],
                         length.out = tknots)[-c(1L, tknots)]
    }
    
    ## step 1: Calculate the supremum value of rate function for each process
    if (! is.null(rho0)) { ## if the baseline rate function, rho0 is specified
        xTime <- seq(from = boundaryKnots[1], to = boundaryKnots[2],
                     length.out = 1e3)
        bsRateFun <- rho0(xTime)
        if (any(bsRateFun < 0)) {
            stop("Baseline rate function must be nonnegative")
        }
    } else if (degree == 0L) { ## if degree == 0L, i.e. piecewise constant
        bsRateFun <- alpha
    } else { ## else splines
        xTime <- seq(from = boundaryKnots[1], to = boundaryKnots[2],
                     length.out = 1e3)
        bsMat <- splines::bs(x = xTime, knots = knots, degree = degree,
                             intercept = TRUE, Boundary.knots = boundaryKnots)
        bsRateFun <- bsMat %*% alpha
    }
    ## rate functions for a single process
    r <- rgamma(n = 1, shape = theta, scale = 1 / theta)
    tempComp <- as.numeric(r * exp(crossprod(beta, x)))
    rho <- tempComp * bsRateFun
    rho_m <- max(rho)

    ## step 2: Simulate W_i every time
    ## estimate the number of W_i before censoring
    iterLim <- round(max(1, tau / qexp(p = 0.1, rate = rho_m)))
    nIter <- 1L
    W <- 0
    eventTime <- NULL
    lastEventTime <- 0
    while (lastEventTime < tau) {
        W <- rexp(n = iterLim, rate = rho_m)
        ## step 3
        eventTime <- c(eventTime, lastEventTime + cumsum(W))
        lastEventTime <- tail(eventTime, 1)
        nIter <- nIter + 1L
        if (nIter > 20) stop("Fix me")
    }
    ## evenTime <- unique(round(eventTime, digits = 0)) 
    eventTime <- eventTime[eventTime < tau]

    ## step 4
    tempn <- length(eventTime)
    U <- runif(n = tempn, min = 0, max = 1)
    
    if (length(eventTime) == 0L) {
        ## no any event 
        rho_t <- 0 # to avoid errors
    } else if (! is.null(rho0)) {
        rho_t <- tempComp * rho0(eventTime)
    } else if (degree == 0L) {
        rho_t <- rho[sapply(eventTime, whereT, bKnots = c(knots, tau))]
    } else {
        bsMat <- splines::bs(x = eventTime, knots = knots, degree = degree,
                             intercept = TRUE, Boundary.knots = boundaryKnots)
        rho_t <- tempComp * (bsMat %*% alpha)
    }
    
    ind <- U <= rho_t / rho_m
    timeout <- c(eventTime[ind], tau)
    eventout <- c(rep(1L, length(timeout) - 1), 0L)
    nRecord <- length(timeout)
    nBeta <- length(beta)
    xNames <- paste("x", seq(nBeta), sep = "")
    x <- matrix(t(x), ncol = nBeta, nrow = nRecord, byrow = TRUE)
    resMat <- cbind(rep(ID, nRecord), timeout, eventout, x)
    attr(resMat, "dimnames")[[2]] <- c("ID", "time", "event", xNames)
    ## return
    resMat
}

xFun <- function(seed) {
    if (! missing(seed)) set.seed(seed)
    c(sample(c(0, 1), size = 1),
      round(rnorm(1, mean = 0, sd = 1), 2))
}

simuDataset <- function (nSubject, ..., boundaryKnots0,
                         x0 = xFun, rho0 = NULL,
                         tau0 = rep(boundaryKnots0[2], nSubject)) {
    
    simuDat <- foreach(i = seq(nSubject), .combine = "rbind") %do% { 
        simuData(ID = i, ..., boundaryKnots = boundaryKnots0,
                 x = x0(), tau = tau0[i], rho0 = rho0)
    }
    ## return
    data.frame(simuDat)   
}

## function that extracts estimates and their se
## by default, 6 pieces' piecewise constant rate function
simuEst <- function (nSubject = 200, beta0 = c(0.5, 0.3), theta0 = 0.5,
                     alpha0 = c(0.06, 0.04, 0.05, 0.03, 0.04, 0.05),
                     knots0 = seq(from = 28, to = 140, by = 28), degree0 = 0,
                     boundaryKnots0 = c(0, 168),
                     rho0 = NULL, ...) {
    ## ... can specifiy x0, rho0, and tau0, in function simuDataset
    
    ## generate sample data
    simuDat <- simuDataset(nSubject, beta0, theta0, alpha0, knots0,
                           degree0, boundaryKnots0 = boundaryKnots0,
                           rho0 = rho0, ...)

    ## model-fitting
    oneFit <- rateReg(Survr(ID, time, event) ~ x1 + x2,
                      data = simuDat, knots = knots0, degree = degree0,
                      control = list(boundaryKnots = boundaryKnots0))
    estBeta <- oneFit@estimates$beta[, 1]
    estTheta <- oneFit@estimates$theta[, 1]
    estAlpha <- oneFit@estimates$alpha[, 1]
    seBeta <- oneFit@estimates$beta[, 2]
    seTheta <- oneFit@estimates$theta[, 2]
    seAlpha <- oneFit@estimates$alpha[, 2]
    estVec <- c(estBeta, estTheta, estAlpha)
    seVec <- c(seBeta, seTheta, seAlpha)
    ## return
    resVec <- c(estVec, seVec)
    names(resVec) <- NULL
    resVec
}

## export class 'rateReg' into clusters
exportClass <- function () {
    setClass(Class = "rateReg", 
             slots = c(call = "call", 
                       formula = "formula", 
                       knots = "numeric",
                       boundaryKnots = "numeric",
                       degree = "integer",
                       df = "integer",
                       estimates = "list",
                       control = "list",
                       start = "list",
                       na.action = "character",
                       xlevels = "list",
                       contrasts = "list",
                       convergCode = "integer",
                       logL = "numeric",
                       fisher = "matrix"))

    setClass(Class = "rateRegMcf", 
         slots = c(formula = "formula",
                   knots = "numeric",
                   degree = "integer",
                   boundaryKnots = "numeric",
                   newdata = "matrix",
                   MCF = "data.frame",
                   level = "numeric", 
                   na.action = "character",
                   control = "list", 
                   multiGroup = "logical"))

    return(NULL)
}

## summary simulation results compared with true model
simuSummary <- function (object, beta0 = c(0.5, 0.3), theta0 = 0.5, 
                         alpha0 = c(0.06, 0.04, 0.05, 0.03, 0.04, 0.05)) {
    est0 <- c(beta0, theta0, alpha0)
    nBeta <- length(beta0)
    nAlpha <- length(alpha0)
    nParam <- nBeta + 1L + nAlpha   # add 1 for theta
    if (ncol(object) != 2L * nParam) {
        stop("Number of parameters does not match.")
    }
    barVec <- colMeans(object)
    barEst <- barVec[seq(nParam)]
    barSe <- barVec[seq(nParam + 1, 2 * nParam)]
    seEst <- apply(object[, seq(nParam)], 2, sd)
    res <- cbind(est0, barEst, seEst, barSe)
    row.names(res) <- c(paste("beta", seq(nBeta), sep = ""),
                        "theta", paste("alpha", seq(nAlpha), sep = ""))
    res
}

## sample general rate function
rho0 <- function (t) {
    0.03 * exp(t / 168) + 0.01 * sin(10 * t / 168)
}

## the integral of sample rate function as sample baseline mcf
mu0t <- function (t) {
    0.03 * 168 * (exp(t / 168) - 1) -
        0.01 * (168 / 10) * (cos(10 * t / 168) - 1)
}

## plot fits for general rate funciton 
plotRate <- function (object, df = 9, lenPara = 12, level = 0.95,
                      knots = seq(28, 140, by = 28), degree = 3) {
    xt <- seq(1, 167, length.out = 1e3)
    alphaMat <- object[, seq(lenPara - df + 1, lenPara)]
    bsMat <- bs(xt, knots = knots, degree = degree,
                intercept = TRUE, Boundary.knots = c(0, 168))
    rateEst <- tcrossprod(bsMat, alphaMat)
    meanRate <- bsMat %*% colMeans(alphaMat)
    lowUpp <- matrix(NA, ncol = 2, nrow = 1e3)
    for (j in seq(1e3)) {
        lowUpp[j, ] <- quantile(rateEst[j, ],
                                probs = c((1 - level) / 2,
                                (1 + level) / 2))
    }
    colnames(lowUpp) <- c("lower", "upper")
    ggWideDat <- data.frame("time" = xt, "rho0" = rho0(xt),
                            "mean" = meanRate, lowUpp, rateEst)
    ggLongDat <- gather(ggWideDat, key = type, value = rate, -time)
    ind1 <- ggLongDat$type %in% levels(ggLongDat$type)[1:2]
    indLower <- ggLongDat$type %in% levels(ggLongDat$type)[3]
    indUpper <- ggLongDat$type %in% levels(ggLongDat$type)[4]
    ind2 <- ! (ind1 | indLower | indUpper)
    ## ggOut <- ggplot(data = ggLongDat[ind3, ], aes(x = time)) +
    ##     geom_line(mapping = aes(y = rate, color = type), color = "gray") +
    ggOut <- ggplot(data = ggLongDat, aes(x = time)) + 
        geom_line(data = ggLongDat[ind1, ], aes(y = rate, color = type)) +
        geom_line(data = ggLongDat[indLower, ],
                  aes(y = rate), linetype = "3313") +
        geom_line(data = ggLongDat[indUpper, ],
                  aes(y = rate), linetype = "3313") +
        theme_bw()
    ggOut
}

## function for simulation study testing implementation of delta method for mcf
simuMcf <- function (data, piecesFit, splineFit, ...) {
    ## note that a lot of settings are fixed for convenience
    if (missing(piecesFit)) {
        piecesFit <- rateReg(Survr(ID, time, event) ~ x1 + x2, data = data,
                             knots = seq(28, 140, by = 28))
    }
    if (missing(splineFit)) {
        splineFit <- rateReg(Survr(ID, time, event) ~ x1 + x2, data = data,
                             knots = c(56, 112), degree = 3)      
    }
    ## compute mcf
    mcf1 <- mcf(piecesFit, ...)@MCF
    mcf2 <- mcf(splineFit, ...)@MCF
    ## return
    cbind(mcf1, mcf2)
}

## function helps summerize mcf test results
sumrzMcf <- function (mcfList, newdata = c(1, 0.2),
                      pieces = TRUE, beta0 = c(0, 0), interval = FALSE, ...) {

    mu0t <- function (t, ...) {
        
        0.05 * 168 * (exp(t / 168) - 1) -
            0.02 * (168 / 10) * (cos(10 * t / 168) - 1)
    }

    covEff <- exp(crossprod(newdata, beta0))
    
    mcfList <- if (pieces) {
        lapply(mcfList, function (inpDat) inpDat[, 1:4])
    } else { # else spline
        lapply(mcfList, function (inpDat) inpDat[, 5:8])
    }
    
    mcfDat <- subset(do.call("rbind", mcfList),
                     subset = (! time %in% c(0, 168)))
    if (! interval) {
        seVec <- apply(mcfDat, 1, function (inpVec) {
            (inpVec["upper"] - inpVec["MCF"]) / qnorm(0.975)
        })
        seDat <- cbind(mcfDat[, 1:2], se = seVec)
        resDat <- ddply(seDat, .(time), function (inpDat) {
            meanEstMcf <- mean(inpDat$MCF)
            seEstMcf <- sd(inpDat$MCF)
            meanEstSe <- mean(inpDat$se)
            ## return
            cbind(meanEstMcf, seEstMcf, meanEstSe)
        })    
    } else {
        resDat <- ddply(mcfDat, .(time), function (inpDat) {
            meanEstMcf <- mean(inpDat$MCF)
            intVec <- matrix(quantile(inpDat$MCF,
                                      probs = c(0.025, 0.975)),
                             nrow = 1)
            colnames(intVec) <- c("empirLower", "empirUpper")
            meanLower <- mean(inpDat$lower)
            meanUpper <- mean(inpDat$upper)
            ## return
            cbind(meanEstMcf, intVec, meanLower, meanUpper)
        })
        
    }
    mcf0 <- mu0t(resDat$time) * covEff
    ## return
    cbind(time = resDat[, 1], mcf0, resDat[, -1])
}
    
