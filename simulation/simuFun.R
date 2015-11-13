################################################################################
## simulation study on performance of spline rate function
## 200 patients with follow-up period: 24 * 7 = 168 days.
## covariate: treatment group, factor with level: treatment and control
## set five internal knots due to 6 visits
################################################################################

### function part ==============================================================
## generate event times for each process (each subject)
## x and beta can be vector of length more than one
## alpha should be vector of length more than one if df > 1
## if knots and degree are set, df will be ignore
## tau can be different for different processes
simuData <- function (ID = 1, beta = 0.3, theta = 0.5, alpha = 0.06,
                     df = NULL, knots = NULL, degree = 0L, 
                     boundaryKnots = c(0, 168), x = 0, tau = 168) {

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
    ## equivalent to calculate the superemum value of baseline rate function
    if (degree == 0L) { ## if degree == 0L, i.e. piecewise constant
        supBaseRate <- max(bsRateFun <- alpha)
    } else { ## else splines
        xTime <- seq(from = boundaryKnots[1], to = boundaryKnots[2],
                     length.out = 1e3)
        bsMat <- splines::bs(x = xTime, knots = knots, degree = degree,
                             intercept = TRUE, Boundary.knots = boundaryKnots)
        bsRateFun <- bsMat %*% alpha
        supBaseRate <- max(bsRateFun)
    }
    ## rate functions for a single process
    r <- rgamma(n = 1, shape = theta, scale = 1 / theta)
    tempComp <- as.numeric(r * exp(crossprod(beta, x)))
    # browser()
    rho <- tempComp * bsRateFun
    rho_m <- tempComp * supBaseRate 

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
    
    if (degree == 0L) {
        rho_t <- rho[apply(as.array(eventTime), 1, 
                           whereT, bKnots = c(knots, tau))]
    } else if (length(eventTime) == 0L) {
        ## no any event 
        rho_t <- 0 # to avoid argument x in splines::bs being NULL
    } else {
        bsMat <- splines::bs(x = eventTime, knots = knots, degree = degree,
                             intercept = TRUE, Boundary.knots = boundaryKnots)
        rho_t <- tempComp * (bsMat %*% alpha)
    }
    
    ind <- U <= rho_t / rho_m
    timeout <- c(eventTime[ind], tau)
    eventout <- c(rep(1, length(timeout) - 1), 0)
    nRecord <- length(timeout)
    nBeta <- length(beta)
    xNames <- paste("x", seq(nBeta), sep = "")
    x <- matrix(t(x), ncol = nBeta, nrow = nRecord, byrow = TRUE)
    resMat <- cbind(rep(ID, nRecord), timeout, eventout, x)
    attr(resMat, "dimnames")[[2]] <- c("ID", "time", "event", xNames)
    ## return
    resMat
}

## function that extracts estimates and their se
## by default, 6 pieces' piecewise constant rate function
simuFit <- function (nSubject = 200, beta0 = c(0.5, 0.3), theta0 = 0.5,
                     alpha0 = c(0.06, 0.04, 0.05, 0.03, 0.04, 0.05),
                     knots0 = seq(from = 28, to = 140, by = 28), degree0 = 0,
                     boundaryKnots0 = c(0, 168),
                     x0 = function() c(sample(c(0, 1), size = 1),
                                       round(rnorm(1, mean = 0, sd = 1), 2)),
                     tau0 = rep(168, nSubject)) {
    ## generate sample data
    simuDat <- foreach(i = seq(nSubject), .combine = "rbind") %do% {
        simuData(ID = i, beta = beta0, theta = theta0, alpha = alpha0,
                knots = knots0, degree = degree0,
                boundaryKnots = boundaryKnots0, x = x0(), tau = tau0[i])
    }
    simuDat <- data.frame(simuDat)
    simuDat$x1 <- factor(simuDat$x1, levels = c(0, 1),
                         labels = c("Treat", "Contr"))
    colnames(simuDat)[4:5] <- c("group", "x1")
    ## model-fitting
    oneFit <- rateReg(Survr(ID, time, event) ~ group + x1,
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
                       degree = "numeric",
                       df = "numeric",
                       estimates = "list",
                       control = "list",
                       start = "list",
                       na.action = "character",
                       xlevels = "list",
                       contrasts = "list",
                       convergence = "integer", 
                       fisher = "matrix"))
    return(NULL)
}

## summary simulation results
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
