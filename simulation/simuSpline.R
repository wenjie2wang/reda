## simulation study on performance of spline rate function
## 200 patients with follow-up period: 24 * 7 = 168 days.
## covariate: treatment group, factor with level: treatment and control
## set five internal knots due to 6 visits

## generate event times for each process (each subject)
## x and beta can be vector of length more than one
## alpha should be vector of length more than one if df > 1
## if knots and degree are set, df will be ignore
## tau can be different for different processes
simuFun <- function (ID = 1, beta = 0.3, theta = 0.5, alpha = 0.06,
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
    ## ind == 2: df is not NULL, while knots is NULL; number of piece <- df
    ## ind == 3: both df and knots are NULL; one-piece constant, df <- 1
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
    nIter <- 1
    W <- 0
    eventTime <- NULL
    lastEventTime <- 0
    while (lastEventTime < tau) {
        W <- rexp(n = iterLim, rate = rho_m)
        ## step 3
        eventTime <- c(eventTime, lastEventTime + cumsum(W))
        lastEventTime <- tail(eventTime, 1)
        nIter <- nIter + 1
        if (nIter == 20) stop("Fix me")
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


## 6 pieces' piecewise constant rate function
## TODO: Generate simulation data to export as example data
set.seed(1216)
nSubject <- 200
simuDat <- foreach(i = seq(nSubject), .combine = "rbind") %do% {
    simuFun(ID = i, beta = c(0.5, 0.3),
            alpha = c(0.06, 0.04, 0.05, 0.03, 0.04, 0.05),
            knots = seq(from = 28, to = 140, by = 28),
            degree = 0, boundaryKnots = c(0, 168),
            x = rbind(ifelse(i <= 100, 0, 1),
                      round(rnorm(1, mean = 0, sd = 1), 2)),
            tau = 168)
}
## simuDat$X2 <- round(simuDat$X2, digits = 2)
## simuDat$X1 <- factor(simuDat$X1, levels = c(0, 1),
##                      labels = c("Treat", "Contr"))
## colnames(simuDat)[4:5] <- c("group", "X1")
## save(simuDat, file = "data/simuDat.RData")


## spline with 2 internal knots and degree 3 and thus df = 6
set.seed(1216)
nSubject <- 200
simuDat <- foreach(i = seq(nSubject), .combine = "rbind") %do% {
    simuFun(ID = i, beta = c(0.5, 0.3),
            alpha = c(0.06, 0.04, 0.05, 0.03, 0.04, 0.05),
            knots = c(56, 112), degree = 3, boundaryKnots = c(0, 168),
            x = rbind(ifelse(i <= 100, 0, 1),
                      round(rnorm(1, mean = 0, sd = 1), 2)),
            tau = 168)
}


