################################################################################
##
##   R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
##   Copyright (C) 2015
##
##   This file is part of the R package reda.
##
##   The R package reda is free software: you can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package reda is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package reda. If not, see <http://www.gnu.org/licenses/>.
##
################################################################################


## collation after class.R
#' @include class.R
NULL


#' Fitting Recurrent Events Regression Model Based on Counts and Rate Function
#'
#' \code{rateReg} returns fitted model results.
#' The default model is a gamma frailty model with a piecewise constant
#' baseline rate function for recurrent event data. 
#'
#' In detail, function \code{rateReg} first constructs the design matrix from
#' the specified arguments: \code{formula}, \code{data}, \code{subset},
#' \code{na.action} and \code{constrasts} before fitting the RECREG model.
#' The constructed design matrix will be checked again to fit in the recurrent
#' event data framework if any observation with missing covariates is removed.
#' (see detail in \code{\link{Survr}} for checking rules).
#' Subject's ID is pinpointed if its observations violate the checking rules.
#' If the ID is not numerical, the appearing order of the subject in the data
#' is pinpointed.
#'
#' The argument \code{start} is an optional list
#' which allows users to specify the initial guess for
#' the parameter values to be estimated.
#' The possible vector elements in the list include
#' \itemize{
#'     \item \code{beta}: coefficient(s) of covariates
#'         set to be 1 by default.
#'     \item \code{theta}: coefficient of random effect
#'         set to be 0.5  by default.
#'     \item \code{alpha}: coefficient(s) of piece-wise constant baseline
#'         rate function set to be 0.15 by default.
#' }
#'
#' The argument \code{control} is an optional list
#' which allows users to control the process of minimization of
#' negative log likelihood function.
#' The possible elements in the list include
#' \itemize{
#'     \item \code{gradtol}: a positive scalar giving the tolerance at
#'         which the scaled gradient is considered close enough to zero
#'         to terminate the algorithm. The default value is 1e-6.
#'     \item \code{stepmax}: a positive scalar which gives the maximum
#'         allowable scaled step length. The default value is 1e5.
#'     \item \code{steptol}: A positive scalar providing the minimum
#'         allowable relative step length. The default value is 1e-6.
#'     \item \code{iterlim}: a positive integer specifying the maximum
#'         number of iterations to be performed before
#'         the program is terminated. The default value is 1e2.
#' }
#' For more details, \code{help(\link[stats]{nlm})}.
#' 
#' @param formula Survr object from function \code{\link{Survr}}. 
#' @param baselinePieces an optional numeric vector consisting of
#' all the right endpoints of baseline pieces.  The default is maximum of time.
#' The default model is of one baseline piece, which  is equivalent to 
#' the negative binomial regression model.
#' @param data an optional data frame, list or environment containing
#' the variables in the model.  If not found in data, the variables are taken 
#' from \code{environment(formula)}, usually the environment from which 
#' function \code{\link{rateReg}} is called.
#' @param subset an optional vector specifying a subset of observations 
#' to be used in the fitting process.
#' @param na.action function which indicates what should the procedure do 
#' if the data contains NAs.  The default is set by the 
#' na.action setting of \code{\link[base]{options}} and is na.fail if that is 
#' not set.  The "factory-fresh" default is \code{\link[stats]{na.omit}}.
#' Another possible value is NULL, no action. 
#' \code{\link[stats]{na.exclude}} can be used as well.
#' \code{help(na.action)} for more details.
#' @param start an optional list of starting values for the parameters
#' to be estimated in the model.
#' @param control an optional list of parameters for controlling the likelihood 
#' function maximization process. For more details, users may \code{help(nlm)}.
#' @param contrasts an optional list, whose entries are values 
#' (numeric matrices or character strings naming functions) to be used 
#' as replacement values for the contrasts replacement function and 
#' whose names are the names of columns of data containing factors.
#' See the \code{contrasts.arg} of 
#' \code{\link[stats]{model.matrix.default}} for more details.
#' @param ... other arguments for future usage.
#' @return a \code{\link{rateReg-class}} object.
#' @references 
#' Fu, Haoda, Junxiang Luo, and Yongming Qu. (2014),
#' "Hypoglycemic Events Analysis via Recurrent Time-to-Event (RECREG) Models," 
#' \emph{Journal of biopharmaceutical statistics}, 2014 Dec 1, Epub 2014 Dec 1.
#' @examples
#' library(reda)
#' ## data(simuDat)
#' regFit <- rateReg(formula = Survr(ID, time, event) ~ x1 + group, 
#'                   data = simuDat, subset = ID %in% 1:100,
#'                   knots = seq(28, 140, length = 5))
#' rateReg(Survr(ID, time, event) ~ x1 + group, 
#'         data = simuDat, subset = ID %in% 75:125)
#' ## str(regFit)
#' show(regFit) # or simply call 'rateRegFit'
#' summary(regFit)
#' coef(regFit)
#' confint(regFit)
#' baseline(regFit)
#' @seealso \code{\link{summary,rateReg-method}}
#' \code{\link{coef,rateReg-method}}
#' \code{\link{confint,rateReg-method}}
#' \code{\link{baseline,rateReg-method}}
#' \code{\link{mcf,rateReg-method}}
#' @importFrom methods new
#' @importFrom stats model.matrix nlm pnorm na.fail na.omit na.exclude na.pass
#' .getXlevels
#' @importFrom splines bs
#' @export
rateReg <- function (formula, df = NULL, knots = NULL, degree = 0L,
                     data, subset, na.action, start = list(), control = list(),
                     contrasts = NULL, ...) {
    ## record the function call to return
    Call <- match.call()

    ## arguments check
    if (missing(formula)) {
        stop("Argument 'formula' is required.")
    } 
    if (missing(data)) {
        data <- environment(formula)
    }
    if (! with(data, inherits(eval(Call[[2]][[2]]), "Survr"))) {
        stop("Response in formula must be a survival recurrent object.")
    }

    ## Prepare data: ID, time, event ~ X(s)
    mcall <- match.call(expand.dots = FALSE)
    mmcall <- match(c("formula", "data", "subset", "na.action"),
                    names(mcall), 0L)
    mcall <- mcall[c(1L, mmcall)]
    ## drop unused levels in factors 
    mcall$drop.unused.levels <- TRUE
    mcall[[1L]] <- quote(stats::model.frame)
    mf <- eval(mcall, parent.frame())
    mt <- attr(mf, "terms")
    mm <- stats::model.matrix(formula, data = mf, contrasts.arg = contrasts)
    ## get data.frame if na.action = na.pass for further data checking 
    mcall$na.action <- na.pass
    mf_na <- eval(mcall, parent.frame())
    mm_na <- stats::model.matrix(formula, data = mf_na,
                                 contrasts.arg = contrasts)
    ## number of covariates excluding intercept
    nBeta <- ncol(mm) - 1 
    ## covariates' names
    covar_names <- colnames(mm)[-1]
    ## data 
    dat <- as.data.frame(cbind(mf[, 1][, 1:3], mm[, -1]))
    colnames(dat) <- c("ID", "time", "event", covar_names)

    ## check the impact caused by missing value
    ## if there is missing value removed
    if (nrow(mm_na) > nrow(dat)) {
        ## recover original ID names for possible pin-point
        idFactor <- with(data, attr(eval(Call[[2]][[2]]), "ID"))
        attr(dat, "ID") <- factor(levels(idFactor)[dat$ID],
                                  levels = levels(idFactor)) 
        message("Observations with missing values on covariates are removed.") 
        message("Checking new data set again ... ", appendLF = FALSE)
        check_Survr(dat)
        message("done.")
    }

    ## 'control' for 'nlm' and 'bs'
    control <- c(control, list(time = dat$time))
    control <- do.call("rateReg_control", control)
    boundaryKnots <- control$boundaryKnots
    indIntercept <- control$intercept

    ## check and reformat 'degree' at the same time
    if ((degree <- as.integer(degree)) < 0) {
        stop("'degree' must be a nonnegative integer.")
    }

    ## generate knots if knots is unspecified
    if (degree == 0L) { ## if piece-wise constant
        templist <- pieceConst(x = dat$time,
                               df = df, knots = knots)
        knots <- templist$knots
        df <- templist$df
    } else { ## else degree > 0, call 'bs' for spline 
        bsMat <- bs(x = dat$time, df = df,
                    knots = knots, degree = degree,
                    intercept = indIntercept, Boundary.knots = boundaryKnots)
        ## update df, knots
        knots <- as.numeric(attr(bsMat, "knots"))
        ## set bKnots as c(knots, last_boundary_knots)
        bKnots <- c(knots, boundaryKnots[2])
        df <- degree + length(knots) + as.integer(indIntercept)
        ## generate bsMat for estimated baseline rate function and mcf
        xTime <- seq(from = min(dat$time), to = max(dat$time),
                     length.out = max(1e3, length(unique(dat$time))))
        xTime <- sort(unique(c(xTime, dat$time[all.equal(dat$event, 0)])))
        bsMat_est <- bs(xTime, knots = knots, degree = degree,
                        intercept = indIntercept,
                        Boundary.knots = boundaryKnots)
    }

    ## set bKnots as c(knots, last_boundary_knots)
    bKnots <- c(knots, boundaryKnots[2])
    alphaName <- nameBases(bKnots = bKnots, degree = degree, df = df, 
                           leftBound = boundaryKnots[1])
                                            
    ## start' values for 'nlm'
    startlist <- c(start, list(nBeta = nBeta, nAlpha = df))
    start <- do.call("rateReg_start", startlist)
    ini <- do.call("c", start)
    length_par <- length(ini)

    ## log likelihood
    fit <- stats::nlm(logL_rateReg, ini, data = dat, 
                      bKnots = bKnots, degree = degree,
                      bsMat = bsMat, bsMat_est = bsMat_est, xTime = xTime,
                      hessian = TRUE,
                      gradtol = control$gradtol, stepmax = control$stepmax,
                      steptol = control$steptol, iterlim = control$iterlim)

    ## estimates for beta
    est_beta <- matrix(NA, nrow = nBeta, ncol = 3)
    colnames(est_beta) <- c("coef", "se", "Pr(>|z|)")
    rownames(est_beta) <- covar_names
    se_vec <- sqrt(diag(solve(fit$hessian)))
    est_beta[, 1] <- fit$estimate[1:nBeta]
    est_beta[, 2] <- se_vec[1:nBeta]
    est_beta[, 3] <- 2 * stats::pnorm(-abs(est_beta[, 1]/est_beta[, 2]))

    ## estimates for theta
    est_theta <- matrix(NA, nrow = 1, ncol = 2)
    colnames(est_theta) <- c("theta", "se")
    est_theta[1, ] <- c(fit$estimate[nBeta + 1], se_vec[nBeta + 1])

    ## estimates for alpha
    est_alpha <- matrix(NA, nrow = df, ncol = 2)
    colnames(est_alpha) <- c("alpha", "se")
    rownames(est_alpha) <- alphaName 
    est_alpha[, 1] <- fit$estimate[(nBeta + 2):length_par]
    est_alpha[, 2] <- se_vec[(nBeta + 2):length_par]

    ## output: na.action
    if (is.null(attr(mf, "na.action"))) {
        na.action <- options("na.action")[[1]]
    } else {
        na.action <- paste("na", class(attr(mf, "na.action")), sep = ".")
    }
    
    ## output: contrasts
    if (is.null(contrasts)) {
        contrasts <- list(contrasts = NA)
    } else {
        contrasts <- attr(mm, "contrasts")    
    }
    
    ## results to return
    results <- methods::new("rateReg", 
                            call = Call, formula = formula, 
                            knots = knots,
                            boundaryKnots = boundaryKnots,
                            degree = degree, df = df,
                            estimates = list(beta = est_beta, 
                                             theta = est_theta, 
                                             alpha = est_alpha),
                            control = control,
                            start = start, 
                            na.action = na.action,
                            xlevels = .getXlevels(mt, mf),
                            contrasts = contrasts,
                            convergCode = fit$code,
                            logL = - fit$minimum,
                            fisher = fit$hessian)
    ## return
    results
}


## internal functions ==========================================================
whereT <- function (tt, bKnots) {
    ## designed to be used inside function 'sapply', 'tt' is length one vector
    ## return the baseline segment number 'tt' belongs to
    ## or the row number of bsMat for censoring time
    min(which(tt <= bKnots))
}

pieceConst <- function (x, df = NULL, knots = NULL) {
    ind <- (is.null(df) + 1) * is.null(knots) + 1
    ## ind == 1: knots is not NULL; df <- length(knots) + 1
    ## ind == 2: df is not NULL, while knots is NULL; number of piece <- df
    ## ind == 3: both df and knots are NULL; one-piece constant, df <- 1
    df <- switch(ind, length(knots) + 1L, as.integer(df), 1L)
    if (ind > 1) {
        tknots <- df + 1L
        quans <- seq.int(from = 0, to = 1,
                         length.out = tknots)[-c(1L, tknots)]
        knots <- as.numeric(stats::quantile(x, quans))
    }
    ## return
    list(df = df, knots = knots)
}

## baseline rate function
rho_0 <- function (par_BaselinePW, Tvec, bKnots, degree, bsMat) {
    ## if piecewise constant, degree == 0
    if (degree == 0) {
        indx <- sapply(Tvec, whereT, bKnots)
        return(par_BaselinePW[indx])  # function ends
    } 
    ## else spline with degree >= 1
    ## return
    bsMat %*% par_BaselinePW
}

## mean cumulative function
mu0 <- function (par_BaselinePW, Tvec, bKnots, degree, bsMat_est, xTime) {
    ## if piecewise constant, degree == 0
    if (degree == 0) {
        ## segement number of each subject
        indx <- sapply(Tvec, whereT, bKnots)
        BL_segments <- c(bKnots[1], diff(bKnots))
        ## The MCF at each time point  
        CumMean_Pieces <- diffinv(BL_segments * par_BaselinePW)[-1]  
        mu_tau <- CumMean_Pieces[indx] -
            (bKnots[indx] - Tvec) * par_BaselinePW[indx]
        return(mu_tau)  # function ends
    }
    ## else spline with degree > 1
    stepTime <- xTime[2] - xTime[1]
    baseRate <- bsMat_est %*% par_BaselinePW
    indx <- sapply(Tvec, whereT, bKnots = xTime)
    mu_tau <- sapply(indx, function (ind) {
        sum(baseRate[seq(ind)]) * stepTime
    })
    ## return
    mu_tau
}

dmu0_dalpha <- function (tt, bKnots, degree, bsMat_est, xTime) {
    ## if baseline rate function is piecewise constant
    if (degree == 0L) {
        indx <- min(which(tt <= bKnots))
        ## BL_segments 
        value <- diff(c(0, bKnots))
        n_pieces <- length(bKnots)
        ## if tt lies in the last segment
        if (indx == n_pieces) {
            value[n_pieces] <- ifelse(n_pieces == 1, tt, 
                                      tt - bKnots[n_pieces - 1])
        } else if (indx > 1) { ## if tt lies in one of the middle segments
            value[(indx + 1) : n_pieces] <- 0
            value[indx] <- tt - bKnots[indx - 1]
        } else { ## if tt lies in the first segment
            value[(indx + 1) : n_pieces] <- 0
            value[indx] <- tt
        }
        ## return and end the function
        return(value)    
    }
    ## else it is spline with degree > 1
    stepTime <- xTime[2] - xTime[1]
    indx <- min(which(tt <= xTime))
    derVec <- sapply(indx, function (ind) {
        colSums(bsMat_est[seq(ind), ]) * stepTime
    })
    ## return
    derVec
}

dl_dalpha_part1 <- function (par_alpha, indx, degree, bsMat) {
    ## if rate function is piecewise constant
    if (degree == 0L) {
        return(1 / par_alpha * table(indx))     
    }
    ## else rate function is spline
    drho0_dbeta_ij <- bsMat
    rho0_ij <- bsMat %*% par_alpha
    ## return
    colSums(t(1 / rho0_ij) %*% drho0_dbeta_ij)
}

## compute negative log likelihood
logL_rateReg <- function (par, data, bKnots, degree, bsMat, bsMat_est, xTime) {
    nBeta <- ncol(data) - 3
    ## par = \THETA in the paper
    par_beta <- par[1 : nBeta]
    par_theta <- par[nBeta + 1]
    par_alpha <- par[(nBeta + 2) : length(par)]
    m <- length(unique(data$ID))
    expXBeta <- exp(as.matrix(data[, 4 : (3 + nBeta)]) %*% as.matrix(par_beta))
    ## index for event and censoring
    ind_event <- data$event == 1
    ind_cens <- data$event == 0
    ## rate function
    rho_0_ij <- rho_0(par_BaselinePW = par_alpha,
                      Tvec = data$time[ind_event],
                      bKnots = bKnots, degree = degree, 
                      bsMat = bsMat[ind_event, ])
    rho_i <- expXBeta[ind_event] * rho_0_ij
    rho_i[rho_i < 1e-100] <- 1e-100
    sum_log_rho_i <- sum(log(rho_i))
    ## n_ij: number of event for each subject
    ## these codes to make sure that the order will not change 
    ## if the patient ID is not ordered
    n_ij <- table(data$ID)[order(unique(data$ID))] - 1  
    ## if there is a subject with 0 event, 
    ## the sequence will not be generated for this subject
    theta_j_1 <- par_theta + sequence(n_ij) - 1  
    theta_j_1[theta_j_1 < 1e-100] <- 1e-100
    sum_log_theta_j_1 <- sum(log(theta_j_1))
    ## integral that involves censoring time tau
    ## baseline mcf
    mu0i <- mu0(par_BaselinePW = par_alpha, Tvec = data$time[ind_cens],
                bKnots = bKnots, degree = degree, bsMat_est = bsMat_est,
                xTime = xTime)
    mui <- mu0i * expXBeta[ind_cens]
    mui_theta <- par_theta + mui
    mui_theta[mui_theta < 1e-100] <- 1e-100
    sum_log_theta_mui <- sum((n_ij + par_theta) * log(mui_theta))
    if (par_theta < 1e-100) {
        par_theta <- 1e-100
    }
    logLH <- m * par_theta * log(par_theta) + sum_log_rho_i + 
        sum_log_theta_j_1 - sum_log_theta_mui
    penal <- ifelse(par_theta < 0 | min(par_alpha) < 0, 1e+50, 0)
    negLH <- -logLH + penal
    ## Calculate the gradient
    xMat_i <- as.matrix(data[ind_cens, 4:(3 + nBeta)])
    dl_dbeta_i <- sweep(x = xMat_i, MARGIN = 1, FUN = "*", 
                        STATS = (n_ij - mui)/(par_theta + mui) * par_theta)
    dl_dbeta <- colSums(dl_dbeta_i)
    dl_dtheta <- m + m * log(par_theta) + 
        sum(1/(par_theta + sequence(n_ij) - 1)) - 
        sum((n_ij + par_theta)/(par_theta + mui)) - sum(log(mui_theta))
    indx <- sapply(data$time[ind_event], whereT, bKnots)
    if (length(unique(indx)) < length(bKnots)) {
        stop("Some segements have zero events!")
    }
    ## reform dimension by 'array' for one-piece baseline 
    dim_n1 <- length(par_alpha)
    dim_n2 <- length(data$time[ind_cens])
    tempPart2 <- array(sapply(data$time[ind_cens], dmu0_dalpha,
                              bKnots, degree, bsMat_est, xTime),
                       c(dim_n1, dim_n2))
    dl_dalpha_part2 <- sweep(t(tempPart2), MARGIN = 1, FUN = "*", 
                             STATS = (n_ij + par_theta) / (par_theta + mui) *
                                 expXBeta[ind_cens])
    tempPart2 <- colSums(dl_dalpha_part2)
    tempPart1 <- dl_dalpha_part1(par_alpha, indx, degree, bsMat[ind_event, ])
    dl_dalpha <-  tempPart1 - tempPart2
    attr(negLH, "gradient") <- -c(dl_dbeta, dl_dtheta, dl_dalpha)
    ## return
    negLH
}

rateReg_control <- function (gradtol = 1e-6, stepmax = 1e5, 
                           steptol = 1e-6, iterlim = 1e2,
                           boundaryKnots = NULL, intercept = TRUE,
                           time) {
    ## controls for function stats::nlm
    if (!is.numeric(gradtol) || gradtol <= 0) {
        stop("value of 'gradtol' must be > 0")
    }
    if (!is.numeric(stepmax) || stepmax <= 0) {
        stop("value of 'stepmax' must be > 0")
    } 
    if (!is.numeric(steptol) || steptol <= 0) {
        stop("value of 'steptol' must be > 0")
    } 
    if (!is.numeric(iterlim) || iterlim <= 0) {
        stop("maximum number of iterations must be > 0")
    }
    if (is.null(boundaryKnots)) {
        boundaryKnots <- c(0, max(time))
    } else {
        boundaryKnots <- sort(boundaryKnots)
        ind1 <- boundaryKnots[1] > min(time)
        ind2 <- boundaryKnots[2] < max(time)
        if (ind1 || ind2) {
            stop("boundary knots should not lie inside the range of visit time")
        }
    }
    ## return
    list(gradtol = gradtol, stepmax = stepmax, 
         steptol = steptol, iterlim = iterlim,
         boundaryKnots = boundaryKnots, intercept = intercept)
}

rateReg_start <- function (beta, theta = 0.5, alpha, nBeta, nAlpha) {
    ## beta = starting value(s) for coefficients of covariates
    ## theta = starting value for random effects
    ## alpha = starting values for piece-wise baseline rate functions
    if (missing(beta)) {
        beta <- rep(1, nBeta)
    } else if (length(beta) != nBeta) {
        stop(paste("number of starting values for coefficients of covariates",
                   "does not match with the specified formula"))
    }
    if (theta <= 0) {
        stop("value of parameter for random effects must be > 0")
    }
    if (missing(alpha)) {
        alpha <- rep(0.15, nAlpha)
    }
    ## return
    list(beta = beta, theta = theta, alpha = alpha)
}

## generate intervals from specified baseline pieces
nameBases <- function (bKnots, degree, df, leftBound) {
    nAlpha <- length(bKnots)
    intervals <- rep(NA, nAlpha)
    if (degree == 0L) {
        intervals[1] <- paste0("(", leftBound, ", ", bKnots[1], "]", sep = "")
        if (nAlpha > 1) {
            for(i in 2:nAlpha){
                intervals[i] <- paste0("(", bKnots[i - 1], ", ", 
                                       bKnots[i], "]", sep = "")
            }
        }
        return(intervals)
    }
    ## else degree > 0 
    ## return
    paste("B-spline", seq(df))
}
