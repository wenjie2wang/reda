##
## R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
## Copyright (C) 2015-2020
##
## This file is part of the R package reda.
##
## The R package reda is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package reda is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##


## collation after class.R
##' @include class.R
NULL


##' Recurrent Events Regression Based on Counts and Rate Function
##'
##' This function fits recurrent event data (event counts) by gamma frailty
##' model with spline rate function. The default model is the gamma frailty
##' model with one piece constant baseline rate function, which is equivalent to
##' negative binomial regression with the same shape and rate parameter in the
##' gamma prior. Spline (including piecewise constant) baseline hazard rate
##' function can be specified for the model fitting.
##'
##' Function \code{\link{Recur}} in the formula response by default first checks
##' the dataset and will report an error if the dataset does not fall into
##' recurrent event data framework.  Subject's ID will be pinpointed if its
##' observation violates any checking rule. See \code{\link{Recur}} for all the
##' checking rules.
##'
##' Function \code{rateReg} first constructs the design matrix from
##' the specified arguments: \code{formula}, \code{data}, \code{subset},
##' \code{na.action} and \code{constrasts} before model fitting.
##' The constructed design matrix will be checked again to
##' fit the recurrent event data framework
##' if any observation with missing covariates is removed.
##'
##' The model fitting process involves minimization of negative log
##' likelihood function, which calls function \code{\link[stats]{constrOptim}}
##' internally. \code{help(constrOptim)} for more details.
##'
##' The argument \code{start} is an optional list
##' that allows users to specify the initial guess for
##' the parameter values for the minimization of
##' negative log likelihood function.
##' The available numeric vector elements in the list include
##' \itemize{
##'     \item \code{beta}: Coefficient(s) of covariates,
##'         set to be all 0.1 by default.
##'     \item \code{theta}: Parameter in Gamma(theta, 1 / theta) for
##'         frailty random effect, set to be 0.5 by default.
##'     \item \code{alpha}: Coefficient(s) of baseline rate function,
##'         set to be all 0.05 by default.
##' }
##' The argument \code{control} is an optional list
##' that allows users to control the process of minimization of
##' negative log likelihood function passed to \code{constrOptim}
##' and specify the boundary knots of baseline rate function.
##' The available options additional to those that can be passed from
##' \code{control} to \code{constrOptim} include
##' \itemize{
##'     \item \code{Boundary.knots}: A length-two numeric vector to specify
##'         the boundary knots for baseline rate funtion. By default,
##'         the left boundary knot is the smallest origin time and
##'         the right one takes the largest censoring time from data.
##'     \item \code{verbose}: A optional logical value with default \code{TRUE}.
##'         Set it to be \code{FALSE} to supress any possible message
##'         from this function.
##' }
##'
##' @param formula \code{Recur} object produced by function \code{\link{Recur}}.
##'     The terminal events and risk-free episodes specified in \code{Recur}
##'     will be ignored since the model does not support them.
##' @param data An optional data frame, list or environment containing the
##'     variables in the model.  If not found in data, the variables are taken
##'     from \code{environment(formula)}, usually the environment from which
##'     function \code{\link{rateReg}} is called.
##' @param subset An optional vector specifying a subset of observations to be
##'     used in the fitting process.
##' @param df An optional nonnegative integer to specify the degree of freedom
##'     of baseline rate function. If argument \code{knots} or \code{degree} are
##'     specified, \code{df} will be neglected whether it is specified or not.
##' @param knots An optional numeric vector that represents all the internal
##'     knots of baseline rate function.  The default is \code{NULL},
##'     representing no any internal knots.
##' @param degree An optional nonnegative integer to specify the degree of
##'     spline bases.
##' @param na.action A function that indicates what should the procedure do if
##'     the data contains \code{NA}s.  The default is set by the na.action
##'     setting of \code{options}.  The "factory-fresh" default is
##'     \code{na.omit}.  Other possible values inlcude
##'     \code{na.fail}, \code{na.exclude}, and
##'     \code{na.pass}.  \code{help(na.fail)} for details.
##' @param spline An optional character that specifies the flavor of splines.
##'     The possible option is \code{bSplines} for B-splines or \code{mSplines}
##'     for M-splines.
##' @param start An optional list of starting values for the parameters to be
##'     estimated in the model.  See more in Section details.
##' @param control An optional list of parameters to control the maximization
##'     process of negative log likelihood function and adjust the baseline rate
##'     function.  See more in Section details.
##' @param contrasts An optional list, whose entries are values (numeric
##'     matrices or character strings naming functions) to be used as
##'     replacement values for the contrasts replacement function and whose
##'     names are the names of columns of data containing factors.  See
##'     \code{contrasts.arg} of \code{model.matrix.default} for details.
##' @param ... Other arguments for future usage.
##'
##' @return A \code{rateReg} object, whose slots include
##' \itemize{
##'     \item \code{call}: Function call of \code{rateReg}.
##'     \item \code{formula}: Formula used in the model fitting.
##'     \item \code{nObs}: Number of observations.
##'     \item \code{spline}: A list contains
##'         \itemize{
##'             \item \code{spline}: The name of splines used.
##'             \item \code{knots}: Internal knots specified for the baseline
##'                 rate function.
##'             \item \code{Boundary.knots}: Boundary knots specified for the
##'                 baseline rate function.
##'             \item \code{degree}: Degree of spline bases specified in
##'                 baseline rate function.
##'             \item \code{df}: Degree of freedom of the model specified.
##'     }
##'     \item \code{estimates}: Estimated coefficients of covariates and
##'         baseline rate function, and estimated rate parameter of
##'         gamma frailty variable.
##'     \item \code{control}: The control list specified for model fitting.
##'     \item \code{start}: The initial guess specified for the parameters
##'         to be estimated.
##'     \item \code{na.action}: The procedure specified to deal with
##'         missing values in the covariate.
##'     \item \code{xlevels}: A list that records the levels in
##'         each factor variable.
##'     \item \code{contrasts}: Contrasts specified and used for each
##'         factor variable.
##'     \item \code{convergCode}: \code{code} returned by function
##'         \code{\link[stats]{optim}}, which is an integer indicating why the
##'         optimization process terminated. \code{help(optim)} for details.
##'     \item \code{logL}: Log likelihood of the fitted model.
##'     \item \code{fisher}: Observed Fisher information matrix.
##' }
##'
##' @references
##'
##' Fu, H., Luo, J., & Qu, Y. (2016).
##' Hypoglycemic events analysis via recurrent time-to-event (HEART) models.
##' \emph{Journal Of Biopharmaceutical Statistics}, 26(2), 280--298.
##'
##' @example inst/examples/ex_rateReg.R
##'
##' @seealso
##' \code{\link{summary,rateReg-method}} for summary of fitted model;
##' \code{\link{coef,rateReg-method}} for estimated covariate coefficients;
##' \code{\link{confint,rateReg-method}} for confidence interval of
##' covariate coefficients;
##' \code{\link{baseRate,rateReg-method}} for estimated coefficients of baseline
##' rate function;
##' \code{\link{mcf,rateReg-method}} for estimated MCF from a fitted model;
##' \code{\link{plot,mcf.rateReg-method}} for plotting estimated MCF.
##'
##' @importFrom splines2 bSpline cSpline ibs iSpline mSpline
##'
##' @importFrom stats .getXlevels constrOptim model.extract na.fail na.omit
##'     na.exclude na.pass predict
##'
##' @export
rateReg <- function(formula, data, subset, df = NULL, knots = NULL, degree = 0L,
                    na.action, spline = c("bSplines", "mSplines"),
                    start = list(), control = list(), contrasts = NULL, ...) {
    ## record the function call to return
    Call <- match.call()
    spline <- match.arg(spline)
    if (missing(formula))
        stop("Argument 'formula' is required.")
    if (missing(data))
        data <- environment(formula)

    ## take care of subset individual for possible non-numeric ID
    if (! missing(subset)) {
        sSubset <- substitute(subset)
        subIdx <- eval(sSubset, data, parent.frame())
        if (! is.logical(subIdx))
            stop("'subset' must be logical")
        subIdx <- subIdx & ! is.na(subIdx)
        data <- data[subIdx, ]
    }

    ## Prepare data: ID, time, event ~ X(s)
    mcall <- match.call(expand.dots = FALSE)
    mmcall <- match(c("formula", "data", "na.action"), names(mcall), 0L)
    mcall <- mcall[c(1L, mmcall)]
    ## re-define data
    ## mcall$data <- substitute(data)
    mcall$data <- data
    ## drop unused levels in factors
    mcall$drop.unused.levels <- TRUE
    mcall[[1L]] <- quote(stats::model.frame)
    mf <- eval(mcall, parent.frame())
    mt <- attr(mf, "terms")
    mm <- stats::model.matrix(formula, data = mf, contrasts.arg = contrasts)

    ## check response constructed from Recur
    resp <- stats::model.extract(mf, "response")
    if (! (is.Recur(resp) || is.Survr(resp)))
        stop("Response in the formula must be an 'Recur' object.")

    ## number of covariates excluding intercept
    nBeta <- ncol(mm) - 1L
    ## if (nBeta == 0L)
    ##     warning("No covariate was specified in formula.")
    ## covariates' names
    covar_names <- colnames(mm)[- 1L]

    ## 'control' for optimization and splines' boundary knots
    control <- do.call("rateReg_control", control)
    control4rateReg <- control$control4rateReg
    control4optim <- control$control4optim

    ## for possible missing values in covaraites
    if (length(na.action <- attr(mf, "na.action"))) {
        ## update if there is missing value removed
        attr(resp, "ID") <- attr(resp, "ID")[- na.action]
        ## check data for possible error caused by removal of missing values
        if (control4rateReg$verbose)
            message("Observations with missing value in covariates ",
                    "are removed.\nChecking the new dataset again...\n",
                    appendLF = FALSE)
        if (is.Recur(resp)) {
            resp <- check_Recur(resp, check = "hard")
        } else {
            resp <- check_Survr(resp, check = TRUE)
        }
        if (control4rateReg$verbose)
            message("Done!")
    }

    ## sorted data by ID, time, and event
    ord <- attr(resp, "ord")
    ## data matrix processed
    xMat <- mm[ord, - 1L, drop = FALSE]
    dat <- if (is.Recur(resp)) {
               as.data.frame(cbind(
                   resp[ord, c("id", "time2", "event", "origin")], xMat
               ))
           } else {
               as.data.frame(cbind(
                   resp[ord, ], xMat
               ))
           }

    colnames(dat) <- c("ID", "time", "event", "origin", covar_names)
    nObs <- nrow(dat)

    ## set up boundary knots
    Boundary.knots <- if (is.null(control4rateReg$Boundary.knots)) {
                          with(dat, c(min(origin, na.rm = TRUE),
                                      max(time, na.rm = TRUE)))
                      } else {
                          control4rateReg$Boundary.knots
                      }

    ## generate knots if knots is unspecified
    if (spline == "bSplines") {
        ## B-spline
        iMat <- splines2::ibs(x = dat$time, df = df, knots = knots,
                              degree = degree, intercept = TRUE,
                              Boundary.knots = Boundary.knots)
        bMat <- attr(iMat, "bsMat")
    } else {
        ## M-spline
        iMat <- splines2::iSpline(x = dat$time, df = df, knots = knots,
                                  degree = degree, intercept = TRUE,
                                  Boundary.knots = Boundary.knots)
        bMat <- attr(iMat, "msMat")
    }
    iMat0 <- predict(iMat, dat$origin)
    iMat <- iMat - iMat0

    ## update df, knots, degree, and Boundary.knots
    knots <- as.numeric(attr(iMat, "knots"))
    df <- (degree <- attr(iMat, "degree")) + length(knots) + 1L
    Boundary.knots <- attr(iMat, "Boundary.knots")
    ## name each basis for alpha output
    alphaName <- nameBases(df = df, spline = spline)

    ## start' values
    startlist <- c(start, list(nBeta_ = nBeta, nAlpha_ = df))
    start <- do.call("rateReg_start", startlist)
    ini <- do.call("c", start)
    length_par <- length(ini)

    ## check whether the knots are reasonable
    if (any(colSums(bMat[dat$event > 0, , drop = FALSE]) == 0)) {
        stop(wrapMessages(
            "Some spline basis does not capture any event time",
            "and thus is possibly redundent.",
            "Please adjust spline knots or degree."
        ))
    }

    ## prepare anything needed in LogL_rateReg but free from parameters
    ## index for event and censoring
    ind_event <- dat$event > 0
    ind_cens <- ! ind_event
    ## basis matrix at event times
    bMat_event <- bMat[ind_event, , drop = FALSE]
    ## n_ij: number of event for each subject
    ## the following code makes sure the order will not change
    ## if the patient ID is not ordered
    ## not needed if data are already sorted by ID
    ## n_ij <- table(dat$ID)[order(unique(dat$ID))] - 1L
    n_ij <- table(dat$ID) - 1L
    seq_n_ij <- sequence(n_ij)
    nSub <- length(n_ij)
    dmu0_dalpha <- iMat[ind_cens, , drop = FALSE]
    xMat_i <- xMat[ind_cens, , drop = FALSE]

    ## log likelihood
    ## theta > 0 and alph >= 0
    fit <- stats::constrOptim(ini, f = logL_rateReg, grad = logL_rateReg_grad,
                              ui = cbind(matrix(0, length_par - nBeta, nBeta),
                                         diag(length_par - nBeta)),
                              ci = c(sqrt(.Machine$double.eps),
                                     rep(0, length_par - nBeta - 1)),
                              hessian = TRUE,
                              nBeta = nBeta,
                              nSub = nSub,
                              xMat = xMat,
                              ind_event = ind_event,
                              ind_cens = ind_cens,
                              bMat_event = bMat_event,
                              n_ij = n_ij,
                              seq_n_ij = seq_n_ij,
                              dmu0_dalpha = dmu0_dalpha,
                              xMat_i = xMat_i,
                              control = control4optim)

    ## estimates for beta
    est_beta <- matrix(NA_real_, nrow = nBeta, ncol = 5L)
    colnames(est_beta) <- c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")
    rownames(est_beta) <- covar_names

    se_vec <- sqrt(diag(solve(fit$hessian)))
    est_beta[, 1L] <- fit$par[seq_len(nBeta)]
    est_beta[, 2L] <- exp(est_beta[, 1L])
    est_beta[, 3L] <- se_vec[seq_len(nBeta)]
    est_beta[, 4L] <- est_beta[, 1L] / est_beta[, 3L]
    est_beta[, 5L] <- 2 * stats::pnorm(- abs(est_beta[, 4L]))

    ## estimates for theta
    est_theta <- matrix(NA_real_, nrow = 1L, ncol = 2L)
    colnames(est_theta) <- c("parameter", "se")
    rownames(est_theta) <- "Frailty"
    est_theta[1L, ] <- c(fit$par[nBeta + 1L], se_vec[nBeta + 1L])

    ## estimates for alpha
    est_alpha <- matrix(NA_real_, nrow = df, ncol = 2)
    colnames(est_alpha) <- c("coef", "se(coef)")
    rownames(est_alpha) <- alphaName
    est_alpha[, 1L] <- fit$par[(tmpIdx <- (nBeta + 2L) : length_par)]
    est_alpha[, 2L] <- se_vec[tmpIdx]

    ## output: na.action
    na.action <- if (is.null(na.action)) {
                     options("na.action")[[1L]]
                 } else {
                     paste0("na.", class(attr(mf, "na.action")))
                 }

    ## output: contrasts
    contrasts <- if (is.null(contrasts)) {
                     list(contrasts = NULL)
                 } else {
                     attr(mm, "contrasts")
                 }

    ## output: df, degree of freefom, including beta and theta
    df <- list(beta = nBeta, theta = 1L, alpha = df)

    ## output: xlevels
    xlevels <- .getXlevels(mt, mf)
    if (is.null(xlevels)) xlevels <- list()

    ## return
    methods::new("rateReg",
                 call = Call,
                 formula = formula,
                 nObs = nObs,
                 spline = list(spline = spline,
                               df = df,
                               knots = knots,
                               degree = degree,
                               Boundary.knots = Boundary.knots),
                 estimates = list(beta = est_beta,
                                  theta = est_theta,
                                  alpha = est_alpha),
                 control = control,
                 start = start,
                 na.action = na.action,
                 xlevels = xlevels,
                 contrasts = contrasts,
                 convergCode = fit$convergence,
                 logL = - fit$value,
                 fisher = fit$hessian)
}


### internal functions =========================================================
## compute negative log likelihood
logL_rateReg <- function(par, nBeta, nSub, xMat, ind_event, ind_cens,
                         bMat_event, n_ij, seq_n_ij, dmu0_dalpha, xMat_i)
{
    ## par = \THETA in the paper
    par_theta <- max(par[nBeta + 1L], .Machine$double.eps)
    par_alpha <- par[(nBeta + 2L) : length(par)]
    expXBeta <-
        if (nBeta) {
            as.numeric(exp(xMat %*% as.matrix(par[seq_len(nBeta)])))
        } else {
            rep(1, nrow(xMat))
        }

    ## baseline rate function
    rho_0_ij <- pmax(as.numeric(bMat_event %*% par_alpha), .Machine$double.eps)
    rho_i <- pmax(expXBeta[ind_event] * rho_0_ij, .Machine$double.eps)
    sum_log_rho_i <- sum(log(rho_i))

    theta_j_1 <- par_theta + seq_n_ij - 1
    sum_log_theta_j_1 <- sum(log(theta_j_1))

    ## baseline mcf, integral of rho_0 that involves censoring time tau
    mu0i <- as.numeric(dmu0_dalpha %*% par_alpha)
    mui <- mu0i * expXBeta[ind_cens]
    mui_theta <- pmax(par_theta + mui, .Machine$double.eps)
    sum_log_theta_mui <- sum((n_ij_theta <- n_ij + par_theta) * log(mui_theta))

    ## log likelihood function
    logLH <- nSub * par_theta * log(par_theta) + sum_log_rho_i +
        sum_log_theta_j_1 - sum_log_theta_mui

    ## return negative log likelihood
    negLH <- - logLH

    ## ## Calculate the gradient
    ## ## on beta, vector
    ## dl_dbeta_i <- sweep(x = as.matrix(xMat_i), MARGIN = 1, FUN = "*",
    ##                     STATS = par_theta * (n_ij - mui) / (par_theta + mui))
    ## dl_dbeta <- colSums(dl_dbeta_i)

    ## ## on theta
    ## dl_dtheta <- nSub + nSub * log(par_theta) +
    ##     sum(1 / (par_theta + seq_n_ij - 1)) -
    ##     sum((n_ij + par_theta) / (par_theta + mui)) - sum(log(mui_theta))

    ## ## on alpha
    ## part1 <- crossprod(1 / rho_0_ij, bMat_event)
    ## dl_dalpha_part2 <- sweep(dmu0_dalpha, MARGIN = 1, FUN = "*",
    ##                          STATS = expXBeta[ind_cens] * n_ij_theta /
    ##                              mui_theta)
    ## part2 <- colSums(dl_dalpha_part2)
    ## dl_dalpha <-  part1 - part2
    ## ## return gradient as one attribute
    ## attr(negLH, "gradient") <- - c(dl_dbeta, dl_dtheta, dl_dalpha)
    negLH
}


## compute negative log likelihood
logL_rateReg_grad <- function(par, nBeta, nSub, xMat, ind_event, ind_cens,
                              bMat_event, n_ij, seq_n_ij, dmu0_dalpha, xMat_i)
{
    ## par = \THETA in the paper
    par_theta <- max(par[nBeta + 1L], .Machine$double.eps)
    n_ij_theta <- n_ij + par_theta
    par_alpha <- par[(nBeta + 2L) : length(par)]
    expXBeta <- as.numeric(exp(xMat %*% as.matrix(par[seq_len(nBeta)])))

    ## baseline rate function
    rho_0_ij <- pmax(as.numeric(bMat_event %*% par_alpha), .Machine$double.eps)

    ## baseline mcf, integral of rho_0 that involves censoring time tau
    mu0i <- as.numeric(dmu0_dalpha %*% par_alpha)
    mui <- mu0i * expXBeta[ind_cens]
    mui_theta <- pmax(par_theta + mui, .Machine$double.eps)

    ## Calculate the gradient
    ## on beta, vector
    dl_dbeta_i <- sweep(x = as.matrix(xMat_i), MARGIN = 1, FUN = "*",
                        STATS = par_theta * (n_ij - mui) / (par_theta + mui))
    dl_dbeta <- colSums(dl_dbeta_i)

    ## on theta
    dl_dtheta <- nSub + nSub * log(par_theta) +
        sum(1 / (par_theta + seq_n_ij - 1)) -
        sum((n_ij + par_theta) / (par_theta + mui)) - sum(log(mui_theta))

    ## on alpha
    part1 <- crossprod(1 / rho_0_ij, bMat_event)
    dl_dalpha_part2 <- sweep(dmu0_dalpha, MARGIN = 1, FUN = "*",
                             STATS = expXBeta[ind_cens] * n_ij_theta /
                                 mui_theta)
    part2 <- colSums(dl_dalpha_part2)
    dl_dalpha <-  part1 - part2
    ## return gradient
    - c(dl_dbeta, dl_dtheta, dl_dalpha)
}


rateReg_control <- function(Boundary.knots = NULL,
                            verbose = TRUE, ...)
{
    if (! isLogicalOne(verbose))
        stop("The option 'verbose' must be a logical value.", call. = FALSE)
    ## return
    list(control4rateReg = list(Boundary.knots = Boundary.knots,
                                verbose = verbose),
         control4optim = list(...))
}


rateReg_start <- function (beta, theta = 0.5, alpha, ..., nBeta_, nAlpha_)
{
    ## beta = starting value(s) for coefficients of covariates
    ## theta = starting value for random effects
    ## alpha = starting values for coefficients of baseline rate bases
    if (missing(beta)) {
        beta <- rep(0.1, nBeta_)
    } else if (length(beta) != nBeta_) {
        stop(wrapMessages(
            "The Number of starting values for covariate coefficients",
            "does not match with the specified formula."
        ), call. = FALSE)
    }
    if (! isNumOne(theta) || theta <= 0)
        stop("The parameter for gamma frailty must be a positive number.",
             call. = FALSE)
    if (missing(alpha))
        alpha <- rep(0.05, nAlpha_)
    ## return
    list(beta = beta, theta = theta, alpha = alpha)
}


## generate intervals from specified baseline pieces
nameBases <- function(df, spline)
{
    if (spline == "bSplines")
        return(paste0("B-spline", seq_len(df)))
    paste0("M-spline", seq_len(df))
}
