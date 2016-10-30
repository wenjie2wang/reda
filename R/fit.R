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


## collation after class.R
##' @include class.R
NULL


##' Fit Recurrent Events Regression Based on Counts and Rate Function
##'
##' The default model is the gamma frailty model with one piece constant
##' baseline rate function, which is equivalent to negative binomial regression
##' of the same shape and rate parameter in gamma prior.
##' Spline and piecewise constant baseline rate function can be
##' specified and applied to model fitting instead.
##' \code{rateReg} returns the fitted model
##' through a \code{\link{rateReg-class}} object.
##'
##' @details
##' Function \code{\link{Survr}} in the formula response first checks
##' the dataset and will report an error if the dataset does not
##' fall into recurrent event data framework.
##' Subject's ID is pinpointed if its observation violates any checking rule.
##' See \code{\link{Survr}} for all the checking rules.
##'
##' Function \code{rateReg} first constructs the design matrix from
##' the specified arguments: \code{formula}, \code{data}, \code{subset},
##' \code{na.action} and \code{constrasts} before model fitting.
##' The constructed design matrix will be checked again to
##' fit the recurrent event data framework
##' if any observation with missing covariates is removed.
##'
##' The model fitting process involves minimization of negative log
##' likelihood function, which calls function \code{\link[stats]{nlm}}
##' from package \pkg{stats} internally.
##' \code{help(nlm)} for more details.
##'
##' The argument \code{start} is an optional list
##' that allows users to specify the initial guess for
##' the parameter values for the minimization of
##' negative log likelihood function.
##' The available numeric vector elements in the list include
##' \itemize{
##'     \item \code{beta}: Coefficient(s) of covariates,
##'         set to be 0.1 by default.
##'     \item \code{theta}: Parameter of frailty random effect,
##'         set to be 0.5 by default.
##'     \item \code{alpha}: Coefficient(s) of baseline rate function,
##'         set to be 0.05 by default.
##' }
##' The argument \code{control} is an optional list
##' that allows users to control the process of minimization of
##' negative log likelihood function and to specify the boundary knots,
##' intercept for baseline rate function.
##' The available elements in the list include
##' \itemize{
##'     \item \code{gradtol}: A positive scalar giving the tolerance at
##'         which the scaled gradient is considered close enough to zero
##'         to terminate the algorithm. The default value is 1e-6.
##'     \item \code{stepmax}: A positive scalar that gives the maximum
##'         allowable scaled step length. The default value is 1e5.
##'     \item \code{steptol}: A positive scalar providing the minimum
##'         allowable relative step length. The default value is 1e-6.
##'     \item \code{iterlim}: A positive integer specifying the maximum
##'         number of iterations to be performed before
##'         the program is terminated. The default value is 1e2.
##'     \item \code{Boundary.knots}: A length-two numeric vector to specify
##'         the boundary knots for baseline rate funtion. By default,
##'         the left boundary knot is zero and the right one takes the
##'         largest censoring time from data.
##'     \item \code{intercept}: A logical value specifying whether
##'         intercept is included in spline baseline rate function.
##'         For piecewise constatn baseline (\code{df}=0), the specified
##'         value would be neglected. The default value
##'         is \code{TRUE}, i.e. the intercept is included.
##' }
##'
##' @param formula \code{Survr} object produced by function \code{\link{Survr}}.
##' @param df An optional nonnegative integer to specify the degree of freedom
##' of baseline rate function. If argument \code{knots} or \code{degree} are
##' specified, \code{df} will be neglected whether it is specified or not.
##' @param knots An optional numeric vector that represents all the internal
##' knots of baseline rate function.
##' The default is \code{NULL}, representing no any internal knots.
##' @param degree An optional nonnegative integer to specify the degree of
##' spline bases.
##' @param data An optional data frame, list or environment containing
##' the variables in the model.  If not found in data, the variables are taken
##' from \code{environment(formula)}, usually the environment from which
##' function \code{\link{rateReg}} is called.
##' @param subset An optional vector specifying a subset of observations
##' to be used in the fitting process.
##' @param na.action A function that indicates what should the procedure
##' do if the data contains \code{NA}s.  The default is set by the
##' na.action setting of \code{\link[base]{options}}.
##' The "factory-fresh" default is \code{\link[stats]{na.omit}}.
##' Other possible values inlcude \code{\link{na.fail}},
##' \code{\link{na.exclude}}, and \code{\link{na.pass}}.
##' \code{help(na.fail)} for details.
##' @param start An optional list of starting values for the parameters
##' to be estimated in the model.  See more in section details.
##' @param control An optional list of parameters to control the
##' maximization process of negative log likelihood function
##' and adjust the baseline rate function.
##' See more in section details.
##' @param contrasts An optional list, whose entries are values
##' (numeric matrices or character strings naming functions) to be used
##' as replacement values for the contrasts replacement function and
##' whose names are the names of columns of data containing factors.
##' See \code{contrasts.arg} of \code{\link[stats]{model.matrix.default}}
##' for details.
##' @param ... Other arguments for future usage.
##' @return A \code{\link{rateReg-class}} object, whose slots include
##' \itemize{
##'     \item \code{call}: Function call of \code{rateReg}.
##'     \item \code{formula}: Formula used in the model fitting.
##'     \item \code{nObs}: Number of observations.
##'     \item \code{knots}: Internal knots specified for the baseline
##'         rate function.
##'     \item \code{Boundary.knots}: Boundary knots specified for the baseline
##'         rate function.
##'     \item \code{degree}: Degree of spline bases specified in baseline
##'         rate function.
##'     \item \code{df}: Degree of freedom of the model specified.
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
##'         \code{\link[stats]{nlm}}, which is an integer indicating why the
##'         optimization process terminated. \code{help(nlm)} for details.
##'     \item \code{logL}: Log likelihood of the fitted model.
##'     \item \code{fisher}: Observed Fisher information matrix.
##' }
##'
##' @references
##' Fu, H., Luo, L., & Qu Y. (2014). Hypoglycemic Events Analysis via
##' Recurrent Time-to-Event (HEART) Models.
##' \emph{Journal of biopharmaceutical statistics}, Epub 2014 Dec 1.
##' @examples
##' library(reda)
##'
##' ## constant rate function
##' constFit <- rateReg(Survr(ID, time, event) ~ group + x1,
##'                     data = simuDat, subset = ID %in% 1:50)
##'
##' ## 6 pieces' piecewise constant rate function
##' piecesFit <- rateReg(Survr(ID, time, event) ~ group + x1,
##'                      data = simuDat, subset = ID %in% 1:50,
##'                      knots = seq(28, 140, by = 28))
##'
##' ## fit rate function with cubic spline
##' splineFit <- rateReg(Survr(ID, time, event) ~ group + x1,
##'                      data = simuDat, subset = ID %in% 1:50,
##'                      knots = c(56, 84, 112), degree = 3)
##'
##' ## brief summary of fitted models
##' constFit
##' piecesFit
##' splineFit
##'
##' ## more specific summary
##' summary(constFit)
##' summary(piecesFit)
##' summary(splineFit)
##'
##' ## model selection based on AIC or BIC
##' AIC(constFit, piecesFit, splineFit)
##' BIC(constFit, piecesFit, splineFit)
##'
##' ## estimated covariate coefficients
##' coef(piecesFit)
##' coef(splineFit)
##'
##' ## confidence intervals for covariate coefficients
##' confint(piecesFit)
##' confint(splineFit, "x1", 0.9)
##' confint(splineFit, 1, 0.975)
##'
##' ## estimated coefficients for baseline rate function
##' baseRate(piecesFit)
##' baseRate(splineFit)
##'
##' ## estimated baseline mean cumulative function (MCF) from a fitted model
##' piecesMcf <- mcf(piecesFit)
##' plotMcf(piecesMcf, conf.int = TRUE, col = "blueviolet") +
##'     ggplot2::xlab("Days") + ggplot2::theme_bw()
##'
##' ## estimated MCF for given new data
##' newDat <- data.frame(x1 = rep(0, 2), group = c("Treat", "Contr"))
##' splineMcf <- mcf(splineFit, newdata = newDat, groupName = "Group",
##'                  groupLevels = c("Treatment", "Control"))
##' plotMcf(splineMcf, conf.int = TRUE, lty = c(1, 5)) +
##'     ggplot2::xlab("Days") + ggplot2::theme_bw()
##'
##' @seealso
##' \code{\link{summary,rateReg-method}} for summary of fitted model;
##' \code{\link{coef,rateReg-method}} for estimated covariate coefficients;
##' \code{\link{confint,rateReg-method}} for confidence interval of
##' covariate coefficients;
##' \code{\link{baseRate,rateReg-method}} for estimated coefficients of baseline
##' rate function;
##' \code{\link{mcf,rateReg-method}} for estimated MCF from a fitted model;
##' \code{\link{plotMcf,rateRegMcf-method}} for plotting estimated MCF.
##' @importFrom splines2 bSpline ibs
##' @importFrom stats na.fail na.omit na.exclude na.pass .getXlevels
##' @export
rateReg <- function(formula, df = NULL, knots = NULL, degree = 0L,
                    data, subset, na.action, start = list(), control = list(),
                    contrasts = NULL, ...) {
    ## record the function call to return
    Call <- match.call()

    ## arguments check
    if (missing(formula))
        stop("Argument 'formula' is required.")
    if (missing(data))
        data <- environment(formula)
    if (! with(data, inherits(eval(Call[[2]][[2]]), "Survr")))
        stop("Response in formula must be a survival recurrent object.")

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
    if ((nBeta <- ncol(mm) - 1L) <= 0)
        stop("Covariates must be specified in formula.")

    ## covariates' names
    covar_names <- colnames(mm)[- 1L]
    ## data
    dat <- as.data.frame(cbind(mf[, 1L][, seq_len(3L)], mm[, - 1L]))
    colnames(dat) <- c("ID", "time", "event", covar_names)
    nObs <- nrow(dat)

    ## check the impact caused by missing value
    ## if there is missing value removed
    if (nrow(mm_na) > nObs) {
        ## recover original ID names for possible pin-point
        idFactor <- with(data, attr(eval(Call[[2L]][[2L]]), "ID"))
        attr(dat, "ID") <- factor(levels(idFactor)[dat$ID],
                                  levels = levels(idFactor))
        message("Observations with missing covariate value are removed.")
        message("Checking the new dataset again ... ", appendLF = FALSE)
        check_Survr(dat)
        message("done.")
    }

    ## 'control' for optimization and splines' boundary knots
    control <- do.call("rateReg_control", control)
    Boundary.knots <- if (is.null(control$Boundary.knots)) {
                          c(0, max(dat$time))
                      } else {
                          control$Boundary.knots
                      }

    ## generate knots if knots is unspecified
    iMat <- splines2::ibs(x = dat$time, df = df, knots = knots,
                          degree = degree, intercept = TRUE,
                          Boundary.knots = Boundary.knots)
    bMat <- attr(iMat, "bsMat")

    ## update df, knots, degree, and Boundary.knots
    knots <- as.numeric(attr(iMat, "knots"))
    df <- (degree <- attr(iMat, "degree")) + length(knots) + 1L
    Boundary.knots <- attr(iMat, "Boundary.knots")

    ## name each basis for alpha output
    alphaName <- nameBases(degree = degree, df = df, knots = knots,
                           Boundary.knots = Boundary.knots)

    ## start' values for 'nlm'
    startlist <- c(start, list(nBeta_ = nBeta, nAlpha_ = df))
    start <- do.call("rateReg_start", startlist)
    ini <- do.call("c", start)
    length_par <- length(ini)

    ## check whether the knots are reasonable
    if (any(colSums(bMat[dat$event == 1, , drop = FALSE]) == 0))
        stop("Some time segement contains zero events.")

    ## log likelihood
    fit <- stats::nlm(logL_rateReg, ini, dat = dat, bMat = bMat, iMat = iMat,
                      hessian = TRUE,
                      gradtol = control$gradtol, stepmax = control$stepmax,
                      steptol = control$steptol, iterlim = control$iterlim)

    ## estimates for beta
    est_beta <- matrix(NA, nrow = nBeta, ncol = 5)
    colnames(est_beta) <- c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")
    rownames(est_beta) <- covar_names

    se_vec <- sqrt(diag(solve(fit$hessian)))
    est_beta[, 1] <- fit$estimate[1:nBeta]
    est_beta[, 2] <- exp(est_beta[, 1])
    est_beta[, 3] <- se_vec[1:nBeta]
    est_beta[, 4] <- est_beta[, 1] / est_beta[, 3]
    est_beta[, 5] <- 2 * stats::pnorm(- abs(est_beta[, 4]))

    ## estimates for theta
    est_theta <- matrix(NA, nrow = 1, ncol = 2)
    colnames(est_theta) <- c("parameter", "se")
    rownames(est_theta) <- "Frailty"
    est_theta[1, ] <- c(fit$estimate[nBeta + 1], se_vec[nBeta + 1])

    ## estimates for alpha
    est_alpha <- matrix(NA, nrow = df, ncol = 2)
    colnames(est_alpha) <- c("coef", "se(coef)")
    rownames(est_alpha) <- alphaName
    est_alpha[, 1] <- fit$estimate[(nBeta + 2):length_par]
    est_alpha[, 2] <- se_vec[(nBeta + 2):length_par]

    ## output: na.action
    na.action <- if (is.null(attr(mf, "na.action"))) {
                     options("na.action")[[1]]
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

    ## results to return
    results <- methods::new("rateReg",
                            call = Call,
                            formula = formula,
                            nObs = nObs,
                            knots = knots,
                            Boundary.knots = Boundary.knots,
                            degree = degree,
                            df = df,
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


### internal functions =========================================================
## compute negative log likelihood
logL_rateReg <- function(par, dat, bMat, iMat) {

    ## number of covariates, possibly zero
    nBeta <- ncol(dat) - 3L

    ## par = \THETA in the paper
    par_theta <- max(par[nBeta + 1L], .Machine$double.eps)
    par_alpha <- par[(nBeta + 2L) : length(par)]
    m <- length(unique(dat$ID))
    xMat <- as.matrix(dat[, 4L : (3L + nBeta)])
    expXBeta <- as.vector(exp(xMat %*% as.matrix(par[seq_len(nBeta)])))

    ## index for event and censoring
    ind_event <- dat$event == 1L
    ind_cens <- ! ind_event

    ## basis matrix at event times
    bMat_event <- bMat[ind_event, , drop = FALSE]

    ## baseline rate function
    rho_0_ij <- pmax(as.vector(bMat_event %*% par_alpha), .Machine$double.eps)
    rho_i <- pmax(expXBeta[ind_event] * rho_0_ij, .Machine$double.eps)
    sum_log_rho_i <- sum(log(rho_i))

    ## n_ij: number of event for each subject
    ## the following code makes sure the order will not change
    ## if the patient ID is not ordered
    n_ij <- table(dat$ID)[order(unique(dat$ID))] - 1L

    ## if there is a subject with 0 event,
    ## the sequence will not be generated for this subject
    ## note that sequence(0L) leads to integer(0)
    seq_n_ij <- sequence(n_ij)
    theta_j_1 <- par_theta + seq_n_ij - 1
    sum_log_theta_j_1 <- sum(log(theta_j_1))

    ## baseline mcf, integral of rho_0 that involves censoring time tau
    mu0i <- as.vector(iMat[ind_cens, , drop = FALSE] %*% par_alpha)
    mui <- mu0i * expXBeta[ind_cens]
    mui_theta <- pmax(par_theta + mui, .Machine$double.eps)
    sum_log_theta_mui <- sum((n_ij_theta <- n_ij + par_theta) * log(mui_theta))

    ## log-likelihood function
    logLH <- m * par_theta * log(par_theta) + sum_log_rho_i +
        sum_log_theta_j_1 - sum_log_theta_mui

    ## or use constrOptim?
    penal <- ifelse(par_theta < 0 || min(par_alpha) < 0, 1e50, 0)
    negLH <- - logLH + penal

    ## Calculate the gradient
    ## on beta, vector
    xMat_i <- xMat[ind_cens, , drop = FALSE]
    dl_dbeta_i <- sweep(x = as.matrix(xMat_i), MARGIN = 1, FUN = "*",
                        STATS = par_theta * (n_ij - mui) / (par_theta + mui))
    dl_dbeta <- colSums(dl_dbeta_i)

    ## on theta
    dl_dtheta <- m + m * log(par_theta) + sum(1 / (par_theta + seq_n_ij - 1)) -
        sum((n_ij + par_theta) / (par_theta + mui)) - sum(log(mui_theta))

    ## on alpha
    part1 <- crossprod(1 / rho_0_ij, bMat_event)
    dmu0_dalpha <- iMat[ind_cens, , drop = FALSE]
    dl_dalpha_part2 <- sweep(dmu0_dalpha, MARGIN = 1, FUN = "*",
                             STATS = expXBeta[ind_cens] * n_ij_theta /
                                 mui_theta)
    part2 <- colSums(dl_dalpha_part2)
    dl_dalpha <-  part1 - part2
    ## return gradient as one attribute
    attr(negLH, "gradient") <- - c(dl_dbeta, dl_dtheta, dl_dalpha)
    negLH
}


rateReg_control <- function(gradtol = 1e-6, stepmax = 1e5,
                            steptol = 1e-6, iterlim = 1e2,
                            Boundary.knots = NULL, ...) {
    ## controls for function stats::nlm
    if (! is.numeric(gradtol) || gradtol <= 0)
        stop("Value of 'gradtol' must be > 0.")
    if (! is.numeric(stepmax) || stepmax <= 0)
        stop("Value of 'stepmax' must be > 0.")
    if (! is.numeric(steptol) || steptol <= 0)
        stop("Value of 'steptol' must be > 0.")
    if (! is.numeric(iterlim) || iterlim <= 0)
        stop("Maximum number of iterations must be > 0.")
    ## return
    list(gradtol = gradtol, stepmax = stepmax,
         steptol = steptol, iterlim = iterlim,
         Boundary.knots = Boundary.knots)
}


rateReg_start <- function (beta, theta = 0.5, alpha, ..., nBeta_, nAlpha_) {
    ## beta = starting value(s) for coefficients of covariates
    ## theta = starting value for random effects
    ## alpha = starting values for coefficients of baseline rate bases
    if (missing(beta)) {
        beta <- rep(0.1, nBeta_)
    } else if (length(beta) != nBeta_) {
        stop(paste("Number of starting values for coefficients of covariates",
                   "does not match with the specified formula."))
    }
    if (theta <= 0)
        stop("Value of parameter for random effects must be > 0.")
    if (missing(alpha))
        alpha <- rep(0.05, nAlpha_)
    ## return
    list(beta = beta, theta = theta, alpha = alpha)
}


## generate intervals from specified baseline pieces
nameBases <- function(degree, df, knots, Boundary.knots) {
    leftVec <- c(Boundary.knots[1L], knots)
    rightVec <- c(knots, Boundary.knots[2L])
    if (! degree)
        return(paste0("(", leftVec, ", ", rightVec, "]"))
    ## else degree > 0
    ## return
    paste0("B-spline.", seq_len(df))
}
