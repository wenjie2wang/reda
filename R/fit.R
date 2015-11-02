################################################################################
##
##   R package reda by Haoda Fu, Jun Yan, and Wenjie Wang
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


#' Fitting HEART Model
#'
#' \code{heart} returns fitted model results.
#' HEART model is a gamma frailty model with a piecewise constant
#' baseline rate function for recurrent event data. 
#' The model is named after the paper title of \emph{Fu et al. (2014)},  
#' Hypoglycemic Events Analysis via Recurrent Time-to-Event (HEART) Models
#'
#' In detail, function \code{heart} first constructs the design matrix from
#' the specified arguments: \code{formula}, \code{data}, \code{subset},
#' \code{na.action} and \code{constrasts} before fitting the HEART model.
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
#' function \code{\link{heart}} is called.
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
#' @return a \code{\link{heart-class}} object.
#' @references 
#' Fu, Haoda, Junxiang Luo, and Yongming Qu. (2014),
#' "Hypoglycemic Events Analysis via Recurrent Time-to-Event (HEART) Models," 
#' \emph{Journal of biopharmaceutical statistics}, 2014 Dec 1, Epub 2014 Dec 1.
#' @examples
#' library(reda)
#' ## data(simuDat)
#' heartFit <- heart(formula = Survr(ID, time, event) ~ x1 + group, 
#'                   data = simuDat, subset = ID %in% 75:125,
#'                   baselinePieces = seq(28, 168, length = 6))
#' ## str(heartFit)
#' show(heartFit) # or simply call 'heartFit'
#' summary(heartFit)
#' coef(heartFit)
#' confint(heartFit)
#' baseline(heartFit)
#' @seealso \code{\link{summary,heart-method}} \code{\link{coef,heart-method}}
#' \code{\link{confint,heart-method}} \code{\link{baseline,heart-method}}
#' \code{\link{mcf,heart-method}}
#' @importFrom methods new
#' @importFrom stats model.matrix nlm pnorm na.fail na.omit na.exclude na.pass
#' .getXlevels
#' @importFrom splines bs
#' @export
heart <- function(formula, df = NULL, knots = NULL, degree = 0,
                  data, subset, na.action, 
                  start = list(), control = list(), contrasts = NULL, ...) {
    ## record the function call to return
    Call <- match.call()
    ## arguments check
    if(missing(formula)) {
        stop("Argument 'formula' is required.")
    } 
    if(missing(data)) {
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
    ## get data.frame if na.action = na.pass
    mcall$na.action <- na.pass
    mf_na <- eval(mcall, parent.frame())
    mm_na <- stats::model.matrix(formula, data = mf_na,
                                 contrasts.arg = contrasts)
    ## number of covariates excluding intercept
    nbeta <- ncol(mm) - 1 
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
    
    ## baselinePieces
    if(missing(baselinePieces)) {
        baselinePieces <- as.numeric(max(dat$time))
    } 
    ## number of baseline pieces or rate functions
    nalpha <- length(baselinePieces)
    if (baselinePieces[nalpha] < max(dat$time)) {
        baselinePieces[nalpha] <- max(dat$time) + 1e-08
        warning("Baseline pieces is extended.")
    }
    ## friendly version of baseline pieces to print out
    print_blpieces <- int_baseline(baselinePieces = baselinePieces)
    attr(baselinePieces, "name") <- print_blpieces
    ## 'control' and 'start' values for 'nlm'
    startlist <- c(start, nbeta = nbeta, nalpha = nalpha)
    start <- do.call("heart_start", startlist)
    control <- do.call("heart_control", control)
    ini <- do.call("c", start)
    length_par <- length(ini)
    ## log likelihood
    fit <- stats::nlm(logL_heart, ini, data = dat, 
                      baselinePieces = baselinePieces, hessian = TRUE,
                      gradtol = control$gradtol, stepmax = control$stepmax,
                      steptol = control$steptol, iterlim = control$iterlim)
    
    est_beta <- matrix(NA, nrow = nbeta, ncol = 3)
    colnames(est_beta) <- c("coef", "se", "Pr(>|z|)")
    rownames(est_beta) <- covar_names
    
    se_vec <- sqrt(diag(solve(fit$hessian)))
    
    est_beta[, 1] <- fit$estimate[1:nbeta]
    est_beta[, 2] <- se_vec[1:nbeta]
    est_beta[, 3] <- 2 * stats::pnorm(-abs(est_beta[, 1]/est_beta[, 2]))
    
    est_theta <- matrix(NA, nrow = 1, ncol = 2)
    colnames(est_theta) <- c("theta", "se")
    est_theta[1, ] <- c(fit$estimate[nbeta + 1], se_vec[nbeta + 1])
    
    est_alpha <- matrix(NA, nrow = nalpha, ncol = 2)
    colnames(est_alpha) <- c("alpha", "se")
    rownames(est_alpha) <- print_blpieces 
    
    est_alpha[, 1] <- fit$estimate[(nbeta + 2):length_par]
    est_alpha[, 2] <- se_vec[(nbeta + 2):length_par]

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
    results <- methods::new("heart", 
                            call = Call, formula = formula, 
                            baselinePieces = baselinePieces,
                            estimates = list(beta = est_beta, 
                                             theta = est_theta, 
                                             alpha = est_alpha),
                            control = control,
                            start = start, 
                            na.action = na.action,
                            xlevels = .getXlevels(mt, mf),
                            contrasts = contrasts,
                            convergence = fit$code, 
                            fisher = fit$hessian)
    ## return
    results
}


## internal functions ==========================================================
whereT <- function(tt, baselinePieces) {
    ## designed to be used inside function 'apply', 'tt' is length one vector
    ## return the baseline segment number 'tt' belongs to
    min(which(tt <= baselinePieces))
}

## baseline rate function
rho_0 <- function(par_BaselinePW, Tvec, bKnots, degree) {
    ## if piecewise constant
    if (degree == 0) {
        indx <- apply(as.array(Tvec), 1, whereT, bKnots)
        return(par_BaselinePW[indx])  # function ends
    } 
    ## else spline with degree >= 1
    knots <- head(bKnots, -1)
    ## B-spline matrix including intercept
    ## with ncol: length(knots) + degree + 1
    bsmat <- bs(x = Tvec, knots = knots, degree = degree,
                intercept = TRUE, Boundary.knots = c(0, tail(bKnots, 1)))
    ## return
    crossprod(t(bsmat), par_BaselinePW)
}

mu0 <- function(par_BaselinePW, baselinePieces, Tvec) {
    indx <- apply(as.array(Tvec), 1, whereT, baselinePieces)
    BL_segments <- c(baselinePieces[1], diff(baselinePieces))
    ## The MCF at each time point  
    CumMean_Pieces <- diffinv(BL_segments * par_BaselinePW)[-1]  
    ## return
    CumMean_Pieces[indx] - (baselinePieces[indx] - Tvec) * par_BaselinePW[indx]
}

dmu0_alpha <- function(tt, baselinePieces) {
    indx <- min(which(tt <= baselinePieces))
    ## BL_segments 
    value <- diff(c(0, baselinePieces))
    n_pieces <- length(baselinePieces)
    if (indx == n_pieces) {
        value[n_pieces] <- ifelse(n_pieces == 1, tt, 
                                  tt - baselinePieces[n_pieces - 1])
    } else if (indx > 1) {
        value[(indx + 1):n_pieces] <- 0
        value[indx] <- tt - baselinePieces[indx - 1]
    } else {
        value[(indx + 1):n_pieces] <- 0
        value[indx] <- tt
    }
    ## return
    value
}

## generate intervals from specified baseline pieces
int_baseline <- function(baselinePieces){
    nalpha <- length(baselinePieces)
    intervals <- rep(NA, nalpha)
    intervals[1] <- paste0("[0, ", baselinePieces[1], "]", sep = "")
    if (nalpha > 1) {
        for(i in 2:nalpha){
            intervals[i] <- paste0("(", baselinePieces[i - 1], ", ", 
                                   baselinePieces[i], "]", sep = "")
        }
    }
    ## return
    intervals
}

## compute negative log likelihood
logL_heart <- function(par, data, baselinePieces) {
    nbeta <- ncol(data) - 3
    ## par = \THETA in the paper
    par_beta <- par[1 : nbeta]
    par_theta <- par[nbeta + 1]
    par_alpha <- par[(nbeta + 2) : length(par)]
    m <- length(unique(data$ID))
    expXBeta <- exp(crossprod(t(as.matrix(data[, 4 : (3 + nbeta)])),
                              as.matrix(par_beta)))
    ## index for event and censoring
    ind_event <- data$event == 1
    ind_cens <- data$event == 0
    ## rate function
    rho_0_ij <- rho_0(par_BaselinePW = par_alpha,
                      baselinePieces = baselinePieces, 
                      data$time[ind_event])
    rho_i <- expXBeta[ind_event] * rho_0_ij
    rho_i[rho_i < 1e-100] <- 1e-100
    sum_log_rho_i <- sum(log(rho_i))
    ## these codes to make sure that the order will not change 
    ## if the patient ID is not ordered
    n_ij <- table(data$ID)[order(unique(data$ID))] - 1  
    ## if there is a subject with 0 event, 
    ## the sequence will not be generated for this subject
    theta_j_1 <- par_theta + sequence(n_ij) - 1  
    theta_j_1[theta_j_1 < 1e-100] <- 1e-100
    sum_log_theta_j_1 <- sum(log(theta_j_1))
    mu0i <- mu0(par_BaselinePW = par_alpha, baselinePieces = baselinePieces, 
                data$time[ind_cens])
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
    dl_dbeta <- apply(diag((n_ij - mui)/(par_theta + mui) * par_theta) %*% 
                      as.matrix(data[ind_cens, 4:(3 + nbeta)]), 2, sum)
    dl_dtheta <- m + m * log(par_theta) + 
        sum(1/(par_theta + sequence(n_ij) - 1)) - 
        sum((n_ij + par_theta)/(par_theta + mui)) - sum(log(mui_theta))
    indx <- apply(as.array(data$time[ind_event]), 1, 
                  whereT, baselinePieces)
    if (length(unique(indx)) < length(par_alpha)) {
        stop("Some segements have zero events!")
    }
    indx_taui <- apply(as.array(data$time[ind_cens]), 1, 
                       whereT, baselinePieces)
    ## reform dimension by 'array' for one-piece baseline 
    dim_n1 <- length(baselinePieces)
    dim_n2 <- length(data$time[ind_cens])
    tempart2 <- array(apply(array(data$time[ind_cens]), 1, 
                            dmu0_alpha, baselinePieces), c(dim_n1, dim_n2))
    dl_dalpha_part2 <- diag((n_ij + par_theta) / (par_theta + mui) * 
                            expXBeta[ind_cens]) %*% t(tempart2)
    dl_dalpha <- 1 / par_alpha * table(indx) - apply(dl_dalpha_part2, 2, sum)
    attr(negLH, "gradient") <- -c(dl_dbeta, dl_dtheta, dl_dalpha)
    ## return
    negLH
}

heart_control <- function (gradtol = 1e-6, stepmax = 1e5, 
                           steptol = 1e-6, iterlim = 1e2) {
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
    ## return
    list(gradtol = gradtol, stepmax = stepmax, 
         steptol = steptol, iterlim = iterlim)
}

heart_start <- function(beta, theta = 0.5, alpha, nbeta, nalpha) {
    ## beta = starting value(s) for coefficients of covariates
    ## theta = starting value for random effects
    ## alpha = starting values for piece-wise baseline rate functions
    if (missing(beta)) {
        beta <- rep(1, nbeta)
    } else if (length(beta) != nbeta) {
        stop(paste("number of starting values for coefficients of covariates",
                   "does not match with the specified formula"))
    }
    if (theta <= 0) {
        stop("value of parameter for random effects must be > 0")
    }
    if (missing(alpha)) {
        alpha <- rep(0.15, nalpha)
    }
    ## return
    list(beta = beta, theta = theta, alpha = alpha)
}

