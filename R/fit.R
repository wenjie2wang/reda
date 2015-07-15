################################################################################
##
##   R package heart by Haoda Fu, Jun Yan, and Wenjie Wang
##   Copyright (C) 2015
##
##   This file is part of the R package heart.
##
##   The R package heart is free software: you can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package heart is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package heart. If not, see <http://www.gnu.org/licenses/>.
##
################################################################################

## functions to export =========================================================
## create S4 Class called "heart" for heart object from function heart
#' @export
setClass(Class = "heart", 
         slots = c(call = "call", 
                   formula = "formula", 
                   baselinepieces = "numeric",
                   estimates = "list",
                   control = "list",
                   start = "list",
                   na.action = "character",
                   contrasts = "list",
                   convergence = "integer", 
                   hessian = "matrix"))

#' Fitting HEART Model
#'
#' \code{heart} returns fitted model results.
#'
#' HEART model is a piece-wise Gamma frailty model for recurrent events. 
#' The model is named after the paper title of \emph{Fu et al. (2014)},  
#' Hypoglycemic Events Analysis via Recurrent Time-to-Event (HEART) Models
#'
#' @param formula Survr object from function \code{Survr} in package survrec.
#' @param baselinepieces an optional numeric vector consisting of
#' all the right endpoints of baseline pieces.
#' @param data an optional data frame, list or environment containing
#' the variables in the model.  If not found in data, the variables are taken 
#' from \code{environment(formula)}, usually the environment from which 
#' \code{heart} is called.
#' @param subset an optional vector specifying a subset of observations 
#' to be used in the fitting process.
#' @param na.action a function which indicates what should the procedure do 
#' if the data contains NAs.  The default is set by the 
#' na.action setting of \code{\link[base]{options}} and is na.fail if that is 
#' not set.  The "factory-fresh" default is \code{\link[stats]{na.omit}}.
#' Another possible value is NULL, no action.  
#' Value \code{\link[stats]{na.exclude}} can be useful. 
#' @param start an optional list of starting values for the parameters
#' to be estimated in the model.
#' @param control an optional list of parameters for controlling the likelihood 
#' function maximization process. For more detail, users may \code{help(nlm)}. 
#' @param contrasts an optional list, whose entries are values 
#' (numeric matrices or character strings naming functions) to be used 
#' as replacement values for the contrasts replacement function and 
#' whose names are the names of columns of data containing factors.
#' See the \code{contrasts.arg} of \code{model.matrix.default} for more detail.
#' @return a heart object.
#' @references 
#' Fu, Haoda, Junxiang Luo, and Yongming Qu. (2014),
#' "Hypoglycemic Events Analysis via Recurrent Time-to-Event (HEART) Models," 
#' \emph{Journal of biopharmaceutical statistics}, 2014 Dec 1, Epub 2014 Dec 1.
#' @examples
#' data(simuDat)
#' heartfit <- heart(formula = Survr(ID, time, event) ~ X1 + group, 
#'                   data = simuDat, baselinepieces = seq(28, 168, length = 6))
#' str(heartfit)
#' show(heartfit) # or simply call heartfit
#' summary(heartfit)
#' coef(heartfit)
#' confint(heartfit)
#' baseline(heartfit)
#' plot_MCF(heartfit)
#' @import survrec
#' @export
heart <- function(formula, baselinepieces, data, subset, na.action, 
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
    stop("'formula' must be a survival recurrent object.")
  }
  ## Prepare data: ID, time, event ~ X(s)
  mcall <- match.call(expand.dots = FALSE)
  mmcall <- match(c("formula", "data", "subset", "na.action"), names(mcall), 0L)
  mcall <- mcall[c(1L, mmcall)]
  ## drop unused levels in factors 
  mcall$drop.unused.levels <- TRUE
  mcall[[1L]] <- quote(stats::model.frame)
  mf <- eval(mcall, parent.frame())
  mm <- stats::model.matrix(formula, data = mf)
  ## number of covariates excluding intercept
  nbeta <- ncol(mm) - 1 
  ## covariates' names
  covar_names <- base::colnames(mm)[-1]
  ## data 
  dat <- as.data.frame(cbind(mf[, 1][, 1:3], mm[, -1]))
  colnames(dat) <- c("ID", "time", "event", covar_names)
  ## baselinepieces
  if(missing(baselinepieces)) {
    warning("Baseline pieces are splitted by median event time.")
    baselinepieces <- as.vector(round(quantile(dat$time, c(0.5, 1)), 
                                digits = 1))
  } 
  ## number of baseline pieces or rate functions
  nalpha <- length(baselinepieces)
  if (baselinepieces[nalpha] < max(dat$time)) {
    baselinepieces[nalpha] <- max(dat$time) + 1e-08
    warning("Baseline pieces is extended.")
  }
  ## friendly version of baseline pieces to print out
  print_blpieces <- int_baseline(baselinepieces = baselinepieces)
  attr(baselinepieces, "name") <- print_blpieces
  ## 'control' and 'start' values for 'nlm'
  startlist <- c(start, nbeta = nbeta, nalpha = nalpha)
  start <- do.call("heart_start", startlist)
  control <- do.call("heart_control", control)
  ini <- do.call("c", start)
  length_par <- length(ini)
  ## log likelihood
  fit <- nlm(logL_heart, ini, data = dat, 
             baselinepieces = baselinepieces, hessian = TRUE,
             gradtol = control$gradtol, stepmax = control$stepmax,
             steptol = control$steptol, iterlim = control$iterlim)

  est_beta <- matrix(NA, nrow = nbeta, ncol = 3)
  colnames(est_beta) <- c("coef", "se", "Pr(>|z|)")
  rownames(est_beta) <- covar_names
  
  se_vec <- sqrt(diag(solve(fit$hessian)))
  
  est_beta[, 1] <- fit$estimate[1:nbeta]
  est_beta[, 2] <- se_vec[1:nbeta]
  est_beta[, 3] <- 2 * pnorm(-abs(est_beta[, 1]/est_beta[, 2]))
  
  est_theta <- matrix(NA, nrow = 1, ncol = 2)
  colnames(est_theta) <- c("theta", "se")
  est_theta[1, ] <- c(fit$estimate[nbeta + 1], se_vec[nbeta + 1])
  
  est_alpha <- matrix(NA, nrow = nalpha, ncol = 2)
  colnames(est_alpha) <- c("alpha", "se")
  rownames(est_alpha) <- print_blpieces 
  
  est_alpha[, 1] <- fit$estimate[(nbeta + 2):length_par]
  est_alpha[, 2] <- se_vec[(nbeta + 2):length_par]
  if (is.null(attr(mf, "na.action"))) {
    na.action <- options("na.action")[[1]]
  } else {
    na.action <- paste("na", class(attr(mf, "na.action")), sep = ".")
  }
  if (is.null(contrasts)){
    contrasts <- options("contrasts")
  }
  ## results to return
  results <- new("heart", 
                 call = Call, formula = formula, 
                 baselinepieces = baselinepieces,
                 estimates = list(beta = est_beta, 
                                  theta = est_theta, 
                                  alpha = est_alpha),
                 control = control,
                 start = start, 
                 na.action = na.action,
                 contrasts = contrasts,
                 convergence = fit$code, 
                 hessian = fit$hessian)
  ## return
  results
}


## internal functions ==========================================================
whereT <- function(tt, baselinepieces) {
  ## return
  min(which(tt <= baselinepieces))
}

rho_0 <- function(par_BaselinePW, baselinepieces, Tvec) {
  indx <- apply(as.array(Tvec), 1, whereT, baselinepieces)
  ## return
  par_BaselinePW[indx]
}

mu0 <- function(par_BaselinePW, baselinepieces, Tvec) {
  indx <- apply(as.array(Tvec), 1, whereT, baselinepieces)
  BL_segments <- c(baselinepieces[1], diff(baselinepieces))
  ## The CMF at each time point  
  CumMean_Pieces <- diffinv(BL_segments * par_BaselinePW)[-1]  
  ## return
  CumMean_Pieces[indx] - (baselinepieces[indx] - Tvec) * par_BaselinePW[indx]
}

dmu0_alpha <- function(tt, baselinepieces) {
  BL_segments <- c(baselinepieces[1], diff(baselinepieces))
  indx <- min(which(tt <= baselinepieces))
  value <- BL_segments
  n_pieces <- length(baselinepieces)
  if (indx == n_pieces) {
    value[n_pieces] <- tt - baselinepieces[n_pieces - 1]
  } else if (indx > 1) {
    value[(indx + 1):n_pieces] <- 0
    value[indx] <- tt - baselinepieces[indx - 1]
  } else {
    value[(indx + 1):n_pieces] <- 0
    value[indx] <- tt
  }
  ## return
  value
}

## generate intervals from specified baseline pieces
int_baseline <- function(baselinepieces){
  nalpha <- length(baselinepieces)
  intervals <- rep(NA, nalpha)
  intervals[1] <- paste0("[0, ", baselinepieces[1], "]", sep = "")
  if (nalpha > 1) {
    for(i in 2:nalpha){
      intervals[i] <- paste0("(", baselinepieces[i - 1], ", ", 
                             baselinepieces[i], "]", sep = "")
    }
  }
  ## return
  intervals
}

## compute log likelihood
logL_heart <- function(par, data, baselinepieces) {
  nbeta <- ncol(data) - 3
  ## par = \THETA in the paper
  par_beta <- par[1:nbeta]
  par_theta <- par[nbeta + 1]
  par_alpha <- par[(nbeta + 2):length(par)]
  m <- length(unique(data$ID))
  expXBeta <- exp(as.matrix(data[, 4:(3 + nbeta)]) %*% (as.matrix(par_beta)))
  rho_0_ij <- rho_0(par_BaselinePW = par_alpha, baselinepieces = baselinepieces, 
                    data$time[data$event == 1])
  rho_i <- expXBeta[data$event == 1] * rho_0_ij
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
  mu0i <- mu0(par_BaselinePW = par_alpha, baselinepieces = baselinepieces, 
              data$time[data$event == 0])
  mui <- mu0i * expXBeta[data$event == 0]
  mui_theta <- par_theta + mui
  mui_theta[mui_theta < 1e-100] <- 1e-100
  sum.log_theta_mui <- sum((n_ij + par_theta) * log(mui_theta))
  if (par_theta < 1e-100) {
    par_theta <- 1e-100
  }
  logLH <- m * par_theta * log(par_theta) + sum_log_rho_i + 
    sum_log_theta_j_1 - sum.log_theta_mui
  penal <- ifelse(par_theta < 0 | min(par_alpha) < 0, 1e+50, 0)
  negLH <- -logLH + penal
  ## Calculate the gradient
  dl_dbeta <- apply(diag((n_ij - mui)/(par_theta + mui) * par_theta) %*% 
                      as.matrix(data[data$event == 0, 4:(3 + nbeta)]), 2, sum)
  dl_dtheta <- m + m * log(par_theta) + 
    sum(1/(par_theta + sequence(n_ij) - 1)) - 
    sum((n_ij + par_theta)/(par_theta + mui)) - sum(log(mui_theta))
  indx <- apply(as.array(data$time[data$event == 1]), 1, 
                whereT, baselinepieces)
  if (length(unique(indx)) < length(par_alpha)) {
    stop("Some segements have zero events!")
  }
  indx_taui <- apply(as.array(data$time[data$event == 0]), 1, 
                     whereT, baselinepieces)
  dl_dalpha_part2 <- diag((n_ij + par_theta) / (par_theta + mui) * 
                            expXBeta[data$event == 0]) %*% 
    t(apply(array(data$time[data$event == 0]), 1, dmu0_alpha, baselinepieces))
  dl_dalpha <- 1 / par_alpha * table(indx) - apply(dl_dalpha_part2, 2, sum)
  attr(negLH, "gradient") <- -c(dl_dbeta, dl_dtheta, dl_dalpha)
  ## return
  negLH
}

heart_control <- function (gradtol = 1e-6, stepmax = 1e5, 
                           steptol = 1e-6, iterlim = 1e2) {
  
  ## controls for function stats::nlm
  m <- match.call(expand.dots = FALSE)
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
    stop("number of starting values for coefficients of covariates 
           does not match with the specified formula")
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


