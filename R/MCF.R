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

#' Empirical Mean Cumulative Function (MCF)
#' 
#' \code{sample_MCF} computes empirical mean cumulative function (MCF) on every
#' time point from recurrent event sample data.  
#' It does not assume any particular underlying model. 
#' @usage sample_MCF(formula, data, subset, na.action, ...)
#' @param formula Survr object from function \code{Survr} in package survrec.
#' The covariate specified in the rhs of the formula can either be 1 or 
#' any one factor variable in the data.  The former computes the overall 
#' empirical MCF from sample.  The latter computes the empirical MCF for each 
#' level of the factor variable specified respectively.
#' @param data an optional data frame, list or environment containing
#' the variables in the model.  If not found in data, the variables are taken 
#' from \code{environment(formula)}, usually the environment from which 
#' the function is called.
#' @param subset an optional vector specifying a subset of observations 
#' to be used in the fitting process.
#' @param na.action a function which indicates what should the procedure do 
#' if the data contains NAs.  The default is set by the 
#' na.action setting of \code{\link[base]{options}} and is na.fail if that is 
#' not set.  The "factory-fresh" default is \code{\link[stats]{na.omit}}.
#' Another possible value is NULL, no action.  
#' Value \code{\link[stats]{na.exclude}} can be useful. 
#' @param ... further arguments
#' @return data.frame
#' @export
sample_MCF <- function(formula, data, subset, na.action, ...) {
  ## record the function call 
  Call <- match.call()
  ## arguments check
  if(missing(formula)) {
    stop("Argument 'formula' is missing.")
  } 
  if(missing(data)) {
    data <- environment(formula)
  }
  if (! with(data, inherits(eval(Call[[2]][[2]]), "Survr"))) {
    stop("'formula' must be a survival recurrent object.")
  }
  ## Prepare data: ID, time, event ~ X
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
  ## compute sample MCF 
  if (nbeta == 0) {
    num_pat <- length(unique(dat$ID))
    unitimes <- sort(unique(dat$time), decreasing = FALSE)
    unitimes <- head(unitimes, length(unitimes) - 1)
    num_at_risk <- apply(array(unitimes), 1, 
                         function(a) {
                           with(base::subset(dat, time <= a), 
                                num_pat - sum(event == 0))
                         })
    num_fails <- apply(array(unitimes), 1, 
                       function(a) {
                         with(base::subset(dat, time <= a), sum(event == 1))
                       })
    increment <- num_fails / num_at_risk
    sMCF <- cumsum(increment)
    ## data.frame to return
    num_at_risk <- c(length(unique(dat$ID)), num_at_risk)
    num_fails <- c(0, num_fails)
    increment <- c(0, increment)
    sMCF <- c(0, sMCF)
    unitimes <- c(0, unitimes)
    return(data.frame(time = unitimes, fails = num_fails, risk = num_at_risk, 
                      increment = increment, MCF = sMCF))
  } else {
    ## argument check
    if (nbeta != 1) {
      stop("The covariate in formula must be a factor or 1.")
    } 
    ## extract the covariate
    fac <- with(data, eval(Call[[2]][[3]]))
    if (! is.factor(fac)){
      stop("The covariate in formula must be a factor or 1.")    
    }
    ## number of levels
    num_levels <- length(levels(fac))
    outdat <- NULL
    for (i in seq(num_levels)) {
      subdat <- base::subset(dat, subset = fac %in% levels(fac)[i])
      num_pat <- length(unique(subdat$ID))
      unitimes <- sort(unique(subdat$time), decreasing = FALSE)
      unitimes <- head(unitimes, length(unitimes) - 1)
      num_at_risk <- apply(array(unitimes), 1, 
                           function(a) {
                             with(base::subset(subdat, time <= a), 
                                  num_pat - sum(event == 0))
                           })
      num_fails <- apply(array(unitimes), 1, 
                         function(a) {
                           with(base::subset(subdat, time <= a), sum(event == 1))
                         })
      increment <- num_fails / num_at_risk
      sMCF <- cumsum(increment)
      ## subdata.frame to return
      num_at_risk <- c(length(unique(subdat$ID)), num_at_risk)
      num_fails <- c(0, num_fails)
      increment <- c(0, increment)
      sMCF <- c(0, sMCF)
      unitimes <- c(0, unitimes)
      outdat <- rbind(outdat, data.frame(unitimes, num_fails, num_at_risk, 
                                         increment, sMCF, levels(fac)[i]))
    }
    colnames(outdat) <- c("time", "fails", "risk", 
                          "increment", "MCF", as.character(Call[[2]][[3]]))
    return(outdat)
  }
}

## Estimated Mean Cumulative Function (MCF) from HEART Model
heart_MCF <- function(object, design, CI = TRUE, level = 0.95) {
  beta <- object@estimates$beta[, 1]
  fcovnames <- as.character(object@call[[2]][[3]])
  covnames <- fcovnames[fcovnames != "+"]
  nbeta <- length(beta)
  if (missing(design)) {
    design <- rep(0, nbeta)
  } else if (length(design) != nbeta) {
    stop("The number of input covariates does not match with 'heart' object")
  } 
  coveff <- as.numeric(exp(crossprod(design, beta)))
  baselinepieces <- object@baselinepieces
  n_xx <- 1000
  n_pieces <- length(baselinepieces)
  BL_segments <- c(baselinepieces[1], diff(baselinepieces))
  xx <- seq(0, max(baselinepieces), length = n_xx)
  indx <- sapply(xx, whereT, baselinepieces = baselinepieces)
  LinCom_M <- t(apply(array(indx), 1, 
                      function(ind_indx) {
                        c(BL_segments[1:ind_indx], rep(0, n_pieces - ind_indx))
                      }))
  CMF_B4_indx <- c(0, baselinepieces)[indx]
  LinCom_M[(indx - 1) * n_xx + 1:n_xx] <- xx - CMF_B4_indx
  n_par <- dim(object@hessian)[1]
  Cov_M <- solve(object@hessian)[c((n_par - n_pieces + 1):n_par), 
                                 c((n_par - n_pieces + 1):n_par)]
  CI_band <- qnorm((1 + level)/2) * 
    sqrt(diag(LinCom_M %*% Cov_M %*% t(LinCom_M))) * coveff
  baseline_mean <- LinCom_M %*% object@estimates$alpha[, 1] * coveff
  lower <- baseline_mean - CI_band
  upper <- baseline_mean + CI_band
  ## return
  data.frame(time = xx, MCF = baseline_mean, lower = lower, upper = upper)
}
