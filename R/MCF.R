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


#' Mean Cumulative Function (MCF)
#' 
#' S4 class generic function to compute mean empirical cumulative function (MCF)
#' from sample data or estimated MCF from HEART model.
#' 
#' @usage MCF(object, ...)
#' @param object an object used to dispatch a method.
#' @param ... further arguments passed to or from other methods.
#' @export
setGeneric(name = "MCF",
           def = function(object, ...) {
             standardGeneric("MCF")
           })


#' An S4 class to represent computed empirical MCF
#' @slot formula formula 
#' @slot MCF data.frame
#' @slot multigroup logical value
#' @importFrom methods setClass
#' @export
setClass(Class = "empirMCF", 
         slots = c(formula = "formula", MCF = "data.frame", 
                   multigroup = "logical"))


#' @describeIn MCF Empirical mean cumulative function (MCF)
#' 
#' \code{MCF} computes empirical mean cumulative function (MCF) on every
#' time point from recurrent event sample data.  
#' It does not assume any particular underlying model. 
#' @param object Survr object from function \code{Survr} in package survrec.
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
#' @param ... further arguments.
#' @return empirMCF object
#' @importFrom utils head 
#' @importFrom methods new
#' @importFrom stats na.fail na.omit na.exclude na.pass
#' @export
setMethod(f = "MCF", signature = "formula", 
          definition = function(object, data, subset, na.action, ...) {
            ## record the function call 
            Call <- match.call()
            Call[[1]] <- as.name("MCF")
            names(Call) <- sub(pattern = "object", 
                               replacement = "formula", names(Call))
            mfnames <- c("formula", "data", "subset", "na.action")
            mcall <- Call[c(1, match(mfnames, names(Call), nomatch = 0L))]
            ## drop unused levels in factors 
            mcall$drop.unused.levels <- TRUE
            ## Prepare data: ID, time, event ~ X
            mcall[[1]] <- quote(stats::model.frame)
            if (is.R()) {
              mm <- eval.parent(mcall)
            } else { 
              mm <- eval(mcall, sys.parent()) 
            }
            if (! inherits(mm[, 1], "Survr")) {
              stop("Response must be a survival recurrent object.")
            }
            Terms <- terms(object)
            ord <- attr(Terms, "order")
            if (length(ord) & any(ord != 1)) {
              stop("Interaction terms are not valid for this function")
            }
            ## number of covariates excluding intercept
            nbeta <- ncol(mm) - 1 
            nsample <- nrow(mm)
            X <- if (nbeta == 0) { 
              factor(rep(1, nsample)) 
            } else {
              mm[, 2]
            }
            ## covariates' names
            covar_names <- attr(Terms, "term.labels")
            ## data 
            dat <- as.data.frame(cbind(mm[, 1], X))
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
                                   with(base::subset(dat, time <= a), 
                                        sum(event == 1))
                                 })
              increment <- num_fails / num_at_risk
              sMCF <- cumsum(increment)
              ## data.frame to return
              num_at_risk <- c(length(unique(dat$ID)), num_at_risk)
              num_fails <- c(0, num_fails)
              increment <- c(0, increment)
              sMCF <- c(0, sMCF)
              unitimes <- c(0, unitimes)
              outdat <- data.frame(time = unitimes, fails = num_fails, 
                                   risk = num_at_risk, increment = increment, 
                                   MCF = sMCF)
              out <- new("empirMCF", formula = object, 
                         MCF = outdat, multigroup = FALSE)
              return(out)
            } else {
              ## argument check
              if (nbeta != 1 && ! is.factor(X)) {
                stop("The covariate in object must be a factor or 1.")
              } 
              ## number of levels
              num_levels <- length(levels(X))
              if (num_levels == 1) {
                warning("The factor covariate has only one level.")
              }
              outdat <- NULL
              for (i in seq(num_levels)) {
                subdat <- base::subset(dat, subset = X %in% levels(X)[i])
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
                                     with(base::subset(subdat, time <= a), 
                                          sum(event == 1))
                                   })
                increment <- num_fails / num_at_risk
                sMCF <- cumsum(increment)
                ## subdata.frame to return
                num_at_risk <- c(length(unique(subdat$ID)), num_at_risk)
                num_fails <- c(0, num_fails)
                increment <- c(0, increment)
                sMCF <- c(0, sMCF)
                unitimes <- c(0, unitimes)
                outdat <- rbind(outdat, data.frame(unitimes, num_fails, 
                                                   num_at_risk, increment, 
                                                   sMCF, levels(X)[i]))
              }
              colnames(outdat) <- c("time", "fails", "risk", 
                                    "increment", "MCF", covar_names)
              out <- methods::new("empirMCF", formula = object, 
                                  MCF = outdat, multigroup = TRUE)
              return(out)
            }
          })


#' An S4 class to represent estimated MCF from HEART model
#' @slot formula formula
#' @slot baselinepieces numeric vector.
#' @slot newdata numeric matrix.
#' @slot MCF data.frame.
#' @slot level a numeric value within (0, 1).
#' @slot na.action a length-one character vector.
#' @slot control list.
#' @slot multigroup logical. 
#' @export
setClass(Class = "heartMCF", 
         slots = c(formula = "formula", baselinepieces = "numeric",
                   newdata = "matrix", MCF = "data.frame", level = "numeric", 
                   na.action = "character", control = "list", 
                   multigroup = "logical"))


#' @describeIn MCF Estimated Mean Cumulative Function (MCF) from HEART Model
#' 
#' @param object heart object.
#' @param newdata data.frame.
#' @param groupname a length-one charactor vector.
#' @param grouplevels a charactor vector.
#' @param level a length-one numeric vector.
#' @param na.action function.
#' @param control list.
#' @param ... further arguments.
#' @importFrom methods new
#' @importFrom stats terms na.fail na.omit na.exclude na.pass qnorm model.matrix
#' model.frame delete.response 
#' @export
setMethod(f = "MCF", signature = "heart", 
          definition = function(object, newdata, groupname, grouplevels, 
                                level = 0.95, na.action, control = list(), 
                                ...) {
            beta <- object@estimates$beta[, 1]
            fcovnames <- as.character(object@call[[2]][[3]])
            covnames <- fcovnames[fcovnames != "+"]
            nbeta <- length(beta)
            baselinepieces <- object@baselinepieces
            controlist <- c(control, 
                            list("baselinepieces" = as.numeric(baselinepieces)))
            control <- do.call("heart_MCF_control", controlist)
            n_xx <- control$numgrid
            n_pieces <- length(baselinepieces)
            BL_segments <- c(baselinepieces[1], diff(baselinepieces))
            xx <- control$grid
            indx <- sapply(xx, whereT, baselinepieces = baselinepieces)
            LinCom_M <- t(apply(array(indx), 1, 
                                function(ind_indx) {
                                  c(BL_segments[1:ind_indx], 
                                    rep(0, n_pieces - ind_indx))
                                }))
            CMF_B4_indx <- c(0, baselinepieces)[indx]
            LinCom_M[(indx - 1) * n_xx + 1:n_xx] <- xx - CMF_B4_indx
            n_par <- nrow(object@hessian)
            Cov_M <- solve(object@hessian)[c((n_par - n_pieces + 1):n_par), 
                                           c((n_par - n_pieces + 1):n_par)]
            ## about newdata
            tt <- terms(object@formula)
            Terms <- delete.response(tt)
            if (missing("na.action")) {
              na.action <- options("na.action")[[1]]
            }
            if (missing(newdata)) {
              X <- matrix(rep(0, nbeta), nrow = 1)
              colnames(X) <- covnames
              rownames(X) <- "1"
            } else {
              mf <- model.frame(Terms, newdata, na.action = na.action, 
                                xlev = object@xlevels)
              if (is.null(attr(mf, "na.action"))) {
                na.action <- options("na.action")[[1]]
              } else {
                na.action <- paste("na", class(attr(mf, "na.action")), 
                                   sep = ".")
              }
              X <- model.matrix(Terms, mf, contrasts.arg = object@contrasts)
              ## remove intercept and deplicated rows
              X <- unique(base::subset(X, select = -`(Intercept)`))
              if (ncol(X) != nbeta) {
                stop("The number of input covariates does not match 
                   with 'heart' object")
              }
            }
            ndesign <- nrow(X)
            multigroup <- ifelse(ndesign == 1, FALSE, TRUE)
            coveff <- as.numeric(exp(crossprod(X, beta)))
            outdat <- NULL
            for (i in seq(ndesign)) {
              CI_band <- qnorm((1 + level)/2) * 
                sqrt(diag(LinCom_M %*% Cov_M %*% t(LinCom_M))) * coveff[i]
              baseline_mean <- LinCom_M %*% 
                object@estimates$alpha[, 1] * coveff[i]
              lower <- baseline_mean - CI_band
              upper <- baseline_mean + CI_band
              outdat <- rbind(outdat, data.frame(time = xx, MCF = baseline_mean, 
                                                 lower = lower, upper = upper))
            }
            if (multigroup) {
              if (missing(grouplevels)) {
                grouplevels <- LETTERS[seq(ndesign)]
              }
              if (missing(groupname)) {
                groupname <- "group"
              }
              tempcol <- factor(rep(grouplevels[seq(ndesign)], each = n_xx))
              outdat <- cbind(outdat, tempcol)
              colnames(outdat)[ncol(outdat)] <- groupname
            }
            ## output
            out <- new("heartMCF", 
                       formula = object@formula, 
                       baselinepieces = object@baselinepieces, 
                       newdata = X, MCF = outdat, level = level,
                       na.action = na.action, control = control, 
                       multigroup = multigroup)
            ## return
            out
          })
          

## internal function ===========================================================
heart_MCF_control <- function (grid, numgrid = 1000, from, to, 
                               baselinepieces) {
  ## controls for function MCF with signiture heart
  if (missing(from)) {
    from <- 0
  }
  if (missing(to)) {
    to <- max(baselinepieces)
  }
  if (! missing(grid)) {
    if (! is.numeric(grid) || is.unsorted(grid)) {
      stop("'grid' specified must be a increasing numeric vector.")
    }
    numgrid <- length(grid)
    from <- min(grid)
    to <- max(grid)  
  } else {
    grid <- seq(from = from, to = to, length.out = numgrid)
  }
  if (min(grid) < 0 || max(grid) > max(baselinepieces)) {
    stop("'grid' must be within the coverage of baseline pieces.")
  }
  ## return
  list(grid = grid, numgrid = numgrid, from = from, to = to)
}
