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


#' Mean Cumulative Function (MCF)
#' 
#' S4 class generic function to compute mean empirical cumulative function (MCF)
#' from sample data or estimated MCF from HEART model.
#' 
#' For formula with \code{\link{Survr}} object as response, 
#' the covariate specified in the rhs of the formula can either be 1 or 
#' any one factor variable in the data.  The former computes the overall 
#' empirical MCF from sample.  The latter computes the empirical MCF for each 
#' level of the factor variable specified respectively.
#' \code{mcf} computes empirical mean cumulative function (MCF) on every unique
#' time point from recurrent event sample data.  
#' It does not assume any particular underlying model. 
#' 
#' For \code{\link{rateReg-class}} object, 
#' \code{mcf} estimates the MCF of baseline rate function.
#' 
#' @param object an object used to dispatch a method.
#' @param ... other arguments for future usage.
#' @param level a optional numeric value 
#' indicating the confidence level required. 
#' For \code{mcf,formula-method}, it must be between 0.5 and 1.
#' The default value is 0.95.
#' @seealso \code{\link{rateReg}} \code{\link{plotMcf}}
#' @examples 
#' library(reda)
#'  
#' ## empirical MCF
#' sampleMCF <- mcf(Survr(ID, time, event) ~ group, data = simuDat)
#' 
#' mcf(Survr(ID, time, event) ~ group, data = simuDat, 
#' subset = ID %in% 100:101, na.action = na.omit)
#' 
#' ## estimated MCF for baseline rate function from HEART model
#' rateRegFit <- rateReg(formula = Survr(ID, time, event) ~ x1 + group, 
#'                   data = simuDat, subset = ID %in% 75:125,
#'                   bKnots = seq(28, 168, length = 6))
#' baselineMCF <- mcf(rateRegFit)
#' 
#' mcf(rateRegFit, level = 0.9, control = list(length.out = 500))
#'  
#' @export
setGeneric(name = "mcf",
           def = function(object, ...) {
               standardGeneric("mcf")
           })


#' @describeIn mcf Empirical mean cumulative function (MCF)
#' 
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
#' @return \code{\link{empirMcf-class}} or \code{\link{rateRegMcf-class}} object
#' @aliases mcf,formula-method 
#' @importFrom utils head 
#' @importFrom methods new
#' @importFrom stats na.fail na.omit na.exclude na.pass
#' @importFrom plyr ddply
#' @export
setMethod(f = "mcf", signature = "formula", 
          definition = function(object, data, subset, na.action, 
                                level = 0.95, ...) {
              ## record the function call 
              Call <- match.call()
              Call[[1]] <- as.name("mcf")
              names(Call) <- sub(pattern = "object", 
                                 replacement = "formula", names(Call))
              mfnames <- c("formula", "data", "subset", "na.action")
              mfind <- match(mfnames, names(Call), nomatch = 0L)
              Call$formula <- eval(object)
              Call$data <- eval(substitute(alist(data)))[[1]]
              Call$subset <- eval(substitute(alist(subset)))[[1]]
              Call$na.action <- eval(substitute(alist(na.action)))[[1]]
              mcall <- Call[c(1, mfind)]
              ## Call to return
              Call <- mcall
              
              ## drop unused levels in factors 
              mcall$drop.unused.levels <- TRUE
              ## Prepare data: ID, time, event ~ X
              mcall[[1]] <- quote(stats::model.frame)
              if (is.R()) {
                  mm <- eval.parent(mcall)
              } else { 
                  mm <- eval(mcall, sys.parent()) 
              }
              if(missing(data)) {
                  data <- environment(object)
              }
              if (! with(data, inherits(eval(object[[2]]), "Survr"))) {
                  stop("Response in formula must be a 'Survr' object.")
              }
              Terms <- terms(object)
              ord <- attr(Terms, "order")
              if (length(ord) & any(ord != 1)) {
                  stop("Interaction terms are not valid for this function")
              }
              ## check for confidence level: 0.5 < level < 1 
              if (level <= 0.5 | level >= 1) {
                  stop("Confidence level specified must be between 0.5 and 1.")
              }
              ## number of covariates excluding intercept
              nBeta <- ncol(mm) - 1L
              nSample <- nrow(mm)
              if (nBeta == 0L) { 
                  X <- covar_names <- NULL
              } else {
                  X <- mm[, 2]
                  ## covariates' names
                  covar_names <- attr(Terms, "term.labels")
              }
              ## data 
              dat <- as.data.frame(cbind(mm[, 1], X))
              colnames(dat) <- c("ID", "time", "event", covar_names)     
              
              if (nBeta == 0L) { ## if no covariates specified
                  outDat <- sMcf(dat, level = level)
                  out <- new("empirMcf", formula = object, 
                             MCF = outDat, multiGroup = FALSE)
                  return(out)
              }
              ## else one covariate specified
              ## argument check
              if (nBeta != 1L && ! is.factor(X)) {
                  stop("The covariate in object must be a factor or '1'.")
              } 
              ## number of levels
              num_levels <- length(levels(X))
              if (num_levels == 1) {
                  warning("The factor covariate has only one level.")
              }
              # browser()
              outDat <- data.frame(matrix(NA, nrow = nSample, ncol = 8))
              for (i in seq(num_levels)) {
                  subdat <- base::subset(dat, subset = X %in% levels(X)[i])
                  rowInd <- if(i == 1L) {
                      seq(nrow(subdat))
                  } else {
                      seq(from = tail(rowInd, 1) + 1, by = 1,
                          length.out = nrow(subdat))                            
                  }
                  outDat[rowInd, ] <- data.frame(sMcf(subdat, level = level),
                                                 levels(X)[i])
              }
              colnames(outDat) <- c(colnames(dat)[1:3], "MCF",
                                    "var", "lower", "upper", covar_names)
              out <- methods::new("empirMcf", call = Call, formula = object, 
                                  MCF = outDat, multiGroup = TRUE)
              ## return
              out
          })


#' @describeIn mcf Estimated Mean Cumulative Function (MCF) from Model Fits
#' 
#' @param newdata An optional data.frame. 
#' @param groupName An optional length-one charactor vector. 
#' @param groupLevels An optional charactor vector.
#' @param control An optional list to specify the time grid 
#' where the MCF are estimated.
#' The possible elements of the control list include 
#' \code{grid}, \code{length.out}, \code{from} and \code{to}.
#' The time grid can be directly specified via element \code{grid}.
#' \code{length.out} represents the length of grid points.
#' The dafault value is 200.
#' \code{from} mean the starting point of grid. It takes 0 as default.
#' \code{to} means the endpoint of grid. 
#' It takes the endpoint of baseline pieces as default.
#' When \code{grid} is missing, the grid will be generated 
#' via \code{\link{seq}} with arguments \code{from}, \code{to} 
#' and \code{length.out}
#' @aliases mcf,rateReg-method
#' @importFrom methods new
#' @importFrom stats terms na.fail na.omit na.exclude na.pass qnorm model.matrix
#' model.frame delete.response 
#' @export
setMethod(f = "mcf", signature = "rateReg", 
          definition = function(object, newdata, groupName, groupLevels, 
                                level = 0.95, na.action, control = list(), 
                                ...) {
              beta <- object@estimates$beta[, 1]
              alpha <- object@estimates$alpha[, 1]
              fcovnames <- as.character(object@call[[2]][[3]])
              covnames <- fcovnames[fcovnames != "+"]
              nBeta <- length(beta)
              knots <- object@knots
              degree <- object@degree
              boundaryKnots <- object@boundaryKnots
              bKnots <- c(knots, boundaryKnots[2])
              controlist <- c(control, list("bKnots" = bKnots)
              control <- do.call("rateReg_mcf_control", controlist)
              n_xx <- control$length.out

              
              ## compute mcf
              mu0(par_BaselinePW = alpha, Tvec = grid, bKnots = bKnots,
                  degree = degree, bsMat_est)

              if (degree == 0L) { ## if piecewise constant
                  n_pieces <- length(bKnots)
                  BL_segments <- c(bKnots[1], diff(bKnots))
                  xx <- control$grid
                  indx <- sapply(xx, whereT, bKnots = bKnots)
                  LinCom_M <- matrix(NA, ncol = n_pieces, nrow = n_xx)
                  i <- 1
                  for (ind_indx in indx) {
                      LinCom_M[i] <- c(BL_segments[1:ind_indx], 
                                       rep(0, n_pieces - ind_indx))
                      i <- i + 1
                  }
                  CMF_B4_indx <- c(boundaryKnots[1], bKnots)[indx]
                  LinCom_M[(indx - 1) * n_xx + seq(n_xx)] <- xx - CMF_B4_indx
              } else { ## spline rate function of degree at least one
                  
              }
              
              
              n_par <- nrow(object@fisher)
              ## covariance matrix of beta and alpha
              Cov_ind <- c(seq(nBeta), (n_par - n_pieces + 1):n_par)
              Cov_par <- solve(object@fisher)[Cov_ind, Cov_ind]

              ## nonsense, just to suppress Note from R CMD check --as-cran
              `(Intercept)` <- NULL
              
              ## about newdata
              tt <- terms(object@formula)
              Terms <- delete.response(tt)
              if (missing("na.action")) {
                  na.action <- options("na.action")[[1]]
              }
              if (missing(newdata)) {
                  X <- matrix(rep(0, nBeta), nrow = 1)
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
                  if (ncol(X) != nBeta) {
                      stop("The number of input covariates does not match 
                   with 'rateReg' object")
                  }
              }
              ndesign <- nrow(X)
              multiGroup <- ifelse(ndesign == 1, FALSE, TRUE)
              coveff <- as.numeric(exp(crossprod(t(X), beta)))
              outDat <- matrix(NA, ncol = 4, nrow = ndesign * length(xx))
              for (i in seq(ndesign)) {
                  ## Delta-method
                  grad <- cbind(alpha %o% X[i, ], diag(rep(1, n_pieces))) * 
                      coveff[i]
                  Cov_M <- tcrossprod(crossprod(t(grad), Cov_par), grad)
                  criti <- diag(tcrossprod(crossprod(t(LinCom_M), Cov_M),
                                           LinCom_M))
                  CI_band <- qnorm((1 + level) / 2) * sqrt(criti)
                  baseline_mean <- crossprod(t(LinCom_M), alpha) * coveff[i]
                  lower <- baseline_mean - CI_band
                  upper <- baseline_mean + CI_band
                  ind <- seq(from = length(xx) * (i - 1) + 1,
                             to = length(xx) * i, by = 1)
                  outDat[ind, 1] <- xx
                  outDat[ind, 2] <- baseline_mean
                  outDat[ind, 3] <- lower
                  outDat[ind, 4] <- upper
              }
              outDat <- as.data.frame(outDat)
              colnames(outDat) <- c("time", "MCF", "lower", "upper")
              if (multiGroup) {
                  if (missing(groupLevels)) {
                      groupLevels <- LETTERS[seq(ndesign)]
                  }
                  if (missing(groupName)) {
                      groupName <- "group"
                  }
                  tempcol <- factor(rep(groupLevels[seq(ndesign)], each = n_xx))
                  outDat <- cbind(outDat, tempcol)
                  colnames(outDat) <- c("time", "MCF", "lower", "upper",
                                        groupName)
              }
              ## output
              out <- new("rateRegMcf", 
                         formula = object@formula, 
                         bKnots = object@bKnots, 
                         newdata = X, MCF = outDat, level = level,
                         na.action = na.action, control = control, 
                         multiGroup = multiGroup)
              ## return
              out
          })


## internal function ===========================================================
rateReg_mcf_control <- function (grid, length.out = 200, from, to, 
                                 bKnots) {
    ## controls for function MCF with signiture rateReg
    if (missing(from)) {
        from <- 0
    }
    if (missing(to)) {
        to <- max(bKnots)
    }
    if (! missing(grid)) {
        if (! is.numeric(grid) || is.unsorted(grid)) {
            stop("'grid' specified must be a increasing numeric vector.")
        }
        length.out <- length(grid)
        from <- min(grid)
        to <- max(grid)  
    } else {
        grid <- seq(from = from, to = to, length.out = length.out)
    }
    if (min(grid) < 0 || max(grid) > max(bKnots)) {
        stop("'grid' must be within the coverage of baseline pieces.")
    }
    ## return
    list(grid = grid, length.out = length.out, from = from, to = to)
}

## compute sample MCF 
sMcf <- function(inpDat, level) {
    outDat <- plyr::ddply(inpDat, "time", function(a){
        a[order(a$event, decreasing = TRUE),
          c("ID", "time", "event")]
    })
    nSample <- nrow(outDat)
    num_pat <- length(unique(outDat$ID))
    num_at_risk <- apply(array(seq(nSample)), 1, 
                         function(a) {
                             with(outDat[seq(a), ], 
                                  num_pat - sum(event == 0))
                         })
    increment <- ifelse(outDat$event == 1, 1 / num_at_risk, 0)
    smcf <- cumsum(increment)
    inc2 <- increment ^ 2
    inc12 <- (1 - increment) ^ 2
    incre <- inc2 * (inc12 + (num_at_risk - 1) * inc2)
    incre_var <- ifelse(outDat$event == 0, 0, incre)
    var_smcf <- cumsum(incre_var)
    wtonexp <- qnorm(level) * sqrt(var_smcf) / smcf
    upper <- smcf * exp(wtonexp)
    lower <- smcf * exp(- wtonexp)
    ## return
    data.frame(outDat, MCF = smcf, var = var_smcf, lower, upper)
}
