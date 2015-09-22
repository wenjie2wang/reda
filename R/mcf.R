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
#' For \code{\link{heart-class}} object, 
#' \code{mcf} estimates the MCF of baseline rate function.
#' 
#' @param object an object used to dispatch a method.
#' @param ... other arguments for future usage.
#' @param level a optional numeric value 
#' indicating the confidence level required. 
#' For \code{mcf,formula-method}, it must be between 0.5 and 1.
#' The default value is 0.95.
#' @seealso \code{\link{heart}} \code{\link{plotMCF}}
#' @examples 
#' library(reda)
#' data(simuDat)
#' 
#' ## empirical MCF
#' sampleMCF <- mcf(Survr(ID, time, event) ~ group, data = simuDat)
#' 
#' mcf(Survr(ID, time, event) ~ group, data = simuDat, 
#' subset = ID %in% 100:101, na.action = na.omit)
#' 
#' ## estimated MCF for baseline rate function from HEART model
#' heartfit <- heart(formula = Survr(ID, time, event) ~ X1 + group, 
#'                   data = simuDat, subset = ID %in% 75:125,
#'                   baselinepieces = seq(28, 168, length = 6))
#' baselineMCF <- mcf(heartfit)
#' 
#' mcf(heartfit, level = 0.9, control = list(length.out = 500))
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
#' @return \code{\link{empirMCF-class}} or \code{\link{heartMCF-class}} object
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
              stop("Response in formula must be a survival recurrent object.")
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
            nbeta <- ncol(mm) - 1 
            nsample <- nrow(mm)
            if (nbeta == 0) { 
              X <- covar_names <-NULL
            } else {
              X <- mm[, 2]
              ## covariates' names
              covar_names <- attr(Terms, "term.labels")
            }
            ## data 
            dat <- as.data.frame(cbind(mm[, 1], X))
            colnames(dat) <- c("ID", "time", "event", covar_names)
            
            ## compute sample MCF 
            ### function part
            sMCF <- function(inpdat) {
              outdat <- plyr::ddply(inpdat, "time", function(a){
                a[order(a$event, decreasing = TRUE), c("ID", "time", "event")]
              })
              nsample <- nrow(outdat)
              num_pat <- length(unique(outdat$ID))
              num_at_risk <- apply(array(seq(nsample)), 1, 
                                   function(a) {
                                     with(outdat[seq(a), ], 
                                          num_pat - sum(event == 0))
                                   })
              increment <- ifelse(outdat$event == 1, 1/num_at_risk, 0)
              smcf <- cumsum(increment)
              incre_var <- ifelse(outdat$event == 0, 0, 
                                  increment^2 * 
                                    ((1 - increment)^2 + 
                                       (num_at_risk - 1) * increment^2))
              var_smcf <- cumsum(incre_var)
              wtonexp <- qnorm(level) * sqrt(var_smcf) / smcf
              upper <- smcf * exp(wtonexp)
              lower <- smcf * exp(-wtonexp)
              # return
              data.frame(outdat, MCF = smcf, var = var_smcf, lower, upper)
            }
            
            if (nbeta == 0) {
              outdat <- sMCF(dat)
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
                resdat <- data.frame(sMCF(subdat), levels(X)[i])
                ## subdata.frame to return
                outdat <- rbind(outdat, resdat)
              }
              colnames(outdat)[ncol(outdat)] <- covar_names
              out <- methods::new("empirMCF", call = Call, formula = object, 
                                  MCF = outdat, multigroup = TRUE)
              return(out)
            }
          })


#' @describeIn mcf Estimated Mean Cumulative Function (MCF) from HEART Model
#' 
#' @param newdata an optional data.frame. 
#' @param groupname an optional length-one charactor vector. 
#' @param grouplevels an optional charactor vector.
#' @param control an optional list to specify the time grid 
#' where the MCF are estimated.
#' The possible elements of the control list include 
#' \code{grid}, \code{length.out}, \code{from} and \code{to}.
#' The time grid can be directly specified via element \code{grid}.
#' \code{length.out} represents the length of grid points.
#' The dafault value is 200.
#' \code{from} means the starting point of grid. It takes 0 as default.
#' \code{to} means the endpoint of grid. 
#' It takes the endpoint of baseline pieces as default.
#' When \code{grid} is missing, the grid will be generated 
#' via \code{\link{seq}} with arguments \code{from}, \code{to} 
#' and \code{length.out}
#' @aliases mcf,heart-method
#' @importFrom methods new
#' @importFrom stats terms na.fail na.omit na.exclude na.pass qnorm model.matrix
#' model.frame delete.response 
#' @export
setMethod(f = "mcf", signature = "heart", 
          definition = function(object, newdata, groupname, grouplevels, 
                                level = 0.95, na.action, control = list(), 
                                ...) {
            beta <- object@estimates$beta[, 1]
            alpha <- object@estimates$alpha[, 1]
            fcovnames <- as.character(object@call[[2]][[3]])
            covnames <- fcovnames[fcovnames != "+"]
            nbeta <- length(beta)
            baselinepieces <- object@baselinepieces
            controlist <- c(control, 
                            list("baselinepieces" = as.numeric(baselinepieces)))
            control <- do.call("heart_mcf_control", controlist)
            n_xx <- control$length.out
            n_pieces <- length(baselinepieces)
            BL_segments <- c(baselinepieces[1], diff(baselinepieces))
            xx <- control$grid
            indx <- sapply(xx, whereT, baselinepieces = baselinepieces)
            LinCom_M <- NULL
            for (ind_indx in indx) {
                 LinCom_M <- rbind(LinCom_M, c(BL_segments[1:ind_indx], 
                                               rep(0, n_pieces - ind_indx)))
            }
            CMF_B4_indx <- c(0, baselinepieces)[indx]
            LinCom_M[(indx - 1) * n_xx + 1:n_xx] <- xx - CMF_B4_indx
            n_par <- nrow(object@fisher)
            ## covariance matrix of beta and alpha
            Cov_ind <- c(seq(nbeta), (n_par - n_pieces + 1):n_par)
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
            coveff <- as.numeric(exp(crossprod(t(X), beta)))
            outdat <- NULL
            for (i in seq(ndesign)) {
              ## Delta-method
              grad <- cbind(alpha %o% X[i, ], diag(rep(1, n_pieces))) * 
                coveff[i]
              Cov_M <- tcrossprod(crossprod(t(grad), Cov_par), grad)
              CI_band <- qnorm((1 + level)/2) * 
                sqrt(diag(tcrossprod(crossprod(t(LinCom_M), Cov_M), LinCom_M))) 
              baseline_mean <- crossprod(t(LinCom_M), alpha) * coveff[i]
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
heart_mcf_control <- function (grid, length.out = 200, from, to, 
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
    length.out <- length(grid)
    from <- min(grid)
    to <- max(grid)  
  } else {
    grid <- seq(from = from, to = to, length.out = length.out)
  }
  if (min(grid) < 0 || max(grid) > max(baselinepieces)) {
    stop("'grid' must be within the coverage of baseline pieces.")
  }
  ## return
  list(grid = grid, length.out = length.out, from = from, to = to)
}

