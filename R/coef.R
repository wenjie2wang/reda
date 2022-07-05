##
## R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
## Copyright (C) 2015-2022
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


##' Estimated Coefficients of Covariates
##'
##' \code{coef,rateReg-method} is an S4 class method that extracts estimated
##' coefficients of covariates from \code{rateReg} object produced by function
##' \code{\link{rateReg}}.
##'
##' @param object A \code{rateReg} object.
##' @param ... Other arguments for future usage.
##'
##' @return A named numeric vector.
##'
##' @aliases coef,rateReg-method
##'
##' @seealso
##' \code{\link{rateReg}} for model fitting;
##' \code{\link{confint,rateReg-method}} for confidence intervals
##' for covariate coefficients;
##' \code{\link{summary,rateReg-method}} for summary of a fitted model.
##'
##' @examples
##' ## See examples given in function rateReg.
##' @importFrom stats coef
##' @export
setMethod(f = "coef", signature = "rateReg",
          definition = function(object, ...) {
              object@estimates$beta[, "coef"]
          })


##' Confidence Intervals for Covariate Coefficients
##'
##' \code{confint,rateReg-method} is an S4 class method for
##' \code{\link{rateReg}} object, which returns approximate confidence intervals
##' for all or specified covariates.
##'
##' Under regularity condition (Shao 2003, Theorem 4.16 and Theorem 4.17, page
##' 287, 290), the approximate confidence intervals are constructed loosely
##' based on Fisher information matrix and estimates of coefficients.
##'
##' @param object A \code{rateReg} object.
##' @param parm A specification of which parameters are to be given confidence
##'     intervals, either a vector of numbers or a vector of names.  If missing,
##'     all parameters are considered.
##' @param level An optional numeric value to specify the confidence level
##'     required.  By default, the value is 0.95, which produces 95\% confidence
##'     intervals.
##' @param ... Other arguments for future usage.
##'
##' @return A numeric matrix with row names and column names.
##'
##' @aliases confint,rateReg-method
##'
##' @seealso
##' \code{\link{rateReg}} for model fitting;
##' \code{\link{coef,rateReg-method}} for point estimates
##' of covariate coefficients;
##' \code{\link{summary,rateReg-method}} for summary of a fitted model.
##'
##' @references
##'
##' Shao, J. (2003), \emph{Mathematical statistics},
##' Springer texts in statistics, New York: Springer, 2nd Edition.
##'
##' @examples
##' ## See examples given in function rateReg.
##' @importFrom stats confint qnorm
##' @export
setMethod(f = "confint", signature = "rateReg",
          definition = function(object, parm, level = 0.95, ...) {
              ## internal function
              format.perc <- function(probs) {
                  paste0(format(100 * probs, trim = TRUE,
                                scientific = FALSE), "%")
              }
              betaMat <- object@estimates$beta
              estCoef <- betaMat[, "coef"]
              pnames <- attr(betaMat, "dimnames")[[1L]]
              if (missing(parm)) {
                  parm <- seq_len(nrow(betaMat))
              } else if (is.numeric(parm)) {
                  parm <- intersect(seq_len(nrow(betaMat)), parm)
              } else if (is.character(parm)) {
                  parm <- match(parm, pnames, nomatch = NULL)
              } else {
                  stop("invalid argument 'parm'.")
              }
              a <- (1 + c(- 1, 1) * level) / 2
              fac <- stats::qnorm(a)
              pct <- format.perc(a)
              ci <- array(NA, dim = c(length(parm), 2L),
                          dimnames = list(parm, pct))
              ses <- betaMat[parm, "se(coef)"]
              ci[] <- estCoef[parm] + ses %o% fac
              rownames(ci) <- pnames[parm]
              ## return
              ci
          })
