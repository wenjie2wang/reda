##
## R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
## Copyright (C) 2015-2025
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


##' Summarizing a Fitted Model
##'
##' Summary of estimated coefficients of covariates, rate function bases,
##' and estimated rate parameter of frailty random variable, etc.
##'
##' \code{summary,rateReg-method} returns a
##' \code{summary.rateReg} object,
##' whose slots include
##' \itemize{
##'     \item \code{covarCoef}: Estimated covariate coefficients.
##'     \item \code{frailtyPar}: Estimated rate parameter of gamma frailty.
##'     \item \code{baseRateCoef}: Estimated coeffcients of baseline
##'         rate function.
##' }
##' For the meaning of other slots, see \code{\link{rateReg}}.
##'
##' @param object A \code{rateReg} object.
##' @param showCall A logic value with dafault \code{TRUE},
##' indicating whether function \code{show}
##' prints out the original call information of \code{rateReg}.
##' It may be helpful for a more concise printout.
##' @param showKnots A logic value with default \code{TRUE},
##' indicating whether function \code{show}
##' prints out the internal and boundary knots.
##' Similar to argument \code{showCall}, It may be helpful
##' for a more concise printout.
##' @param ... Other arguments for future usage.
##'
##' @return \code{summary.rateReg} object
##'
##' @aliases summary,rateReg-method
##'
##' @examples
##' ## See examples given in function rateReg.
##'
##' @seealso
##' \code{\link{rateReg}} for model fitting;
##' \code{\link{coef,rateReg-method}} for point estimates of
##' covariate coefficients;
##' \code{\link{confint,rateReg-method}} for confidence intervals
##' of covariate coeffcients;
##' \code{\link{baseRate,rateReg-method}} for coefficients of baseline
##' rate function.
##' @export
setMethod(f = "summary", signature = "rateReg",
          definition = function(object, showCall = TRUE,
                                showKnots = TRUE, ...) {
              Call <- object@call
              attr(Call, "show") <- showCall
              knots <- object@spline$knots
              Boundary.knots <- object@spline$Boundary.knots
              attr(knots, "show") <- showKnots
              beta <- object@estimates$beta
              theta <- object@estimates$theta
              alpha <- object@estimates$alpha
              ## check on object validity by 'new', validObject(results)
              results <- new("summary.rateReg",
                             call = Call,
                             spline = object@spline$spline,
                             knots = knots,
                             Boundary.knots = Boundary.knots,
                             covarCoef = beta,
                             frailtyPar = theta,
                             degree = object@spline$degree,
                             baseRateCoef = alpha,
                             logL = object@logL)
              ## return
              results
          })


##' Summarize an \code{Recur} object
##'
##' @param object An \code{Recur} object.
##' @param ... Other arguments not used.
##'
##' @return \code{summary.Recur} object.
##'
##' @importFrom stats median
##'
##' @export
setMethod(f = "summary", signature = "Recur",
          definition = function(object, ...) {
              Call <- object@call
              n <- length(object@first_idx)
              d0 <- object@.Data[object@.Data[, "event"] == 0, ]
              y <- d0[, "time2"] - d0[, "origin"]
              d <- d0[, "terminal"]
              oy <- order(y)
              d <- d[oy]
              y <- y[oy]
              r <- n - rank(y, ties.method = "min") + 1
              is_first_y <- ! duplicated(y)
              s <- cumprod(1 - (d / r)[is_first_y])
              medTem <- if (s[length(s)] > .5) {
                            NA_real_
                        } else {
                            as.numeric(y[is_first_y][which.max(s - .5 < 0)])
                        }
              new("summary.Recur",
                  call = Call,
                  sampleSize = as.integer(n),
                  reSize = as.integer(sum(object@.Data[, "event"] > 0)),
                  avgReSize = sum(object@.Data[, "event"]) / n,
                  propTem = sum(object@.Data[, "terminal"]) / n,
                  medFU = median(y),
                  medTem = medTem)
          })
