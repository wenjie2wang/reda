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


#' Summarizing Model Fits
#'
#' \code{summary} returns summary of estimated coefficients of covariates,
#' rate function bases, and estimated parameter of frailty variable.
#'
#' Technitically, \code{summary} returns a
#' \code{\link{summaryHeart-class}} object,
#' which can be printed by
#' \code{\link{show,summaryHeart-method}}. 
#'
#' @param object rateReg object from \code{rateReg}.
#' @param showCall a logic value with dafault as TRUE,
#' indicating whether method \code{\link{show,summaryHeart-method}} prints out 
#' the call information of original call of \code{rateReg}.
#' @param showKnots a logic value with default as TRUE, 
#' indicating whether method \code{\link{show,summaryHeart-method}} prints out 
#' the internal and boundary knots.
#' @param ... other arguments for future usage.
#' @return summaryHeart-class object
#' @aliases summary,rateReg-method
#' @seealso \code{\link{rateReg}} \code{\link{coef,rateReg-method}}
#' \code{\link{confint,rateReg-method}} \code{\link{baseRate,rateReg-method}}
#' @importFrom methods new
#' @export
setMethod(f = "summary", signature = "rateReg",
          definition = function(object, showCall = TRUE,
                                showKnots = TRUE, ...) {
              Call <- object@call
              attr(Call, "show") <- showCall
              knots <- object@knots
              boundaryKnots <- object@boundaryKnots
              attr(knots, "show") <- showKnots
              beta <- object@estimates$beta
              theta <- object@estimates$theta
              alpha <- object@estimates$alpha
              results <- new("summaryHeart", 
                             call = Call,
                             knots = knots,
                             boundaryKnots = boundaryKnots,
                             covariateCoef = beta,
                             frailtyPar = theta,
                             degree = object@degree,
                             df = object@df, 
                             baseRateCoef = alpha,
                             logL = object@logL)
              ## return
              results
          })

