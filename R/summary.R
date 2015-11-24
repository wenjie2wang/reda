################################################################################
##
##   R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
##   Copyright (C) 2015
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
#' @include class.R 
NULL


#' Summarizing a Fitted Model
#'
#' \code{summary} mainly returns summary of estimated coefficients of
#' covariates, rate function bases, and estimated parameter of frailty variable.
#'
#' Technitically, \code{summary} returns a
#' \code{\link{summaryRateReg-class}} object,
#' which can be printed by
#' \code{\link{show,summaryRateReg-method}}. 
#'
#' @param object rateReg object from \code{rateReg}.
#' @param showCall A logic value with dafault as TRUE,
#' indicating whether method \code{\link{show,summaryRateReg-method}}
#' prints out the call information of original call of \code{rateReg}.
#' @param showKnots A logic value with default as TRUE, 
#' indicating whether method \code{\link{show,summaryRateReg-method}}
#' prints out the internal and boundary knots.
#' @param ... Other arguments for future usage.
#' @return summaryRateReg-class object
#' @aliases summary,rateReg-method
#' @examples
#' ## See examples given in \code{\link{rateReg}}.
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
              results <- new("summaryRateReg", 
                             call = Call,
                             knots = knots,
                             boundaryKnots = boundaryKnots,
                             covariateCoef = beta,
                             frailtyPar = theta,
                             degree = object@degree,
                             baseRateCoef = alpha,
                             logL = object@logL)
              ## return
              results
          })

