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


#' Summarizing HEART Model Fits
#'
#' \code{summary} returns summary of estimates from HEART model.
#'
#' To be specific, \code{summary} returns a \code{\link{summary.heart-class}}
#' object which can be printed by \code{\link{show,summary.heart-method}}. 
#'
#' @param object heart object from \code{heart}.
#' @param showcall a logic value with dafault as TRUE,
#' indicating whether method \code{\link{show,summary.heart-method}} prints out 
#' the call information of original call of \code{heart}.
#' @param showpieces a logic value with default as TRUE, 
#' indicating whether method \code{\link{show,summary.heart-method}} prints out 
#' the baseline pieces.
#' @param digits, an interger specifying the desired number of decimal places 
#' (round) for estimates. Negative values are allowed 
#' (\code{help(round)} for more details).
#' @return summary.heart-class object
#' @aliases summary,heart-method
#' @seealso \code{\link{heart}} \code{\link{coef,heart-method}}
#' \code{\link{confint,heart-method}} \code{\link{baseline,heart-method}}
#' \code{\link{mcf}}
#' @importFrom methods new
#' @export
setMethod(f = "summary", signature = "heart",
          definition = function(object, showcall = TRUE, showpieces = TRUE, 
                                digits = 3) {
              Call <- object@call
              attr(Call, "show") <- showcall
              blpieces <- object@baselinepieces
              attr(blpieces, "show") <- showpieces
              beta <- round(object@estimates$beta, digits = digits)
              theta <- round(object@estimates$theta, digits = digits)
              alpha <- round(object@estimates$alpha, digits = digits)
              colnames(beta)[1] <- colnames(theta)[1] <- 
                  colnames(alpha)[1] <- "estimates"
              results <- new("summary.heart", 
                             call = Call,
                             baselinepieces = blpieces,
                             coefficients = beta,
                             theta = theta, 
                             baseline = alpha)
              ## return
              results
          })

