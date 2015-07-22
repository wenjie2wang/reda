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


## create S4 Class called "summary.heart" for summary.heart object from summary
#' An S4 class to represent summary of heart object
#' @slot call call
#' @slot baselinepieces numeric vector
#' @slot coefficients numeric matrix
#' @slot theta numeric matrix
#' @slot baseline numeric matrix 
#' @export
setClass(Class = "summary.heart", 
         slots = c(call = "call", 
                   baselinepieces = "numeric",
                   coefficients = "matrix",
                   theta = "matrix",
                   baseline = "matrix"))


#' Summarizing HEART Model Fits
#'
#' \code{summary} returns summary of estimates from HEART model.
#'
#' These are details.
#'
#' @param object heart object from \code{heart}.
#' @param showcall logic value (TRUE or FALSE) with dafault as TRUE,
#' indicating whether method \code{show} for object summary.heart prints out 
#' the call information of original call of \code{heart}.
#' @param showpieces logic value (TRUE or FALSE) with default as TRUE, 
#' indicating whether method \code{show} for object summary.heart prints out 
#' the baseline pieces.
#' @param digits, an interger specifying the desired number of decimal places 
#' (round) for estimates. Negative values are allowed 
#' (\code{help(round)} for more details).
#' @return summary.heart object
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
