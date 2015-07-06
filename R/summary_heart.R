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

## internal functions that do not export =======================================

## create S4 Class called "summary.heart" for object from summary to show
setClass(Class = "summary.heart", 
         slots = c(call = "call", 
                   baselinepieces = "numeric",
                   coefficients = "matrix",
                   theta = "matrix",
                   baseline = "matrix"))



#' Fitting Heart Model: Piece-wise Gamma Frailty Model for Recurrent Events. 
#' The model is named after the paper title of \emph{Fu et al. (2014)},  
#' Hypoglycemic Events Analysis via Recurrent Time-to-Event (HEART) Models
#'
#' \code{functionname} returns fitted model results.
#'
#' This is a test Roxygen comments
#'
#' @param ... Numeric, complex, or logical vectors.
#' @param na.rm A logical scalar. 
#' @return If all inputs are integer and logical, then the output
#'   will be an integer. If integer overflow
#'   \url{http://en.wikipedia.org/wiki/Integer_overflow} occurs, the output
#'   will be NA with a warning. Otherwise it will be a length-one numeric or
#'   complex vector.
#' @examples
#' sum(1:10)
#' sum(1:5, 6:10)
#' sum(F, F, F, T, T)
## function summary for heart object
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
