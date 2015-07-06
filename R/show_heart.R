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
## function show for heart object
setMethod(f = "show", signature = "heart",
          definition = function(object) {
            beta <- round(object@estimates$beta[, "coef"], digits = 3)
            names(beta) <- rownames(object@estimates$beta)
            theta <- round(object@estimates$theta[, "theta"], digits = 3)
            names(theta) <- NULL
            alpha <- round(object@estimates$alpha[, "alpha"], digits = 3)
            names(alpha) <- rownames(object@estimates$alpha)
            cat("\ncall: \n")
            print(object@call)
            cat("\nbaseline pieces: \n")
            cat(attr(object@baselinepieces, "name"), "\n")
            cat("\ncoefficients: \n") 
            print(beta)
            cat("\ntheta: ", theta, "\n")
            cat("\nbaseline rate functions: \n")
            print(alpha)
          })

## function show for summary.heart object
setMethod(f = "show", signature = "summary.heart",
          definition = function(object) {
            if (attr(object@call, "show")) {
              Call <- object@call
              attr(Call, "show") <- NULL
              cat("\ncall: \n")
              print(Call)
            }
            if (attr(object@baselinepieces, "show")) {
              cat("\nbaseline pieces: \n")
              cat(attr(object@baselinepieces, "name"), "\n")
            }
            cat("\ncoefficients: \n") 
            printCoefmat(object@coefficients)
            theta <- as.data.frame(object@theta)
            cat("\ntheta: \n")
            print(theta, row.names = FALSE)
            cat("\nbaseline rate functions: \n")
            print(object@baseline)
          })
