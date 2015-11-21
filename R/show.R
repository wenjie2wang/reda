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


#' Show an object.
#' 
#' An S4 class generic function to display certain object.
#' 
#' \itemize{
#'   \item For \code{\link{rateReg-class}} object, 
#'   it prints brief summary of the fitted HEART model.
#'   \item For \code{\link{summaryHeart-class}} object, 
#'   it prints summary of the fitted HEART model.
#'   \item For \code{\link{empirMcf-class}} object,
#'   it prints formula and the first 100 rows of the computed MCF data frame.
#'   \item For \code{\link{rateRegMcf-class}} object,
#'   it prints formula, baseline pieces and the first 100 rows of the estimated
#'   MCF data frame.
#' }
#' 
#' @param object Certain R object produced by package reda.
#' @name show
#' @seealso \code{\link{rateReg}} \code{\link{summary,rateReg-method}}
#' \code{\link{mcf}}
NULL


#' @rdname show 
#' @aliases show,rateReg-method
#' @importFrom methods show
#' @export
setMethod(f = "show", signature = "rateReg",
          definition = function(object) {
              beta <- object@estimates$beta[, 1]
              theta <- object@estimates$theta[, 1]
              names(theta) <- NULL
              alpha <- object@estimates$alpha[, 1]
              cat("Call: \n")
              print(object@call)
              cat("\nCoefficients of covariates: \n") 
              print(beta)
              cat("\nFrailty parameter: ", theta, "\n")
              if (length(object@knots) > 0) {
              cat("\nInternal knots: \n") 
              cat(object@knots, sep = ", ", fill = TRUE)
              }
              cat("\nBoundary knots: \n")
              cat(object@boundaryKnots, sep = ", ", fill = TRUE)
              if (object@degree > 0) {
                  cat("\nCoefficients of spline bases:\n")
                  print(alpha)    
              } else {
                  cat("\nCoefficients of pieces:\n")
                  print(alpha)
              }
          })


#' @rdname show 
#' @aliases show,summaryHeart-method
#' @importFrom methods show
#' @importFrom stats printCoefmat
#' @export
setMethod(f = "show", signature = "summaryHeart",
          definition = function(object) {
              if (attr(object@call, "show")) {
                  Call <- object@call
                  attr(Call, "show") <- NULL
                  cat("Call: \n")
                  print(Call)
              }
              cat("\nCoefficients of covariates: \n") 
              printCoefmat(object@covariateCoef)
              cat("\nParameter of frailty: \n")
              print(object@frailtyPar)
              if (attr(object@knots, "show")) {
                  cat("\nInternal knots: \n")
                  cat(object@knots, sep = ", ", fill = TRUE)
                  cat("\nBoundary knots:\n")
                  cat(object@boundaryKnots, sep = ", ", fill = TRUE)
              }
              if (object@degree > 0) {
                  cat("\nDegree of spline bases:", object@degree, "\n")
                  cat("\nCoefficients of spline bases:\n")
                  printCoefmat(object@baseRateCoef)    
              } else {
                  cat("\nCoefficients of pieces:\n")
                  printCoefmat(object@baseRateCoef)    
              }
              cat("\nLoglikelihood: ", object@logL, "\n")
          })


#' @rdname show 
#' @aliases show,empirMcf-method 
#' @importFrom methods show 
#' @importFrom utils head
#' @export
setMethod(f = "show", signature = "empirMcf",
          definition = function(object) {
              cat("Call: \n")
              print(object@call)
              cat("\nFormula:\n")
              print(object@formula)
              cat("\nMCF:\n")
              print(object@MCF)
          })


#' @rdname show 
#' @aliases show,rateRegMcf-method
#' @importFrom methods show 
#' @importFrom utils head
#' @export
setMethod(f = "show", signature = "rateRegMcf",
          definition = function(object) {
              cat("Formula:\n")
              print(object@formula)
              cat("\nNew data:\n")
              print(object@newdata)
              cat("\nConfidence level:",
                  paste(format(100 * object@level,
                               trim = TRUE, scientific = FALSE),
                        "%", sep = ""), "\n")
              cat("\nMCF:\n")
              print(object@MCF)
          })
