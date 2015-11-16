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
#' @param object certain R object generated from rateReg package.
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
              beta <- round(object@estimates$beta[, "coef"], digits = 3)
              names(beta) <- rownames(object@estimates$beta)
              theta <- round(object@estimates$theta[, "theta"], digits = 3)
              names(theta) <- NULL
              alpha <- round(object@estimates$alpha[, "alpha"], digits = 3)
              names(alpha) <- rownames(object@estimates$alpha)
              cat("call: \n")
              print(object@call)
              cat("\ncoefficients of covariates: \n") 
              print(beta)
              cat("\ntheta: ", theta, "\n")
              if (length(object@knots) > 0) {
              cat("\nknots: \n", object@knots, "\n")
              }
              if (object@degree > 0) {
                  cat("\ncoefficients of spline bases:\n")
                  print(alpha)    
              } else {
                  cat("\ncoefficients of pieces:\n")
                  print(alpha)
              }
          })


#' @rdname show 
#' @aliases show,summaryHeart-method
#' @importFrom methods show
#' @export
setMethod(f = "show", signature = "summaryHeart",
          definition = function(object) {
              if (attr(object@call, "show")) {
                  Call <- object@call
                  attr(Call, "show") <- NULL
                  cat("call: \n")
                  print(Call)
              }
              if (attr(object@knots, "show")) {
                  cat("\nbaseline pieces: \n")
                  cat(attr(object@knots, "name"), "\n")
              }
              cat("\ncoefficients: \n") 
              printCoefmat(object@coefficients)
              theta <- as.data.frame(object@theta)
              cat("\ntheta: \n")
              print(theta, row.names = FALSE)
              cat("\nbaseline rate functions: \n")
              print(object@baseline)
          })


#' @rdname show 
#' @aliases show,empirMcf-method 
#' @importFrom methods show 
#' @importFrom utils head
#' @export
setMethod(f = "show", signature = "empirMcf",
          definition = function(object) {
              cat("call: \n")
              print(object@call)
              cat("\nformula:\n")
              print(object@formula)
              cat("\nMCF:\n")
              if (nrow(object@MCF) <= 100) {
                  print(object@MCF)
                  cat("\n")
              } else {
                  print(head(object@MCF, 100))
                  cat("...\n\n")
                  cat("Only the first 100 rows are printed.\n")
              }
          })


#' @rdname show 
#' @aliases show,rateRegMcf-method
#' @importFrom methods show 
#' @importFrom utils head
#' @export
setMethod(f = "show", signature = "rateRegMcf",
          definition = function(object) {
              cat("formula:\n")
              print(object@formula)
              cat("\nbaseline pieces:\n")
              print(attr(object@knots, "name"))
              cat("\nMCF:\n")
              if (nrow(object@MCF) <= 100) {
                  print(object@MCF)
                  cat("\n")
              } else {
                  print(head(object@MCF, 100))
                  cat("...\n\n")
                  cat("Only the first 100 rows are printed.\n")
              }
          })
