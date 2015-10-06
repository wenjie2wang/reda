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


#' Show an object.
#' 
#' An S4 class generic function to display certain object.
#' 
#' \itemize{
#'   \item For \code{\link{heart-class}} object, 
#'   it prints brief summary of the fitted HEART model.
#'   \item For \code{\link{summaryHeart-class}} object, 
#'   it prints summary of the fitted HEART model.
#'   \item For \code{\link{empirMcf-class}} object,
#'   it prints formula and the first 100 rows of the computed MCF data frame.
#'   \item For \code{\link{heartMcf-class}} object,
#'   it prints formula, baseline pieces and the first 100 rows of the estimated
#'   MCF data frame.
#' }
#' 
#' @param object certain R object generated from heart package.
#' @name show
#' @seealso \code{\link{heart}} \code{\link{summary,heart-method}}
#' \code{\link{mcf}}
NULL


#' @rdname show 
#' @aliases show,heart-method

#' @importFrom methods show
#' @export
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
              cat(attr(object@baselinePieces, "name"), "\n")
              cat("\ncoefficients: \n") 
              print(beta)
              cat("\ntheta: ", theta, "\n")
              cat("\nbaseline rate functions: \n")
              print(alpha)
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
                  cat("\ncall: \n")
                  print(Call)
              }
              if (attr(object@baselinePieces, "show")) {
                  cat("\nbaseline pieces: \n")
                  cat(attr(object@baselinePieces, "name"), "\n")
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
              cat("\ncall: \n")
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
#' @aliases show,heartMcf-method
#' @importFrom methods show 
#' @importFrom utils head
#' @export
setMethod(f = "show", signature = "heartMcf",
          definition = function(object) {
              cat("formula:\n")
              print(object@formula)
              cat("\nbaseline pieces:\n")
              print(attr(object@baselinePieces, "name"))
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
