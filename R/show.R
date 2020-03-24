##
## R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
## Copyright (C) 2015-2020
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


##' Show an object.
##'
##' S4 class methods that display objects produced from this package (similar to
##' S3 class \code{print} methods).
##'
##' @param object An object used to dispatch a method.
##' @name show-method
##' @importFrom methods show
NULL


##' @rdname show-method
##' @aliases show,Recur-method
##' @export
setMethod(f = "show", signature = "Recur",
          definition = function(object) {
              invisible(print(as.character(object), quote = FALSE))
          })


##' @rdname show-method
##' @aliases show,rateReg-method
##' @export
setMethod(f = "show", signature = "rateReg",
          definition = function(object) {
              beta <- object@estimates$beta[, "coef"]
              names(beta) <- row.names(object@estimates$beta)
              theta <- object@estimates$theta[, "parameter"]
              names(theta) <- NULL
              alpha <- object@estimates$alpha[, "coef"]
              names(alpha) <- row.names(object@estimates$alpha)
              cat("Call: \n")
              print(object@call)
              if (length(beta) > 0) {
                  cat("\nCoefficients of covariates: \n")
                  print(beta)
              }
              cat("\nFrailty parameter: ", theta, "\n")
              knots <- object@spline$knots
              Boundary.knots <- object@spline$Boundary.knots
              if (length(knots)) {
                  cat("\nInternal knots: \n")
                  cat(knots, sep = ", ", fill = TRUE)
              }
              cat("\nBoundary knots: \n")
              cat(Boundary.knots, sep = ", ", fill = TRUE)
              if (object@spline$degree) {
                  cat("\nCoefficients of spline bases:\n")
                  print(alpha)
              } else {
                  cat("\nCoefficients of pieces:\n")
                  print(alpha)
              }
              ## invisible return
              invisible(object)
          })


##' @rdname show-method
##' @aliases show,summaryRateReg-method
##' @importFrom stats printCoefmat
##' @export
setMethod(f = "show", signature = "summary.rateReg",
          definition = function(object) {
              if (attr(object@call, "show")) {
                  Call <- object@call
                  attr(Call, "show") <- NULL
                  cat("Call: \n")
                  print(Call)
              }
              if (nrow(object@covarCoef) > 0) {
                  cat("\nCoefficients of covariates: \n")
                  printCoefmat(object@covarCoef)
              }
              cat("\nParameter of frailty: \n")
              print(object@frailtyPar)
              ## on knots
              if (attr(object@knots, "show")) {
                  if (length(object@knots)) {
                      cat("\nInternal knots: \n")
                      cat(object@knots, sep = ", ", fill = TRUE)
                  }
                  cat("\nBoundary knots:\n")
                  cat(object@Boundary.knots, sep = ", ", fill = TRUE)
              }
              ## baseline rate function
              cat("\nDegree of spline bases:", object@degree, "\n")
              cat("\nCoefficients of spline bases:\n")
              printCoefmat(object@baseRateCoef)
              cat("\nLoglikelihood: ", object@logL, "\n")
              ## invisible return
              invisible(object)
          })


##' @rdname show-method
##' @aliases show,mcf.formula-method
##' @export
setMethod(f = "show", signature = "mcf.formula",
          definition = function(object) {
              cat("Formula:\n")
              print(object@formula)
              cat("\nMCF:\n")
              print(object@MCF)
              ## invisible return
              invisible(object)
          })


##' @rdname show-method
##' @aliases show,mcf.rateReg-method
##' @export
setMethod(
    f = "show",
    signature = "mcf.rateReg",
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
        ## invisible return
        invisible(object)
    }
)


##' @rdname show-method
##' @aliases show,simEvent-method
##' @export
setMethod(
    f = "show",
    signature = "simEvent",
    definition = function(object) {
        cat("'simEvent' S4 class object:\n")
        print(object@.Data)
        ## invisible return
        invisible(object)
    }
)


##' @rdname show-method
##' @aliases show,mcfDiff-method
##' @export
setMethod(
    f = "show",
    signature = "mcfDiff",
    definition = function(object) {
        cat("Call: \n")
        print(object@call)
        if (object@test@testVariance != "none") {
            cat("\nTwo-Sample Pseudo-Score Tests:\n")
            printCoefmat(object@test@.Data)
            cat("\nVariance Estimator:", object@test@testVariance, "\n")
        }
        ## invisible return
        invisible(object)
    }
)


##' @rdname show-method
##' @aliases show,mcfDiff.test-method
##' @export
setMethod(
    f = "show",
    signature = "mcfDiff.test",
    definition = function(object) {
        cat("Two-Sample Pseudo-Score Tests:\n")
        printCoefmat(object@.Data)
        cat("\nVariance Estimator:", object@testVariance, "\n")
        ## invisible return
        invisible(object)
    })
