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


#' Estimated Coefficients of Baseline Rate Function
#' 
#' An S4 class generic function to extract the estimated coefficients
#' of baseline rate function. It returns either coefficients
#' of pieceswise constant rate function
#' or coefficients of B-spline bases
#' 
#' @param object rateReg-class object.
#' @param ... Other arguments for future usage.
#' @return A named numeric vector, either pieceswise constants
#' or coefficients of B-spline bases.
#' @aliases BaseRate,rateReg-method
#' @examples
#' ## See examples given in \code{\link{rateReg}}
#' @seealso \code{\link{rateReg}} \code{\link{summary,rateReg-method}}
#' @importFrom methods setGeneric
#' @export
setGeneric(name = "baseRate",
           def = function(object, ...) {
               standardGeneric("baseRate")
           })


#' @describeIn baseRate Extract estiamted coefficients of
#' baseline rate function from \code{rateReg-class} object.
#' @importFrom methods setMethod
#' @export
setMethod(f = "baseRate", signature = "rateReg",
          definition = function(object, ...) {
              object@estimates$alpha[, 1]
          })

