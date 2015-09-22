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


#' Estimated Baseline Rate Function
#' 
#' An S4 class generic function to extract the estimated baseline rate function 
#' from HEART model. 
#' 
#' @param object heart-class object.
#' @param digits an integer indicating the number of decimal places to be used. 
#' Negative values are allowed (see 'Details' of \code{\link{round}}).
#' The default value is 3.
#' @param ... other arguments for future usage.
#' @return a named vector.
#' @aliases baseline,heart-method
#' @examples 
#' library(heart)
#' data(simuDat)
#' heartfit <- heart(formula = Survr(ID, time, event) ~ X1 + group, 
#'                   data = simuDat, subset = ID %in% 75:125,
#'                   baselinepieces = seq(28, 168, length = 6))
#' baseline(heartfit)
#' @seealso \code{\link{heart}} \code{\link{summary,heart-method}}
#' @export
setGeneric(name = "baseline",
           def = function(object, ...) {
             standardGeneric("baseline")
           })


#' @describeIn baseline Extract estiamted baseline rate function 
#' from heart-class object.
#' @export
setMethod(f = "baseline", signature = "heart",
          definition = function(object, digits = 3, ...) {
            alpha <- round(object@estimates$alpha[, "alpha"], digits = digits)
            names(alpha) <- rownames(object@estimates$alpha)
            ## return
            alpha
          })



