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
#' @param object rateReg-class object.
#' @param ... other arguments for future usage.
#' @return a named vector.
#' @aliases baseline,rateReg-method
#' @examples 
#' library(reda)
#' rateRegFit <- rateReg(formula = Survr(ID, time, event) ~ x1 + group, 
#'                   data = simuDat, subset = ID %in% 75:125,
#'                   baselinePieces = seq(28, 168, length = 6))
#' baseline(rateRegFit)
#' @seealso \code{\link{rateReg}} \code{\link{summary,rateReg-method}}
#' @export
setGeneric(name = "baseline",
           def = function(object, ...) {
               standardGeneric("baseline")
           })


#' @describeIn baseline Extract estiamted baseline rate function 
#' from rateReg-class object.
#' @export
setMethod(f = "baseline", signature = "rateReg",
          definition = function(object, ...) {
              alpha <- object@estimates$alpha[, "alpha"]
              names(alpha) <- rownames(object@estimates$alpha)
              ## return
              alpha
          })

