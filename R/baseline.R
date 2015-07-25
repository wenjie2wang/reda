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


## collation after class.R
#' @include class.R 
NULL


#' Estimated Baseline Rate Function
#' 
#' An S4 class generic function to extract the estimated baseline rate function 
#' from HEART model. 
#' 
#' @param object heart object.
#' @param ... further arguments.
#' @return a named vector.
#' @export
setGeneric(name = "baseline",
           def = function(object, ...) {
             standardGeneric("baseline")
           })


#' @describeIn baseline Extract estiamted baseline rate function 
#' from heart object.
#' @export
setMethod(f = "baseline", signature = "heart",
          definition = function(object, ...) {
            alpha <- round(object@estimates$alpha[, "alpha"], digits = 3)
            names(alpha) <- rownames(object@estimates$alpha)
            ## return
            alpha
          })



