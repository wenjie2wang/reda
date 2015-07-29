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


## create S4 Class called "heart" for heart object from function heart
#' An S4 Class to Represent a Fitted HEART Model
#' 
#' heart-class is an S4 class to represent a HEART model fits. 
#' \code{\link{heart}} produces objects of this class.  
#' See ``Slots'' for details.
#' 
#' @slot call function call.
#' @slot formula formula.
#' @slot baselinepieces a numeric vector.
#' @slot estimates list.
#' @slot control list.
#' @slot start list.
#' @slot na.action a length-one character vector.
#' @slot xlevels list.
#' @slot contrasts list.
#' @slot convergence an integer.
#' @slot hessian a numeric matrix.
#' @aliases heart-class
#' @seealso \code{\link{heart}}
#' @export
setClass(Class = "heart", 
         slots = c(call = "call", 
                   formula = "formula", 
                   baselinepieces = "numeric",
                   estimates = "list",
                   control = "list",
                   start = "list",
                   na.action = "character",
                   xlevels = "list",
                   contrasts = "list",
                   convergence = "integer", 
                   hessian = "matrix"))


## create S4 Class called "summary.heart" for summary.heart object from summary
#' An S4 Class to Represent Summary of heart-class Object
#' 
#' summary.heart-class is an S4 class with selective slots 
#' of heart-class object.  See ``Slots'' for details.  
#' \code{\link{summary}} produces objects of this class. 
#'  
#' @slot call function call.
#' @slot baselinepieces a numeric vector.
#' @slot coefficients a numeric matrix.
#' @slot theta numeric a matrix.
#' @slot baseline a numeric matrix.
#' @aliases summary.heart-class
#' @seealso \code{\link{summary,heart-method}} 
#' @export
setClass(Class = "summary.heart", 
         slots = c(call = "call", 
                   baselinepieces = "numeric",
                   coefficients = "matrix",
                   theta = "matrix",
                   baseline = "matrix"))


#' An S4 Class to Represent Computed Empirical MCF
#' 
#' An S4 class to represent computed empirical mean cumulative function (MCF).
#' \code{\link{MCF}} produces objects of this class.  
#' @slot call function call
#' @slot formula formula. 
#' @slot MCF a data.frame.
#' @slot multigroup a logical value.
#' @aliases empirMCF-class
#' @seealso \code{\link{MCF,formula-method}}
#' @importFrom methods setClass
#' @export
setClass(Class = "empirMCF", 
         slots = c(call = "call", formula = "formula", MCF = "data.frame", 
                   multigroup = "logical"))


#' An S4 Class to Represent Estimated MCF from HEART Model
#' 
#' An S4 class to represent estimated mean cumulative function (MCF) 
#' from HEART Model.
#' \code{\link{MCF}} produces objects of this class.  
#' 
#' @slot formula formula.
#' @slot baselinepieces a numeric vector.
#' @slot newdata a numeric matrix.
#' @slot MCF a data.frame.
#' @slot level a numeric value between 0 and 1.
#' @slot na.action a length-one character vector.
#' @slot control list.
#' @slot multigroup a logical value. 
#' @aliases heartMCF-class
#' @seealso \code{\link{MCF,heart-method}}
#' @export
setClass(Class = "heartMCF", 
         slots = c(formula = "formula", baselinepieces = "numeric",
                   newdata = "matrix", MCF = "data.frame", level = "numeric", 
                   na.action = "character", control = "list", 
                   multigroup = "logical"))

