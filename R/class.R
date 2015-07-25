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
#' An S4 class to represent a fitted HEART model.
#' @slot call call
#' @slot formula formula
#' @slot baselinepieces a numeric vector
#' @slot estimates list 
#' @slot control list
#' @slot start list
#' @slot na.action a length-one character vector
#' @slot xlevels list
#' @slot contrasts list
#' @slot convergence an integer
#' @slot hessian a numeric matrix
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
#' An S4 class to represent summary of heart object
#' @slot call call
#' @slot baselinepieces numeric vector
#' @slot coefficients numeric matrix
#' @slot theta numeric matrix
#' @slot baseline numeric matrix 
#' @export
setClass(Class = "summary.heart", 
         slots = c(call = "call", 
                   baselinepieces = "numeric",
                   coefficients = "matrix",
                   theta = "matrix",
                   baseline = "matrix"))


#' An S4 class to represent computed empirical MCF
#' @slot formula formula 
#' @slot MCF data.frame
#' @slot multigroup logical value
#' @importFrom methods setClass
#' @export
setClass(Class = "empirMCF", 
         slots = c(formula = "formula", MCF = "data.frame", 
                   multigroup = "logical"))


#' An S4 class to represent estimated MCF from HEART model
#' @slot formula formula
#' @slot baselinepieces numeric vector.
#' @slot newdata numeric matrix.
#' @slot MCF data.frame.
#' @slot level a numeric value within (0, 1).
#' @slot na.action a length-one character vector.
#' @slot control list.
#' @slot multigroup logical. 
#' @export
setClass(Class = "heartMCF", 
         slots = c(formula = "formula", baselinepieces = "numeric",
                   newdata = "matrix", MCF = "data.frame", level = "numeric", 
                   na.action = "character", control = "list", 
                   multigroup = "logical"))

