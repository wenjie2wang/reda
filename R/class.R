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


#' An S3 Class to Represent Formula Response for Recurrent Event Data
#' 
#' \code{Survr} is an S3 class to represent 
#' formula response for recurrent event data
#' modeled by methods based on counts and rate function.
#' The last letter 'r' in 'Survr' represents 'rate'.
#'
#' This is a similar function to \code{\link[survrec]{Survr}} in package
#' \code{survrec}, but with a better checking procedure for recurrent event
#' data. The checking rules include that
#' \itemize{
#'     \item Identificator of each subject cannot be missing.
#'     \item Event indicator must be coded as 0 (censoring) or 1 (event).
#'     \item Event time and censoring time cannot be missing.
#'     \item Each subject must have one and only one censoring time.
#'     \item Event time cannot not be later than censoring time.
#' }
#'  
#' @param ID identificator of each subject. 
#' @param time time of reccurence. For each subject the last time are censored.
#' @param event the status indicator, 
#' 0 = censored, 1 = event. 
#' @aliases Survr
#' @seealso \code{\link{rateReg}}
#' @importFrom plyr ddply
#' @export
Survr <- function (ID, time, event) {
    inpDat <- data.frame(ID, time, event)
    dat <- check_Survr(inpDat)
    outDat <- with(dat, as.matrix(cbind(ID, time, event)))
    attr(outDat, "ID") <- attr(dat, "ID")
    oldClass(outDat) <- "Survr"
    invisible(outDat)
}


#' An S4 Class to Represent a Fitted Model
#' 
#' \code{rateReg-class} is an S4 class to represent a model fits. 
#' \code{\link{rateReg}} produces objects of this class.  
#' See ``Slots'' for details.
#' 
#' @slot call Function call.
#' @slot formula Formula.
#' @slot knots A numeric vector.
#' @slot boundaryKnots A numeric vector.
#' @slot degree An integer. 
#' @slot estimates List.
#' @slot control List.
#' @slot start List.
#' @slot na.action A length-one character vector.
#' @slot xlevels List.
#' @slot contrasts List.
#' @slot convergCode An integer.
#' @slot logL A numeric value.
#' @slot fisher A numeric matrix.
#' @aliases rateReg-class
#' @seealso \code{\link{rateReg}}
#' @export
setClass(Class = "rateReg", 
         slots = c(call = "call", 
                   formula = "formula", 
                   knots = "numeric",
                   boundaryKnots = "numeric",
                   degree = "integer",
                   df = "integer",
                   estimates = "list",
                   control = "list",
                   start = "list",
                   na.action = "character",
                   xlevels = "list",
                   contrasts = "list",
                   convergCode = "integer",
                   logL = "numeric",
                   fisher = "matrix"))


## create S4 Class called "summaryHeart" for summaryHeart object from summary
#' An S4 Class to Represent Summary of rateReg-class Object
#' 
#' summaryHeart-class is an S4 class with selective slots 
#' of rateReg-class object.  See ``Slots'' for details.  
#' \code{\link{summary}} produces objects of this class. 
#'  
#' @slot call function call.
#' @slot baselinePieces a numeric vector.
#' @slot coefficients a numeric matrix.
#' @slot theta numeric a matrix.
#' @slot baseline a numeric matrix.
#' @aliases summaryHeart-class
#' @seealso \code{\link{summary,rateReg-method}} 
#' @export
setClass(Class = "summaryHeart", 
         slots = c(call = "call", 
                   knots = "numeric",
                   boundaryKnots = "numeric",
                   covariateCoef = "matrix",
                   frailtyPar = "matrix",
                   degree = "integer",
                   df = "integer",
                   baseRateCoef = "matrix",
                   logL = "numeric"))


#' An S4 Class to Represent Computed Empirical MCF
#' 
#' An S4 class to represent computed empirical mean cumulative function (MCF).
#' \code{\link{mcf}} produces objects of this class.  
#' @slot call function call
#' @slot formula formula. 
#' @slot MCF a data.frame.
#' @slot multiGroup a logical value.
#' @aliases empirMcf-class
#' @seealso \code{\link{mcf,formula-method}}
#' @importFrom methods setClass
#' @export
setClass(Class = "empirMcf", 
         slots = c(call = "call", formula = "formula", MCF = "data.frame", 
                   multiGroup = "logical"))


#' An S4 Class to Represent Estimated MCF from HEART Model
#' 
#' An S4 class to represent estimated mean cumulative function (MCF) 
#' from HEART Model.
#' \code{\link{mcf}} produces objects of this class.  
#' 
#' @slot formula formula.
#' @slot baselinePieces a numeric vector.
#' @slot newdata a numeric matrix.
#' @slot MCF a data.frame.
#' @slot level a numeric value between 0 and 1.
#' @slot na.action a length-one character vector.
#' @slot control list.
#' @slot multiGroup a logical value. 
#' @aliases rateRegMcf-class
#' @seealso \code{\link{mcf,rateReg-method}}
#' @export
setClass(Class = "rateRegMcf", 
         slots = c(formula = "formula",
                   knots = "numeric",
                   degree = "integer",
                   boundaryKnots = "numeric",
                   newdata = "matrix",
                   MCF = "data.frame",
                   level = "numeric", 
                   na.action = "character",
                   control = "list", 
                   multiGroup = "logical"))

