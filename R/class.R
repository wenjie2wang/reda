################################################################################
##
##   R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
##   Copyright (C) 2015
##
##   This file is part of the R package reda.
##
##   The R package reda is free software: You can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   any later version (at your option). See the GNU General Public License
##   at <http://www.gnu.org/licenses/> for details.
##
##   The R package reda is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
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
#' \code{survrec} but with a better checking procedure for recurrent event
#' data modeled by methods based on counts and rate function.
#' The checking rules include that
#' \itemize{
#'     \item Identificator of each subject cannot be missing.
#'     \item Event indicator must be coded as 0 (censoring) or 1 (event).
#'     \item Event time and censoring time cannot be missing.
#'     \item Each subject must have one and only one censoring time.
#'     \item Event time cannot not be later than censoring time.
#' }
#'  
#' @param ID Identificator of each subject. 
#' @param time Time of reccurence event or censoring.
#' @param event The status indicator, 0 = censored, 1 = event. 
#' @aliases Survr
#' @seealso \code{\link{rateReg}}
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
#' @slot df List.
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
                   df = "list",
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
#' \code{summaryHeart-class} is an S4 class with selective slots 
#' of \code{rateReg-class} object.  See ``Slots'' for details.  
#' \code{\link{summary}} produces objects of this class. 
#'  
#' @slot call Function call.
#' @slot knots A numeric vector.
#' @slot boundaryKnots A numeric vector.
#' @slot covariateCoef A numeric matrix.
#' @slot frailtyPar A numeric matrix.
#' @slot degree An integer.
#' @slot baseRateCoef A numeric matrix.
#' @slot logL A numeric value.
#' @aliases summaryHeart-class
#' @seealso \code{\link{summary,rateReg-method}} 
#' @export
setClass(Class = "summaryRateReg", 
         slots = c(call = "call", 
                   knots = "numeric",
                   boundaryKnots = "numeric",
                   covariateCoef = "matrix",
                   frailtyPar = "matrix",
                   degree = "integer",
                   baseRateCoef = "matrix",
                   logL = "numeric"))


#' An S4 Class to Represent Sample MCF
#' 
#' An S4 class to represent sample mean cumulative function (MCF).
#' \code{\link{mcf}} produces objects of this class.  
#' @slot call Function call
#' @slot formula Formula. 
#' @slot MCF A data frame.
#' @slot multiGroup A logical value.
#' @aliases sampleMcf-class
#' @seealso \code{\link{mcf,formula-method}}
#' @importFrom methods setClass
#' @export
setClass(Class = "sampleMcf", 
         slots = c(call = "call",
                   formula = "formula",
                   MCF = "data.frame", 
                   multiGroup = "logical"))


#' An S4 Class to Represent Estimated MCF from HEART Model
#' 
#' An S4 class to represent estimated mean cumulative function (MCF) 
#' from HEART Model.
#' \code{\link{mcf}} produces objects of this class.  
#' 
#' @slot call Function call.
#' @slot formula Formula.
#' @slot knots A numeric vector.
#' @slot degree An integer.
#' @slot boundaryKnots A numeric vector.
#' @slot newdata A numeric matrix.
#' @slot MCF A data frame.
#' @slot level A numeric value between 0 and 1.
#' @slot na.action A length-one character vector.
#' @slot control List.
#' @slot multiGroup A logical value. 
#' @aliases rateRegMcf-class
#' @seealso \code{\link{mcf,rateReg-method}}
#' @export
setClass(Class = "rateRegMcf", 
         slots = c(call = "call",
                   formula = "formula",
                   knots = "numeric",
                   degree = "integer",
                   boundaryKnots = "numeric",
                   newdata = "matrix",
                   MCF = "data.frame",
                   level = "numeric", 
                   na.action = "character",
                   control = "list", 
                   multiGroup = "logical"))

