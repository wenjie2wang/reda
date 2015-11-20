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


#' Estimated Coefficients of Covariates
#'
#' \code{coef} is a S4 class method which extracts
#' estimated coefficients of covariates
#' from the \code{\link{rateReg-class}} object 
#' produced by function \code{\link{rateReg}}.
#'
#' @param object rateReg-class object.
#' @param ... Other arguments for future usage.
#' @return A named numeric vector.
#' @aliases coef,rateReg-method
#' @seealso \code{\link{rateReg}} \code{\link{summary,rateReg-method}}
#' @examples
#' ## Please see examples given in \code{\link{rateReg}}.
#' @importFrom methods setMethod
#' @importFrom stats coef
#' @export
setMethod(f = "coef", signature = "rateReg",
          definition = function(object, ...) {
              object@estimates$beta[, 1]
          })


#' Confidence Intervals for Coefficients of Covariates
#' \code{confint} is a S4 class method for rateReg object, 
#' which computes confidence intervals 
#' for all or selected covariates. 
#'
#' Under regularity condition (Shao, 2003, 
#' Theorem 4.16 and Theorem 4.17, page 287, 290), 
#' the confidence intervals are constructed loosely 
#' based on Fisher information matrix and estimates of coefficients. 
#' See \emph{Fu et al. (2014)} for more details.
#' 
#' @param object rateReg-class object.
#' @param parm A specification of which parameters are 
#' to be given confidence intervals, 
#' either a vector of numbers or a vector of names. 
#' If missing, all parameters are considered.
#' @param level A optional numeric value to specify
#' the confidence level required.
#' By default, the value is 0.95,
#' which specifies 95\% confidence intervals.
#' @param ... other arguments for future usage.
#' @return a numeric matrix with rownames and colnames.
#' @aliases confint,rateReg-method
#' @seealso \code{\link{rateReg}} \code{\link{coef,rateReg-method}}
#' @references 
#' Fu, Haoda, Junxiang Luo, and Yongming Qu. (2014),
#' "Hypoglycemic Events Analysis via Recurrent Time-to-Event (HEART) Models,"
#'  \emph{Journal of biopharmaceutical statistics}, 2014 Dec 1, Epub 2014 Dec 1.
#' 
#' Shao, J. (2003), 
#' \emph{Mathematical statistics}, Springer texts in statistics, 
#' New York: Springer, 2nd edition.
#' @examples
#' ## Pleas see examples given in \code{\link{rateReg}}.
#' @importFrom stats confint qnorm 
#' @export
setMethod(f = "confint", signature = "rateReg",
          definition = function(object, parm, level = 0.95, ...) {
              ## internal function
              format.perc <- function (probs){
                  paste(format(100 * probs, trim = TRUE,
                               scientific = FALSE), "%", sep = "")
              }
              betaMat <- object@estimates$beta
              estCoef <- betaMat[, 1]
              pnames <- attr(betaMat, "dimnames")[[1]]
              if (missing(parm)) {
                  parm <- seq(nrow(betaMat))
              } else if (is.numeric(parm)) { 
                  parm <- intersect(seq(nrow(betaMat)), parm) 
              } else if (is.character(parm)) {
                  parm <- match(parm, pnames, nomatch = NULL)
              } else {
                  stop("invalid argument param")
              }
              a <- (1 + c(-1, 1) * level)/2
              fac <- qnorm(a)
              pct <- format.perc(a)
              ci <- array(NA, dim = c(length(parm), 2L), 
                          dimnames = list(parm, pct))
              ses <- betaMat[parm, 2]
              ci[] <- estCoef[parm] + ses %o% fac
              rownames(ci) <- pnames[parm]
              ## return
              ci
          })



