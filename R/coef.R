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

#' Extract coefficients estiamtes from HEART model.
#'
#' \code{coef} is a S4 class generic function 
#' which extracts model coefficients from the heart object returned by modeling 
#' function \code{heart}. 
#'
#' @usage 
#' coef(object)
#' @param object heart object.
#' @param ... other arguments.
#' @return A named numeric vector.
#' @examples
#' data(simuDat)
#' heartfit <- heart(formula = Survr(ID, time, event) ~ X1 + group, 
#'                  data = simuDat, baselinepieces = seq(28, 168, length = 5))
#' coef(heartfit)
#' @export
setMethod(f = "coef", signature = "heart",
          definition = function(object, ...) {
            beta <- round(object@estimates$beta[, "coef"], digits = 3)
            names(beta) <- rownames(object@estimates$beta)
            ## return
            beta
          })



#' Confidence Intervals for HEART Model Coefficients
#'
#' \code{confint} is a S4 class generic function for heart object, 
#' which computes confidence intervals 
#' for all or selected coefficients in a fitted HEART model. 
#'
#' Under regularity condition (Shao, 2003, 
#' Theorem 4.16 and Theorem 4.17, page 287, 290), 
#' the confidence intervals are constructed loosely 
#' based on Fisher information matrix and estimates of coefficients. 
#' See \emph{Fu et al. (2014)} for more details.
#' 
#' @usage 
#' confint(object, parm, level = 0.95, ...)
#' @param object heart object.
#' @param parm a specification of which parameters are 
#' to be given confidence intervals, 
#' either a vector of numbers or a vector of names. 
#' If missing, all parameters are considered.
#' @param level the confidence level required.
#' @param ... additional argument(s).
#' @return A numeric matrix with rownames and colnames.
#' @references 
#' Shao, J. (2003), 
#' \emph{Mathematical statistics}, Springer texts in statistics, 
#' New York: Springer, 2nd edition.
#' 
#' Fu, Haoda, Junxiang Luo, and Yongming Qu. (2014),
#' "Hypoglycemic Events Analysis via Recurrent Time-to-Event (HEART) Models," 
#' \emph{Journal of biopharmaceutical statistics}, 2014 Dec 1, Epub 2014 Dec 1.
#' @examples
#' data(simuDat)
#' heartfit <- heart(formula = Survr(ID, time, event) ~ X1 + group, 
#'                   data = simuDat, baselinepieces = seq(28, 168, length = 5))
#' confint(heartfit)
#' @export
setMethod(f = "confint", signature = "heart",
          definition = function(object, parm, level = 0.95, ...) {
            ## internal function
            format.perc <- function (probs, digits){
              paste(format(100 * probs, trim = TRUE, scientific = FALSE, 
                           digits = digits), "%")
            }
            betamat <- object@estimates$beta
            cf <- betamat[, 1]
            pnames <- attr(betamat, "dimnames")[[1]]
            if (missing(parm)) {
              parm <- seq(nrow(betamat))
            } else if (is.numeric(parm)) { 
              parm <- intersect(seq(nrow(betamat)), parm) 
            } else if (is.character(parm)) {
              parm <- match(parm, pnames, nomatch = NULL)
            } else {
              stop("invalid argument param")
            }
            a <- (1 - level)/2
            a <- c(a, 1 - a)
            fac <- qnorm(a)
            pct <- format.perc(a, 3)
            ci <- array(NA, dim = c(length(parm), 2L), 
                        dimnames = list(parm, pct))
            ses <- betamat[parm, 2]
            ci[] <- cf[parm] + ses %o% fac
            ci <- round(ci, digits = 3)
            rownames(ci) <- pnames[parm]
            ## return
            ci
          })

