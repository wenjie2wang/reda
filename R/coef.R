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

#' Extract Coefficients from HEART Model
#'
#' coef is a S4 class generic function 
#' which extracts model coefficients from objects returned by modeling 
#' functions. 
#'
#' These are details
#' 
#' @usage 
#' coef(object)
#' 
#' @param object heart object
#' @param ... other arguments
#' @return Coefficients extracted from the model object object. 
#' For heart object, it will be a named numeric vector.
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

#' function confint for heart object
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


