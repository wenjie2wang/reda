################################################################################
##
##   R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
##   Copyright (C) 2015-2017
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


## collation after class.R
##' @include class.R
NULL


##' Mean Cumulative Function (MCF)
##'
##' An S4 class generic function that estimates mean cumulative function (MCF)
##' from a fitted model or computing the sample nonparametric MCF
##' (also called Nelson-Aalen estimator) from data.
##'
##' For \code{formula} object with \code{\link{Survr}} object as response,
##' the covariate specified at the right hand side of the formula
##' should be either 1 or any "linear" conbination of factor variable
##' in the data.
##' The former computes the overall sample MCF.
##' The latter computes the sample MCF for each level of the combination of
##' the factor variable(s) specified, respectively.
##' The sample MCF is also called Nelson-Aalen nonparametric estimator
##' (Nelson, 2003) and computed on each time point from sample data.
##' The point estimate of sample MCF at each time point does not
##' assume any particular underlying model. The variance estimates
##' at each time point is given by Poisson process method (by default)
##' or Lawless and Nadeau method (LawLess and Nadeau, 1995).
##' The approximate confidence intervals are provided as well,
##' which are constructed based on the asymptotic normality
##' of logarithm of MCF (by default) or MCF itself directly.
##'
##' For \code{\link{rateReg-class}} object,
##' \code{mcf} estimates the baseline MCF and its confidence interval
##' at each time grid if argument \code{newdata} is not specified.
##' Otherwise, \code{mcf} estimates MCF and its confidence interval
##' for the given newdata based on Delta-method.
##'
##' @param object An object used to dispatch a method.
##' @param na.action A function that indicates what should the procedure do
##' if the data contains \code{NA}s.  The default is set by the
##' na.action setting of \code{\link[base]{options}}.
##' The "factory-fresh" default is \code{\link[stats]{na.omit}}.
##' Other possible values inlcude \code{\link[stats]{na.fail}},
##' \code{\link[stats]{na.exclude}}, and \code{\link[stats]{na.pass}}.
##' \code{help(na.fail)} for details.
##' @param level An optional numeric value
##' indicating the confidence level required. The default value is 0.95.
##' @param control An optional list to specify the time grid
##' where the MCF is estimated for \code{rateReg-class} object and other options
##' for the formula method.
##' The available named elements include
##' \itemize{
##'    \item \code{grid}: The time grid where MCF is estimated. A dense grid is
##'        suggested for further using the plot method.
##'    \item \code{length.out}: The length of grid points. The dafault value
##'        is 1,000.
##'    \item \code{from}: The starting point of grid. The default value is the
##'        left boundary knots (for \code{rateReg-class} object).
##'    \item \code{to}: The endpoint of grid. The default value is the right
##'        boundary knots (for \code{rateReg-class} object).
##'    \item \code{B}: The number of bootstrap replicates for using bootstrap
##'        method for variance estimates of sample MCF estimates. The default
##'        value is 1,000.
##'    \item \code{se.method}: The method used for SE estimates for bootstrap.
##'        The default method is \code{"sampleSE"}, which takes the
##'        sample SE of point estimates from bootstrap samples.
##'    \item \code{ci.method}: The method used for confidence interval for
##'        bootstrap. The default method is the normal method.
##' }
##' The option \code{length.out}, \code{from}, \code{to} will be ignored if
##' \code{grid} is specified directly. Otherwise, the grid will be generated
##' by function \code{\link[base]{seq.int}} with specified \code{from},
##' \code{to} and \code{length.out}.
##' @param ... Other arguments for future usage.
##' @return
##' \code{\link{sampleMcf-class}} or \code{\link{rateRegMcf-class}} object.
##' Their slots include
##' \itemize{
##'     \item \code{level}: Confidence level specified.
##'     \item \code{MCF}: Mean cumulative function at each time point.
##'     \item \code{multiGroup}: A logical value indicating whether MCF
##'         is estimated for different specified group.
##'     \item \code{newdata}: Given dataset used to estimate MCF.
##' }
##' For the meaning of other slots, see \code{\link{rateReg}}.
##' @references
##' Lawless, J. F. and Nadeau, C. (1995). Some Simple Robust Methods for the
##' Analysis of Recurrent Events. \emph{Technometrics}, 37, 158--168.
##'
##' Nelson, W. B. (2003). \emph{Recurrent events data analysis for product
##' repairs, disease recurrences, and other applications} (Vol. 10). SIAM.
##'
##' @seealso
##' \code{\link{rateReg}} for model fitting;
##' \code{\link{plot-method}} for plotting MCF.
##' @examples
##' library(reda)
##'
##' ### Example 1. valve-seat data
##' valveMcf <- mcf(Survr(ID, Days, No.) ~ 1, data = valveSeats)
##'
##' ## plot sample MCF
##' plot(valveMcf, conf.int = TRUE, mark.time = TRUE) + ggplot2::xlab("Days")
##'
##' ### Example 2. sample simulated data
##' simuMcf <- mcf(Survr(ID, time, event) ~ group + gender,
##'                data = simuDat, ID %in% 1 : 50, logConfInt = FALSE)
##'
##' ## plot sample MCF
##' plot(simuMcf, conf.int = TRUE, lty = 1 : 4,
##'      legendName = "Treatment & Gender")
##'
##' ## For estimated MCF from a fitted model,
##' ## see examples given in function rateReg.
##' @export
setGeneric(name = "mcf",
           def = function(object, ...) {
               standardGeneric("mcf")
           })
