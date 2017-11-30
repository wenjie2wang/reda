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
##' An S4 class generic function that returns the mean cumulative function (MCF)
##' estimates from a fitted model or returns the nonparametric MCF estimates
##' (also called the Nelson-Aalen estimator) from the sample data.
##'
##' For \code{formula} object with \code{\link{Survr}} object as response, the
##' covariate specified at the right hand side of the formula should be either
##' \code{1} or any "linear" conbination of categorical variable in the data.
##' The former computes the overall sample MCF.  The latter computes the sample
##' MCF for each level of the combination of the categorical variable(s)
##' specified, respectively.  The sample MCF is also called Nelson-Aalen
##' nonparametric estimator (Nelson 2003) and computed on each time point from
##' sample data.  The point estimate of sample MCF at each time point does not
##' assume any particular underlying model. The variance estimates at each time
##' point is computed following the Lawless and Nadeau method (LawLess and
##' Nadeau 1995), the Poisson process method, or the bootstrap methods.  The
##' approximate confidence intervals are provided as well, which are constructed
##' based on the asymptotic normality of the MCF itself (by default) or the
##' logarithm of MCF.
##'
##' For \code{rateReg} object, \code{mcf} estimates the baseline
##' MCF and its confidence interval at each time grid if argument \code{newdata}
##' is not specified.  Otherwise, \code{mcf} estimates MCF and its confidence
##' interval for the given newdata based on Delta-method.
##'
##' @param object An object used to dispatch a method.
##' @param na.action A function that indicates what should the procedure do if
##'     the data contains \code{NA}s.  The default is set by the na.action
##'     setting of \code{\link[base]{options}}.  The "factory-fresh" default is
##'     \code{\link[stats]{na.omit}}.  Other possible values inlcude
##'     \code{\link[stats]{na.fail}}, \code{\link[stats]{na.exclude}}, and
##'     \code{\link[stats]{na.pass}}.  \code{help(na.fail)} for details.
##' @param level An optional numeric value indicating the confidence level
##'     required. The default value is 0.95.
##' @param control An optional named list specifying other options.  For
##'     \code{rateReg} object, it can be used to specify the time grid where the
##'     MCF is estimated. The available named elements are given as follows:
##'     \itemize{
##'         \item \code{grid}: The time grid where MCF is estimated. A dense
##'             grid is suggested for further using the plot method.
##'         \item \code{length.out}: The length of grid points. The dafault
##'             value is 1,000.
##'         \item \code{from}: The starting point of grid. The default value is
##'             the left boundary knots (for \code{rateReg} object).
##'         \item \code{to}: The endpoint of grid. The default value is the
##'             right boundary knots (for \code{rateReg} object).
##'     }
##'     The option \code{length.out}, \code{from}, \code{to} will be ignored if
##'     \code{grid} is specified directly. Otherwise, the grid will be generated
##'     by function \code{\link[base]{seq.int}} with specified \code{from},
##'     \code{to} and \code{length.out}.
##'
##'     For formula method, the available named elements are given as follows:
##'     \itemize{
##'       \item \code{B}: The number of bootstrap replicates for using
##'             bootstrap method for variance estimates of sample MCF estimates.
##'             The default value is 200.
##'         \item \code{se.method}: The method used for SE estimates for
##'             bootstrap. The available methods include \code{"sample.se"}
##'             (the default) and \code{"normality"}. The former takes the
##'             sample SE of point estimates from bootstrap samples; The latter
##'             estimates SE based on interquantile and normality assumption.
##'         \item \code{ci.method}: The method used for confidence interval (CI)
##'             for bootstrap. The available options include \code{"normality"}
##'             (the default) and \code{"percentile"}. The former estimates the
##'             CI based on SE estimates and normality assumption; The latter
##'             takes percentiles of the bootstrap estimates.
##'       \item \code{keep.data}: A logical value specifying whether to keep the
##'             processed data in the output object. If \code{TRUE}
##'             (the default), the processed data will be kept in the output and
##'             available for later usage. Otherwise, an empty data
##'             frame object will be returned in the \code{data} slot.
##'             \code{FALSE} may be set when the memory consumption is
##'             of concern and we only need MCF estimates. For example, the
##'             function \code{mcfDiff} and \code{mcfDiff.test} will not be
##'             applicable for the \code{mcf.formula} object with an empty
##'             \code{data} slot.
##'       \item \code{verbose}: A logical value. The default value is
##'             \code{TRUE}. If \code{FALSE}, possible data checking messages
##'              (not including warnings or errors) will be suppressed.
##'     }
##'
##' @param ... Other arguments for future usage.
##'
##' @return
##' A \code{mcf.formula} or \code{mcf.rateReg} object.
##'
##' A brief description of the slots of a \code{mcf.formula} object is given as
##' follows:
##' \itemize{
##'     \item \code{formula}: Model Formula.
##'     \item \code{data}: Processed data based on the model formula or an
##'         empty data frame if \code{keep.data} is set to be \code{FALSE}.
##'     \item \code{MCF}: A data frame containing estimates for sample MCF.
##'     \item \code{origin}: Time origins.
##'     \item \code{multiGroup}: A logical value indicating whether MCF
##'         is estimated for different groups respectively.
##'     \item \code{logConfInt}: A logical value indicating whether the
##'         variance estimates are based on the normality of logarithm of
##'         the MCF estimates.
##'     \item \code{level}: Confidence level specified.
##' }
##'
##' Most slots of a \code{mcf.rateReg} object are inherited from the input
##' \code{rateReg} object. A brief description of other slots is given as
##' follows:
##' \itemize{
##'     \item \code{newdata}: Given dataset used to estimate MCF.
##'     \item \code{MCF}: A data frame containing MCF estimates.
##'     \item \code{level}: Confidence level specified.
##'     \item \code{na.action}: The way handling missing values.
##'     \item \code{control}: The control list.
##'     \item \code{multiGroup}: A logical value indicating whether MCF
##'         is estimated for different groups respectively.
##' }
##'
##' @references
##'
##' Lawless, J. F. and Nadeau, C. (1995). Some Simple Robust Methods for the
##' Analysis of Recurrent Events. \emph{Technometrics}, 37, 158--168.
##'
##' Nelson, W. B. (2003). \emph{Recurrent Events Data Analysis for Product
##' Repairs, Disease Recurrences, and Other Applications} (Vol. 10). SIAM.
##'
##' @seealso
##' \code{\link{rateReg}} for model fitting;
##' \code{\link{mcfDiff}} for comparing two-sample MCFs.
##' \code{\link{plot-method}} for plotting MCF.
##' @examples
##' library(reda)
##'
##' ### sample MCF
##' ## Example 1. valve-seat data
##' ## the default variance estimates by Lawless and Nadeau (1995) method
##' valveMcf0 <- mcf(Survr(ID, Days, No.) ~ 1, data = valveSeats)
##' plot(valveMcf0, conf.int = TRUE, mark.time = TRUE, addOrigin = TRUE) +
##'     ggplot2::xlab("Days") + ggplot2::theme_bw()
##'
##' ## variance estimates following Poisson process model
##' valveMcf1 <- mcf(Survr(ID, Days, No.) ~ 1,
##'                  data = valveSeats, variance = "Poisson")
##' ## variance estimates by bootstrap method (with 1,000 bootstrap samples)
##' valveMcf2 <- mcf(Survr(ID, Days, No.) ~ 1,
##'                  data = valveSeats, variance = "bootstrap",
##'                  control = list(B = 2e2))
##'
##' ## comparing the variance estimates from different methods
##' library(ggplot2)
##' ciDat <- rbind(cbind(valveMcf0@MCF, Method = "Lawless & Nadeau"),
##'                cbind(valveMcf1@MCF, Method = "Poisson"),
##'                cbind(valveMcf2@MCF, Method = "Bootstrap"))
##' ggplot(ciDat, aes(x = time, y = se)) +
##'     geom_step(aes(color = Method, linetype = Method)) +
##'     xlab("Days") + ylab("SE estimates") + theme_bw()
##'
##' ## comparing the confidence interval estimates from different methods
##' ggplot(ciDat, aes(x = time)) +
##'     geom_step(aes(y = MCF)) +
##'     geom_step(aes(y = lower, color = Method, linetype = Method)) +
##'     geom_step(aes(y = upper, color = Method, linetype = Method)) +
##'     xlab("Days") + ylab("Confidence intervals") + theme_bw()
##'
##'
##' ## Example 2. the simulated data
##' simuMcf <- mcf(Survr(ID, time, event) ~ group + gender,
##'                data = simuDat, ID %in% 1 : 50)
##' plot(simuMcf, conf.int = TRUE, lty = 1 : 4,
##'      legendName = "Treatment & Gender")
##'
##' ### estimate MCF difference between two groups
##' ## one sample MCF object of two groups
##' mcf0 <- mcf(Survr(ID, time, event) ~ group, data = simuDat)
##' ## two-sample pseudo-score tests
##' mcfDiff.test(mcf0)
##' ## difference estimates over time
##' mcf0_diff <- mcfDiff(mcf0, testVariance = "none")
##' plot(mcf0_diff)
##'
##' ## or explicitly ask for the difference of two sample MCF
##' mcf1 <- mcf(Survr(ID, time, event) ~ 1, data = simuDat,
##'             subset = group %in% "Contr")
##' mcf2 <- mcf(Survr(ID, time, event) ~ 1, data = simuDat,
##'             subset = group %in% "Treat")
##' ## perform two-sample tests and estimate difference at the same time
##' mcf12_diff1 <- mcfDiff(mcf1, mcf2)
##' mcf12_diff2 <- mcf1 - mcf2   # or equivalently using the `-` method
##' stopifnot(all.equal(mcf12_diff1, mcf12_diff2))
##' mcf12_diff1
##' plot(mcf12_diff1)
##'
##' ### For estimated MCF from a fitted model,
##' ### see examples given in function rateReg.
##' @export
setGeneric(name = "mcf",
           def = function(object, ...) {
               standardGeneric("mcf")
           })
