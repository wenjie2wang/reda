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


### collation after class.R and mcf-generic.R
##' @include class.R
NULL


##' Estimate Difference Between Two MCF curves
##'
##' This function estimates MCF difference between two groups. Both point
##' estimates and confidence intervals are returned.
##'
##'
##'
##' @param mcf1 A \code{sampleMcf} object representing the MCF for one or two
##'     groups.
##' @param mcf2 An optional second \code{sampleMcf} object.
##' @param level An optional numeric value indicating the confidence level
##'     required. The default value is 0.95.
##' @param ... Other arguments for future usage.
##'
##' @return
##' A \code{mcfDiff} object that contains the following slots:
##' \itemize{
##'     \item \code{MCF}: Estimated Mean cumulative function Difference
##'         at each time point.
##'     \item \code{origin}: Time origins of the two groups.
##'     \item \code{variance}: The method used for variance estimates.
##'     \item \code{logConfInt}: A logical value indicating whether normality
##'         is assumed for \code{log(MCF)} instead of MCF itself.
##'         For \code{mcfDiff} object, it is always \code{FALSE}.
##'     \item \code{level}: Confidence level specified.
##' }
##'
##' @references
##'
##' Doganaksoy, N., & Nelson, W. (1998). A Method to Compare Two Samples of
##' Recurrence Data. \emph{Lifetime Data Analysis}, 4(1), 51--63.
##'
##' @examples
##' ## See examples given for function mcf.
##' @importFrom stats qnorm stepfun
##' @export
mcfDiff <- function(mcf1, mcf2 = NULL, level = 0.95, ...)
{
    ## quick checks
    if (! is.sampleMcf(mcf1) ||
        (! is.sampleMcf(mcf2) && ! is.null(mcf2))) {
        stop(wrapMessages(
            "'mcf1' must be 'sampleMcf' object and",
            "'mcf2' must be eithor 'sampleMcf' object or 'NULL'."
        ))
    }
    if (! isNumOne(level) || level <= 0 || level >= 1)
        stop("Confidence level must be between 0 and 1.")

    ## simple version
    ## if mcf1 contains only one level, mcf2 cannot be null
    ## if mcf1 contains more than one level, mcf2 will be ignored
    ## if mcf2 contains more than one level, error

    ## hard-coded const
    mcfCols <- c("time", "numRisk", "instRate",
                 "MCF", "se", "lower", "upper")
    nColMcf <- length(mcfCols)
    mcfColInd <- seq_len(nColMcf)

    if (mcf1@multiGroup) {
        nLevel <- length(mcf1@origin)
        if (nLevel > 2)
            stop("There were more than two groups.")
        uniLevs <- names(mcf1@origin)
        if (diff(mcf1@origin) != 0)
            warning("Time origins of two groups were not the same.")
        paste_ <- function(...) paste(..., sep = "_")
        datLevs <- sapply(seq_len(NROW(mcf1@MCF)), function (i) {
            paste_(mcf1@MCF[i, - mcfColInd])
        })
        ## first group
        mcfDat1 <- base::subset(mcf1@MCF, subset = datLevs %in% uniLevs[1L],
                                select = mcfCols)
        ## second group
        mcfDat2 <- base::subset(mcf1@MCF, subset = datLevs %in% uniLevs[2L],
                                select = mcfCols)
        ## origin
        originVec <- mcf1@origin
    } else {
        if (is.null(mcf2))
            stop(wrapMessages(
                "The object 'mcf2' cannot be missing",
                "if the object 'mcf1' only contains MCF for one group."
            ))
        if (mcf2@multiGroup)
            stop(wrapMessages(
                "The object 'mcf2' should contain MCF for only one group."
            ))
        if (mcf1@variance != mcf2@variance)
            warning(wrapMessages(
                "The method used for variance estimates were not consistent",
                "between the two 'sampleMcf' objects."
            ))
        ## first group
        mcfDat1 <- mcf1@MCF[, mcfCols]
        ## second group
        mcfDat2 <- mcf2@MCF[, mcfCols]
        ## origin
        originVec <- c(mcf1@origin, mcf2@origin)
    }
    ## first group
    mcfFun1 <- with(mcfDat1, stats::stepfun(time, c(0, MCF)))
    seFun1 <- with(mcfDat1, stats::stepfun(time, c(0, se)))
    ## second group
    mcfFun2 <- with(mcfDat2, stats::stepfun(time, c(0, MCF)))
    seFun2 <- with(mcfDat2, stats::stepfun(time, c(0, se)))
    ## unique time points
    uniTime <- sort(unique(c(mcfDat1$time, mcfDat2$time)))
    ## point difference
    mcfDiffVec <- mcfFun1(uniTime) - mcfFun2(uniTime)
    ## se estimates
    seVec <- sqrt((seFun1(uniTime)) ^ 2 + (seFun2(uniTime)) ^ 2)
    criVec <- stats::qnorm(0.5 + level / 2) * seVec
    mcfDiffDat <- data.frame(time = uniTime,
                             MCF = mcfDiffVec,
                             se = seVec,
                             lower = mcfDiffVec - criVec,
                             upper = mcfDiffVec + criVec)
    ## prepare for output
    methods::new("mcfDiff",
                 MCF = mcfDiffDat,
                 origin = originVec,
                 variance = mcf1@variance,
                 logConfInt = FALSE,
                 level = level)
}
