##
## R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
## Copyright (C) 2015-2019
##
## This file is part of the R package reda.
##
## The R package reda is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package reda is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##


### collation after class.R and mcf-generic.R
##' @include class.R
NULL


##' Comparing Two-Sample MCFs
##'
##' This function estimates the sample MCF difference between two groups. Both
##' the point estimates and the confidence intervals are computed (Lawless and
##' Nadeau 1995). The two-sample pseudo-score test proposed by Cook, Lawless,
##' and Nadeau (1996) is also performed by default.
##'
##' The function \code{mcfDiff} estimates the two-sample MCFs' difference and
##' internally calls function \code{mcfDiff.test} to perform the pseudo-score
##' tests by default. A \code{-} method is available as a simple wrapper for the
##' function \code{mcfDiff} for comparing two-sample MCFs from two
##' \code{mcf.formula} objects. For instance, suppose \code{mcf1} and
##' \code{mcf2} are \code{mcf.formula} objects, each of which represents the
##' sample MCF estimates for one group. The function call \code{mcf1 - mcf2} is
##' equivalent to \code{mcfDiff(mcf1, mcf2)}.
##'
##' The null hypothesis of the two-sample pseudo-score test is that there is no
##' difference between the two sample MCFs, while the alternative hypothesis
##' suggests a difference.  The test is based on a family of test statistics
##' proposed by Lawless and Nadeau (1995). The argument \code{testVariance}
##' specifies the method for computing the variance estimates of the test
##' statistics under different model assumption. See the document of argument
##' \code{testVariance} for all applicable options.  For the variance estimates
##' robust to departures from Poisson process assumption, both constant weight
##' and the linear weight function (with scaling) suggested in Cook, Lawless,
##' and Nadeau (1996) are implemented. The constant weight is powerful in cases
##' where the two MCFs are approximately proportional to each other. The linear
##' weight function is originally \code{a(u) = t - u}, where \code{u} represents
##' the time variable and \code{t} is the first time point when the risk set of
##' either group becomes empty. It is further scaled by \code{1 / t} for test
##' statistics invariant to the unit of measurement of the time variable.  The
##' linear weight function puts more emphasis on the difference at earily times
##' than later times and is more powerful for cases where the MCFs are no longer
##' proportional to each other, but not crossing. Also see Cook and Lawless
##' (2007, Section 3.7.5) for more details.
##'
##' @aliases mcfDiff
##'
##' @usage
##' mcfDiff(mcf1, mcf2 = NULL, level = 0.95, ...)
##'
##' @param mcf1 A \code{mcf.formula} object representing the MCF for one or two
##'     groups.
##' @param mcf2 An optional second \code{mcf.formula} object or \code{NULL}.
##' @param level A numeric value indicating the confidence level required. The
##'     default value is 0.95.
##' @param ... Other arguments passed to \code{mcfDiff.test}.
##'
##' @return
##'
##' The function \code{mcfDiff} returns a \code{mcfDiff} object (of S4 class)
##' that contains the following slots:
##' \itemize{
##' \item \code{call}: Function call.
##' \item \code{MCF}: Estimated Mean cumulative function Difference at each time
##' point.
##' \item \code{origin}: Time origins of the two groups.
##' \item \code{variance}: The method used for variance estimates.
##' \item \code{logConfInt}: A logical value indicating whether normality is
##' assumed for \code{log(MCF)} instead of MCF itself.  For \code{mcfDiff}
##' object, it is always \code{FALSE}.
##' \item \code{level}: Confidence level specified.
##' \item \code{test}: A \code{mcfDiff.test} object for the hypothesis test
##' results.
##' }
##'
##' The function \code{mcfDiff.test} returns a \code{mcfDiff.test} object (of S4
##' class) that contains the following slots:
##' \itemize{
##' \item \code{.Data}: A numeric matrix (of two rows and five columns) for
##' hypothesis testing results.
##' \item \code{testVariance}: A character string (or vector of length one)
##' indicating the method used for the variance estimates of the test statistic.
##' }
##'
##' @references
##'
##' Lawless, J. F., & Nadeau, C. (1995). Some Simple Robust Methods for the
##' Analysis of Recurrent Events. \emph{Technometrics}, 37(2), 158--168.
##'
##' Cook, R. J., Lawless, J. F., & Nadeau, C. (1996). Robust Tests for Treatment
##' Comparisons Based on Recurrent Event Responses. \emph{Biometrics}, 52(2),
##' 557--571.
##'
##' Cook, R. J., & Lawless, J. (2007). \emph{The Statistical Analysis of
##' Recurrent Events}.  Springer Science & Business Media.
##'
##' @examples
##' ## See examples given for function mcf.
##' @importFrom stats pchisq qnorm stepfun
##' @export
mcfDiff <- function(mcf1, mcf2 = NULL, level = 0.95, ...)
{
    ## record function call
    Call <- match.call()

    ## quick checks
    mcfDiff_check(mcf1, mcf2)
    if (! isNumOne(level) || is.na(level) || level <= 0 || level >= 1)
        stop("Confidence level must be between 0 and 1.", call. = FALSE)

    ## simple version
    ## if mcf1 contains only one level, mcf2 cannot be null
    ## if mcf1 contains more than one level, mcf2 will be ignored
    ## if mcf2 contains more than one level, error

    ## define some constants
    mcfCols <- c("time", "numRisk", "instRate",
                 "MCF", "se", "lower", "upper")
    nColMcf <- length(mcfCols)
    mcfColInd <- seq_len(nColMcf)

    if (mcf1@multiGroup) {
        nLevel <- length(mcf1@origin)
        if (nLevel > 2)
            stop("There were more than two groups in 'mcf1'.", call. = FALSE)
        if (! is.null(mcf2))
            warning("Only the 'mcf1' object was used.", call. = FALSE)
        uniLevs <- names(mcf1@origin)
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
        ## origins
        originVec <- mcf1@origin
        ## variance
        variance <- mcf1@variance
    } else {
        if (is.null(mcf2))
            stop(wrapMessages(
                "The object 'mcf2' cannot be missing",
                "if the object 'mcf1' only contains MCF for one group."
            ), call. = FALSE)
        if (mcf2@multiGroup)
            stop(wrapMessages(
                "The object 'mcf2' should contain MCF for only one group."
            ), call. = FALSE)
        if (mcf1@variance != mcf2@variance)
            warning(wrapMessages(
                "The methods used for variance estimates were not consistent",
                "between the two 'mcf.formula' objects!"
            ), call. = FALSE)
        ## first group
        mcfDat1 <- mcf1@MCF[, mcfCols]
        ## second group
        mcfDat2 <- mcf2@MCF[, mcfCols]
        ## origins
        originVec <- c(mcf1@origin, mcf2@origin)
        ## variance
        variance <- unique(c(mcf1@variance, mcf2@variance))
    }
    ## warning if the earliest time origins of two groups were different
    if (diff(originVec) != 0)
        warning(wrapMessages(
            "The earliest time origins of the two groups were not the same!"
        ), call. = FALSE)
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
    ## testing part
    testRes <- mcfDiff.test(mcf1, mcf2, ...)

    ## prepare for output
    methods::new("mcfDiff",
                 call = Call,
                 MCF = mcfDiffDat,
                 origin = originVec,
                 variance = variance,
                 logConfInt = FALSE,
                 level = level,
                 test = testRes)
}


##' @rdname mcfDiff
##' @aliases `-`,mcf.formula-method
##'
##' @param e1 The first \code{mcf.formula} object, \code{mcf1}.
##' @param e2 The second \code{mcf.formula} object, \code{mcf2}.
##'
##' @export
setMethod(
    f = "-",
    signature = c("mcf.formula", "mcf.formula"),
    definition = function(e1, e2)
    {
        out <- mcfDiff(e1, e2)
        ## substitute for function call
        out@call[[2L]] <- substitute(e1)
        out@call[[3L]] <- substitute(e2)
        out
    }
)


##' @rdname mcfDiff
##' @aliases mcfDiff.test
##'
##' @usage
##' mcfDiff.test(mcf1, mcf2 = NULL,
##'              testVariance = c("robust", "Poisson", "none"), ...)
##'
##' @param testVariance A character string specifying the method for computing
##'     the variance estimate for the pseudo-score test statistic proposed by
##'     Cook, Lawless, and Nadeau (1996). The applicable options include
##'     \code{"robust"} (default) for an estimate robust to departures from
##'     Poisson assumptions, \code{"Poisson"} for an estimate for Poisson
##'     process, and \code{"none"} for not performing any test (if only the
##'     difference estimates are of interest in \code{mcfDiff}).
##'
##' @importFrom stats stepfun
##' @export
mcfDiff.test <- function(mcf1, mcf2 = NULL,
                         testVariance = c("robust", "Poisson", "none"),
                         ...)
{
    ## reference: Cook, R. J., Lawless, J. F., & Nadeau, C. (1996)
    ## consider two types of weight functions a(u) for mcf difference
    ## constant and linear weight function
    ## and two variance estimates for the test statistic
    testVariance <- match.arg(testVariance)

    ## quick checks
    mcfDiff_check(mcf1, mcf2)

    if (testVariance == "none")
        return(methods::new("mcfDiff.test"))
    ## if no process data in mcf1
    if (nrow(mcf1@data) == 0L) {
        warning("No processed data is available from 'mcf1'.", call. = FALSE)
        return(methods::new("mcfDiff.test"))
    }

    ## define some constants
    mcfCols <- c("time", "numRisk", "instRate",
                 "MCF", "se", "lower", "upper")
    nColMcf <- length(mcfCols)
    mcfColInd <- seq_len(nColMcf)
    datCols <- c("ID", "time", "event", "origin")
    nColDat <- length(datCols)
    datColInd <- seq_len(nColDat)

    ## non-sense to pass R CMD check for "no visible binding"
    time <- NULL

    ## FIXME
    ## temp solution: convert output of Recur to Survr
    conv2Survr <- function(mcf_obj) {
        dat <- mcf_obj@data
        first_idx <- ! duplicated(dat$id)
        origin_vec <- rep(dat$time1[first_idx], table(factor(dat$id)))
        out <- data.frame(ID = dat$id, time = dat$time2,
                          event = dat$event, origin = origin_vec)
        dat$time1 <- dat$time2 <- dat$id <- dat$event <- dat$terminal <- NULL
        mcf_obj@data <- cbind(out, dat)
        mcf_obj
    }
    mcf1 <- conv2Survr(mcf1)
    if (! is.null(mcf2)) {
        mcf2 <- conv2Survr(mcf2)
    }

    if (mcf1@multiGroup) {
        nLevel <- length(mcf1@origin)
        if (nLevel > 2)
            stop("There were more than two groups in 'mcf1'.", call. = FALSE)
        uniLevs <- names(mcf1@origin)
        if (! is.null(mcf2))
            warning("Only the 'mcf1' object was used.", call. = FALSE)

        getLevs <- function(dat, colInd) {
            paste_ <- function(...) paste(..., sep = "_")
            sapply(seq_len(nrow(dat)), function (i, colInd) {
                paste_(dat[i, - colInd])
            }, colInd)
        }
        datLevs <- getLevs(mcf1@data, datColInd)
        mcfDatLevs <- getLevs(mcf1@MCF, mcfColInd)
        ## first group
        eventDat1 <- base::subset(mcf1@data, datLevs %in% uniLevs[1L])
        mcfDat1 <- base::subset(mcf1@MCF, mcfDatLevs %in% uniLevs[1L],
                                select = mcfCols)
        ## second group
        eventDat2 <- base::subset(mcf1@data, datLevs %in% uniLevs[2L])
        mcfDat2 <- base::subset(mcf1@MCF, mcfDatLevs %in% uniLevs[2L],
                                select = mcfCols)
        ## origins
        originVec <- mcf1@origin
    } else {
        if (is.null(mcf2))
            stop(wrapMessages(
                "The object 'mcf2' cannot be missing",
                "if the object 'mcf1' only contains MCF for one group."
            ), call. = FALSE)
        if (mcf2@multiGroup)
            stop(wrapMessages(
                "The object 'mcf2' should contain MCF for only one group."
            ), call. = FALSE)
        ## if no process data in mcf2
        if (nrow(mcf2@data) == 0L) {
            warning("No processed data is available from 'mcf2'.",
                    call. = FALSE)
            return(methods::new("mcfDiff.test"))
        }
        ## first group
        eventDat1 <- mcf1@data
        mcfDat1 <- mcf1@MCF[, mcfCols]
        ## second group
        eventDat2 <- mcf2@data
        mcfDat2 <- mcf2@MCF[, mcfCols]
        ## origins
        originVec <- c(mcf1@origin, mcf2@origin)
    }
    ## warning if the earliest time origins of two groups were different
    if (diff(originVec) != 0)
        warning(wrapMessages(
            "The earliest time origins of the two groups were not the same!"
        ), call. = FALSE)
    ## internal function
    ## Y_k.(u) for k = 1, 2
    y_k <- function(eventDat, newTimes) {
        endIdx <- ! duplicated(eventDat$ID, fromLast = TRUE)
        endTimeVec <- eventDat$time[endIdx]
        originVec <- eventDat$origin[endIdx]
        ## if origins are not all the same
        if (length(unique(originVec)) > 1) {
            sapply(newTimes, function(a) {
                sum(a >= originVec & a <= endTimeVec)
            })
        } else {
            sapply(newTimes, function(a) {
                sum(a <= endTimeVec)
            })
        }
    }
    ## dLambda_k(u) for k = 1, 2
    dLambda_k <- function(mcfDat, newTimes) {
        tmpDat <- base::subset(mcfDat, time %in% newTimes)
        res <- rep(0, length(newTimes))
        names(res) <- as.character(newTimes)
        res[as.character(tmpDat$time)] <- tmpDat$instRate
        res
    }

    ## unique times
    uniTimeVec <- sort(unique(c(mcfDat1$time, mcfDat2$time)))
    yVec1 <- y_k(eventDat1, uniTimeVec)
    yVec2 <- y_k(eventDat2, uniTimeVec)
    ## only need to compute for non-zero elements
    nonZeroIdx <- yVec1 > 0 & yVec2 > 0
    uniTimeVec <- uniTimeVec[nonZeroIdx]
    yVec1 <- yVec1[nonZeroIdx]
    yVec2 <- yVec2[nonZeroIdx]

    ## linear weight function
    ## using the same scaling with sas (hmmm... interesting)
    ## not documented by sas but figured out based on toy example
    ## similar but not sure whether exactly equivalant to the linear weight
    ## suggested in Cook. et.al. (1996).
    tau <- max(uniTimeVec)
    a_u <- (tau - uniTimeVec) / tau

    ## constant weight
    w_const <- (yVec1 * yVec2) / (yVec1 + yVec2)
    max_w_const <- max(w_const)
    w_linear <- a_u * w_const

    dLambda1 <- dLambda_k(mcfDat1, uniTimeVec)
    dLambda2 <- dLambda_k(mcfDat2, uniTimeVec)
    jumpDiff <- dLambda1 - dLambda2
    ## re-scale to avoid numerical issues
    max_jumpDiff <- max(abs(jumpDiff))
    if (max_jumpDiff > 0) {
        re_jumpDiff <- jumpDiff / max_jumpDiff
        reFactor <- exp(log(max_jumpDiff) + log(max_w_const))
    } else {
        re_jumpDiff <- 0
        reFactor <- 1
    }
    ## test statistics U(t)
    testU_const <- sum(w_const / max_w_const * re_jumpDiff) * reFactor
    testU_linear <- sum(w_linear / max_w_const * re_jumpDiff) * reFactor

    ## variance part
    if (testVariance == "Poisson") {
        varComp <- ifelse(yVec1 > 0, dLambda1 / yVec1, 0) +
            ifelse(yVec2 > 0, dLambda2 / yVec2, 0)
        varU_const <- sum(w_const ^ 2 * varComp)
        varU_linear <- sum(w_linear ^ 2 * varComp)
        varU <- c(varU_const, varU_linear)
    } else {
        ## function for each process
        var_comp <- function(subDat, yVec, dLambda,
                             uniTimeVec, w_const, w_linear,
                             max_w_const) {
            ## only compute for time points where w(u) > 0
            tmpDat <- base::subset(subDat, time %in% uniTimeVec)
            ## compute delta_i(tj) as delta_i_tj
            origin <- subDat$origin[1L]
            cenTime <- with(subDat, time[event <= 0])
            y_i_tj <- as.numeric(uniTimeVec >= origin &
                                 uniTimeVec <= cenTime)
            ## compute n_i(t_j) as n_i_tj
            n_i_tj <- rep(0, length(uniTimeVec))
            names(n_i_tj) <- as.character(uniTimeVec)
            ## if multiple events exists at the same time
            tj_i <- unique(tmpDat$time)
            n_i_tj[as.character(tj_i)] <-
                if (any(duplicated(tmpDat$time))) {
                    with(tmpDat, tapply(event, time, sum))
                } else {
                    tmpDat$event
                }
            ## apply formula (3.3) in Cook. et al. (1996)
            res_ij <- ifelse(yVec > 0,
                             y_i_tj / yVec * (n_i_tj - dLambda),
                             0)
            ## rescaling
            max_res_ij <- max(abs(res_ij))
            if (max_res_ij > 0) {
                re_res_ij <- res_ij / max_res_ij
                reFactor <- exp(log(max_res_ij) + log(max_w_const))
            } else {
                re_res_ij <- 0
                reFactor <- 1
            }
            res_const <- (w_const / max_w_const) * re_res_ij
            res_linear <- (w_linear / max_w_const) * re_res_ij
            ## return
            c(const = (sum(res_const) * reFactor) ^ 2 ,
              linear = (sum(res_linear) * reFactor) ^ 2)
        }
        ## apply var_comp to each process
        ## k = 1
        varList1 <- by(eventDat1[, c("ID", "time", "event", "origin")],
                       eventDat1$ID, FUN = var_comp,
                       yVec = yVec1, dLambda = dLambda1,
                       uniTimeVec = uniTimeVec,
                       w_const = w_const,
                       w_linear = w_linear,
                       max_w_const = max_w_const,
                       simplify = FALSE)
        ## k = 2
        varList2 <- by(eventDat2[, c("ID", "time", "event", "origin")],
                       eventDat2$ID, FUN = var_comp,
                       yVec = yVec2, dLambda = dLambda2,
                       uniTimeVec = uniTimeVec,
                       w_const = w_const,
                       w_linear = w_linear,
                       max_w_const = max_w_const,
                       simplify = FALSE)

        ## just in case dimension was dropped
        varU_1 <-
            if (nrow(mcfDat1) > 1) {
                colSums(do.call(rbind, varList1))
            } else {
                do.call(c, varList1)
            }
        varU_2 <-
            if (nrow(mcfDat2) > 1) {
                colSums(do.call(rbind, varList2))
            } else {
                do.call(c, varList2)
            }
        varU <- varU_1 + varU_2
    }

    ## summarize the results in a table
    outTab <- methods::new("mcfDiff.test")
    outTab[, 1L] <- c(testU_const, testU_linear)
    outTab[, 2L] <- varU
    outTab[, 3L] <- outTab[, 1L] ^ 2 / outTab[, 2L]
    outTab[, 4L] <- 1
    outTab[, 5L] <- stats::pchisq(outTab[, 3L],
                                  df = outTab[, 4L],
                                  lower.tail = FALSE)
    outTab@testVariance <- testVariance
    ## return
    outTab
}


### internal functions =========================================================
mcfDiff_check <- function(mcf1, mcf2) {
    ## quick check
    if (! is.mcf.formula(mcf1) ||
        (! is.mcf.formula(mcf2) && ! is.null(mcf2))) {
        stop(wrapMessages(
            "'mcf1' must be 'mcf.formula' object, and",
            "'mcf2' must be eithor 'mcf.formula' object or 'NULL'."
        ), call. = FALSE)
    }
}
