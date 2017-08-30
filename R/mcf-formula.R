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


## collation after class.R and mcf-generic.R
##' @include class.R
##' @include mcf-generic.R
NULL


##' @describeIn mcf Sample MCF from data.
##'
##' @param data An optional data frame, list or environment containing
##' the variables in the model.  If not found in data, the variables are taken
##' from \code{environment(formula)}, usually the environment from which
##' the function is called.
##' @param subset An optional vector specifying a subset of observations
##' to be used in the fitting process.
##' @param variance An optional character specifying the variance estimator.
##' The available options are "Poisson" (default) for Poisson process method,
##' and "LawlessNadeau" for Lawless and Nadeau (1995) method. (A simple example
##' is available at Reliawiki, 2012.)
##' Partial matching on the names is allowed.
##' @param logConfInt An optional logical value. If \code{TRUE} (default),
##' the confidence interval of given level will be constucted based on the
##' normality of logarithm of the MCF function. (Otherwise, the confidence
##' interval will be constructed based on the normality of MCF function.)
##' @aliases mcf,formula-method
##' @importFrom stats na.fail na.omit na.exclude na.pass qnorm
##' @export
setMethod(
    f = "mcf", signature = "formula",
    definition = function(object, data, subset, na.action,
                          variance = c("Poisson", "LawlessNadeau"),
                          logConfInt = TRUE, level = 0.95, ...)
    {
        ## specify the type of variance
        variance <- match.arg(variance)
        ## basic check on inputs
        logConfInt <- logConfInt[1L]
        if (! is.logical(logConfInt))
            stop("'logConfint' has to be either 'TRUE' or 'FALSE'.")
        ## check on level specified
        level <- level[1L]
        if (level <= 0 || level >= 1)
            stop("Confidence level must be between 0 and 1.")

        if (missing(data))
            data <- environment(object)
        ## take care of subset individual for possible non-numeric ID
        if (! missing(subset)) {
            sSubset <- substitute(subset)
            subIdx <- eval(sSubset, data, parent.frame())
            if (! is.logical(subIdx))
                stop("'subset' must be logical.")
            subIdx <- subIdx & ! is.na(subIdx)
            data <- data[subIdx, ]
        }

        ## Prepare data: ID, time, event, origin ~ 1 or X
        mcall <- match.call(expand.dots = FALSE)
        mcall[[1]] <- as.name("mcf")
        names(mcall) <- sub("object", "formula", names(mcall))
        mfnames <- c("formula", "data", "na.action")
        mfind <- match(mfnames, names(mcall), nomatch = 0L)
        mcall$formula <- eval(object)
        ## re-define data
        mcall$data <- substitute(data)
        ## match mcall
        mmcall <- match(c("formula", "data", "na.action"), names(mcall), 0L)
        mcall$na.action <- eval(substitute(alist(na.action)))[[1]]
        mcall <- mcall[c(1L, mmcall)]
        ## drop unused levels in factors
        mcall$drop.unused.levels <- TRUE
        mcall[[1L]] <- quote(stats::model.frame)
        mf <- eval(mcall, parent.frame())
        mt <- attr(mf, "terms")
        mm <- stats::model.matrix(object, data = mf)

        ## check response constructed from Survr
        resp <- stats::model.extract(mf, "response")
        if (! inherits(resp, "Survr"))
            stop("Response in formula must be a survival recurrent object.")
        ## check covariate in formula
        if (! NCOL(mm))
            stop("'1' or covariates must be specified in formula.")
        Terms <- stats::terms(object)
        ord <- attr(Terms, "order")

        if (length(ord) & any(ord != 1))
            stop("Interaction term is not supported for this function.")

        ## for possible missing values in covaraites
        if (length(na.action <- attr(mf, "na.action"))) {
            ## update if there is missing value removed
            attr(resp, "ord") <- order(resp[, "ID"], resp[, "time"])
            attr(resp, "ID") <- attr(resp, "ID")[- na.action]
            check_Survr(resp, check = TRUE)
        }

        ## number of covariates
        nBeta <- length(mf) - 1L

        ## data processed
        dat <- as.data.frame(mf[[1L]])
        if (nBeta) {
            ## covariates' names
            covar_names <- names(mf)[- 1L]
            dat1 <- data.frame(lapply(mf[- 1L], factor))
            dat <- cbind(dat, dat1)
        } else {
            covar_names <- NULL
        }
        const_colnames <- c("ID", "time", "event", "origin")
        colnames(dat) <- c(const_colnames, covar_names)
        ## revert subject ID
        dat$ID <- attr(resp, "ID")

        ## output: na.action
        na.action <- if (is.null(attr(mm, "na.action"))) {
                         options("na.action")[[1]]
                     } else {
                         paste0("na.", class(attr(mm, "na.action")))
                     }

        ## if no covariates specified
        if (! nBeta) {
            outDat <- sMcf0(inpDat = dat, variance = variance,
                            logConfInt = logConfInt, level = level)
            ## remove all censoring rows? probably no for the plot method
            ## outDat <- base::subset(outDat, event == 1)
            rownames(outDat) <- NULL
            out <- new("sampleMcf",
                       formula = object,
                       MCF = outDat,
                       multiGroup = FALSE,
                       na.action = na.action,
                       variance = variance,
                       logConfInt = logConfInt,
                       level = level)
            return(out)
        }

        ## else at least one covariate are specified
        ## get the levels for each covaraite in form of grid
        xGrid <- unique(dat1)
        levs <- apply(xGrid, 1, paste, collapse = "_")
        datLevs <- apply(dat1, 1, paste, collapse = "_")

        ## number of levels
        num_levels <- NROW(xGrid)
        if (num_levels == 1L)
            warning("There is only one level in the covariates.")

        mcf_colnames <- c(const_colnames, "numRisk", "instRate",
                          "MCF", "se", "lower", "upper")
        const_ncol <- length(mcf_colnames)
        outDat <- data.frame(matrix(NA, nrow = nrow(mm),
                                    ncol = const_ncol + nBeta))
        colIdx <- seq_len(const_ncol)
        for (i in seq(num_levels)) {
            subDat <- dat[datLevs %in% levs[i], ]
            rowLen <- nrow(subDat)
            rowInd <- if(i == 1L) {
                          seq_len(rowLen)
                      } else {
                          seq(from = rowInd[length(rowInd)] + 1L, by = 1L,
                              length.out = rowLen)
                      }
            outDat[rowInd, colIdx] <- sMcf0(subDat, variance,
                                            logConfInt, level)
            outDat[rowInd, - colIdx] <- xGrid[i, ]
        }
        colnames(outDat) <- c(mcf_colnames, covar_names)

        ## factorize covariates
        outCol <- ncol(outDat)
        for (j in seq_len(nBeta)) {
            outDat[, outCol + 1 - j] <-
                factor(levels(dat1[[nBeta + 1 - j]])[outDat[, outCol + 1 - j]])
        }
        rownames(outDat) <- NULL

        out <- methods::new("sampleMcf",
                            formula = object,
                            MCF = outDat,
                            multiGroup = TRUE,
                            na.action = na.action,
                            variance = variance,
                            logConfInt = logConfInt,
                            level = level)
        ## return
        out
    })

### internal function ==========================================================
## compute sample MCF for all zero origin
sMcf0 <- function(inpDat, variance, logConfInt, level)
{
    ## if time ties, put event time before censoring time
    inpDat <- inpDat[with(inpDat, base::order(time, 1 - event)), seq_len(4L)]

    num_pat <- length(unique(inpDat$ID))
    num_at_risk <- num_pat - cumsum(inpDat$event != 1)
    increment <- ifelse(inpDat$event, 1 / num_at_risk, 0)

    ## index of unique event and censoring time
    ## indx <- ! duplicated(inpDat$time, fromLast = TRUE)

    ## sample mcf at each time point
    smcf <- cumsum(increment)
    inc2 <- increment ^ 2

    ## sample variance estimator
    if (variance == "Poisson") {
        var_smcf <- cumsum(inc2)
        se_smcf <- sqrt(var_smcf)
    } else if (variance == "LawlessNadeau") {
        inc12 <- (1 - increment) ^ 2
        incre <- inc2 * (inc12 + (num_at_risk - 1) * inc2)
        incre_var <- ifelse(inpDat$event < 1, 0, incre)
        var_smcf <- cumsum(incre_var)
        se_smcf <- sqrt(var_smcf)
    } else {
        stop("Invalid 'variance' specified.")
    }

    ## Confidence interval for log(MCF) or MCF
    if (logConfInt) {
        criVal <- ifelse(smcf > 0,
                         stats::qnorm(0.5 + level / 2) * se_smcf / smcf,
                         0)
        wtonexp <- exp(criVal)
        upper <- smcf * wtonexp
        lower <- smcf / wtonexp
    } else {
        criVal <- stats::qnorm(0.5 + level / 2) * se_smcf
        upper <- smcf + criVal
        lower <- pmax(smcf - criVal, 0)
    }

    ## return
    data.frame(inpDat,
               numRisk = num_at_risk,
               instRate = increment,
               MCF = smcf,
               se = se_smcf,
               lower = lower,
               upper = upper)
}
