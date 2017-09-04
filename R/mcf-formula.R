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
##' @param data An optional data frame, list or environment containing the
##'     variables in the model.  If not found in data, the variables are taken
##'     from \code{environment(formula)}, usually the environment from which the
##'     function is called.
##' @param subset An optional vector specifying a subset of observations to be
##'     used in the fitting process.
##' @param variance An optional character specifying the method for variance
##'     estimates.  The available options are \code{"Poisson"} (default) for
##'     Poisson process method, \code{"LawlessNadeau"} for Lawless and Nadeau
##'     (1995) method, and \code{"bootstrap"} for bootstrapping method. (A
##'     simple example is available at Reliawiki, 2012.)  Partial matching on
##'     the names is allowed.
##' @param logConfInt An optional logical value. If \code{TRUE} (default), the
##'     confidence interval of given level will be constucted based on the
##'     normality of logarithm of the MCF function. (Otherwise, the confidence
##'     interval will be constructed based on the normality of MCF function.)
##'
##' @aliases mcf,formula-method
##' @importFrom stats na.fail na.omit na.exclude na.pass qnorm quantile sd
##' @export
setMethod(
    f = "mcf", signature = "formula",
    definition = function(object, data, subset, na.action,
                          variance = c("Poisson", "LawlessNadeau", "bootstrap"),
                          logConfInt = TRUE, level = 0.95,
                          control = list(), ...)
    {
        ## specify the type of variance
        variance <- match.arg(variance)
        ## basic check on inputs
        if (length(logConfInt) > 1 || ! is.logical(logConfInt))
            stop("'logConfint' has to be either 'TRUE' or 'FALSE'.")
        ## check on level specified
        if (length(level) > 1 || ! is.numeric(level) ||
            level <= 0 || level >= 1)
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
            outDat <- sMcf(dat,
                           variance = variance,
                           logConfInt = logConfInt,
                           level = level,
                           control = control)
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
            outDat[rowInd, colIdx] <- sMcf(subDat, variance,
                                           logConfInt, level, control)
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
## simple wrapper function for sMcf_point and addVar_sMcf
sMcf <- function(dat, variance, logConfInt, level, control) {
    sMcfDat <- sMcf_point(dat)
    addVar_sMcf(dat, sMcfDat, variance, logConfInt, level, control)
}

## sample MCF point estimates
sMcf_point <- function(inpDat)
{
    ## throw out warning if no any event
    if (all(! inpDat$event))
        warning("No event found.")

    ## origin all the same?
    origin_idx <- length(unique(inpDat$origin)) > 1
    if (origin_idx) {
        ## already sorted by ID, time, and (1 - event) by Survr
        dup_idx <- duplicated(inpDat$ID)
        originVec <- inpDat$origin[! dup_idx]
        rightVec <- inpDat[! inpDat$event, "time"]
    }
    ## if time ties, put event time before censoring time
    inpDat <- inpDat[with(inpDat, order(time, 1 - event)), seq_len(4L)]
    num_at_risk <- if (origin_idx) {
                       ## for non-zero origin
                       sapply(inpDat$time, function(xTime) {
                           sum(xTime >= originVec & xTime <= rightVec)
                       })
                   } else {
                       ## for all zero origin
                       length(unique(inpDat$ID)) - cumsum(inpDat$event != 1)
                   }
    increment <- ifelse(inpDat$event, 1 / num_at_risk, 0)

    ## index of unique event and censoring time
    ## indx <- ! duplicated(inpDat$time, fromLast = TRUE)

    ## sample mcf at each time point
    smcf <- cumsum(increment)

    ## return
    data.frame(inpDat,
               numRisk = num_at_risk,
               instRate = increment,
               MCF = smcf,
               se = NA,
               lower = NA,
               upper = NA)
}

## add variance estimates for sample MCF point estimates
addVar_sMcf <- function(dat, sMcfDat, variance, logConfInt, level, control)
{
    increment <- sMcfDat$instRate
    num_at_risk <- sMcfDat$numRisk
    smcf <- sMcfDat$MCF
    inc2 <- increment ^ 2

    ## sample variance estimator
    if (identical(variance, "Poisson")) {
        var_smcf <- cumsum(inc2)
        se_smcf <- sqrt(var_smcf)
    } else if (identical(variance, "LawlessNadeau")) {
        inc12 <- (1 - increment) ^ 2
        incre <- inc2 * (inc12 + (num_at_risk - 1) * inc2)
        incre_var <- ifelse(sMcfDat$event < 1, 0, incre)
        var_smcf <- cumsum(incre_var)
        se_smcf <- sqrt(var_smcf)
    } else if (identical(variance, "bootstrap")) {
        control <- do.call(addVar_sMcf_control, control)
        sMcfBootMat <- sMcf_boot(dat, sMcfDat, B = control$B)
        se_smcf <- seBoot_normal(sMcfBootMat, upperQuan = 0.75)
    } else
        stop("Unknown variance type.")

    ## confidence interval
    if (identical(variance, "bootstrap") && control$ci.method != "normal") {
        ciMat <- switch(control$ci.method,
                        percentile = ciBoot_percentile(sMcfBootMat, level))
        lower <- ciMat[1L, ]
        upper <- ciMat[2L, ]
    } else {
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
    }

    ## update mcf data
    sMcfDat$se <- se_smcf
    sMcfDat$lower <- lower
    sMcfDat$upper <- upper
    ## return
    sMcfDat
}

## function preparing control
addVar_sMcf_control <- function(B = 1e3,
                                se.method = c("sampleSE", "normal"),
                                ci.method = c("normal", "percentile"))
{
    se.method <- match.arg(se.method)
    ci.method <- match.arg(ci.method)
    list(B = B, se.method = se.method, ci.method = ci.method)
}

## warpper function computing sample MCFs for bootstrap sample generated
sMcf_boot <- function(dat, sMcfDat, B)
{
    ## set up some constant variable invariant to bootstrap replicates
    idTab <- table(dat$ID)
    uID <- names(idTab)
    numSub <- length(idTab)
    rowIndList <- tapply(seq_len(nrow(dat)), dat$ID, c)
    ## reset row.names for later usage by sMcf_boot_one
    row.names(dat) <- NULL
    grid <- sMcfDat$time
    ## loop of sMcf_boot_one
    replicate(B, sMcf_boot_one(dat, idTab, uID, numSub, rowIndList, grid))
}

## generate bootstrap samples and compute its sample MCF
sMcf_boot_one <- function(dat, idTab, uID, numSub, rowIndList, grid, ...)
{
    subInd <- sample(numSub, replace = TRUE)
    subTab <- table(subInd)
    ## recover subject ID with zero count in bootstrap sample
    zeroIdx <- ! match(uID, names(subTab), nomatch = 0)
    zeroTab <- rep(0, sum(zeroIdx))
    names(zeroTab) <- uID[zeroIdx]
    compTab <- c(subTab, zeroTab)
    aID <- factor(names(compTab), levels = levels(dat$ID))
    compTab <- compTab[order(aID)]
    bootList <- rep(rowIndList, compTab)
    bootVec <- do.call(c, bootList)
    ## bootstrap sample data
    bootDat <- dat[bootVec, ]
    ## take advantage of updated rownames indicating duplicates
    rowNames <- row.names(bootDat)
    tmp_m <- regexpr("\\.[0-9]+", rowNames)
    dupNames <- rep("", nrow(bootDat))
    dupNames[tmp_m > 0] <- regmatches(rowNames, tmp_m)
    ## update subject ID in bootstrap sample
    bootDat$ID <- paste0(bootDat$ID, "_boot", dupNames)

    ## apply sMcf_point to bootDat
    mcfDat <- sMcf_point(bootDat)
    stepMcf <- with(mcfDat, stepfun(time, c(0, MCF)))

    ## variance estimate for inst rate?
    stepMcf(grid)
}

## estimate SE based on normality assumption
seBoot_normal <- function(sMcfBootMat, upperQuan = 0.75)
{
    if (! is.numeric(upperQuan) || upperQuan >= 1 || upperQuan <= 0.5)
        stop("The upper quantile ('upperQuan') has to be between 0.5 and 1.")
    quanDiff <- qnorm(upperQuan) - qnorm(1 - upperQuan)
    tmp <- apply(sMcfBootMat, 1L, quan_fun, level = 2 * upperQuan - 1)
    (tmp[2L, ] - tmp[1L, ]) / quanDiff
}

## estimate SE by sample SE
seBoot_sample <- function(sMcfBootMat) {
    apply(sMcfBootMat, 1L, stats::sd)
}

## estimate CI by percentile method
ciBoot_percentile <- function(sMcfBootMat, level = 0.95)
{
    apply(sMcfBootMat, 1L, quan_fun)
}

## estimate CI by BCa method (accelation and bias-correction)
ciBoot_BCa <- function(sMcfBootMat, level) {
    ## TODO
}

## extract required quantiles
quan_fun <- function(x, level = 0.95)
{
    half_level <- level / 2 * c(- 1, 1)
    quan2 <- 0.5 + half_level
    stats::quantile(x, probs = quan2)
}
