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
##' @include mcf-generic.R
NULL


##' @describeIn mcf Sample MCF from data.
##'
##' @param data A data frame, list or environment containing the variables in
##'     the model.  If not found in data, the variables are taken from
##'     \code{environment(formula)}, usually the environment from which the
##'     function is called.
##' @param subset An optional vector specifying a subset of observations to be
##'     used in the fitting process.
##' @param variance A character specifying the method for variance estimates.
##'     The available options are \code{"LawlessNadeau"} (the default) for
##'     Lawless and Nadeau (1995) method, \code{"Poisson"} for Poisson process
##'     method, and \code{"bootstrap"} for bootstrap method.  Partial matching
##'     on the names is allowed.
##' @param logConfInt A logical value. If \code{FALSE} (the default), the
##'     confidence interval are constructed based on the normality of the MCF
##'     estimates. Otherwise, the confidence intervals of given level are
##'     constucted based on the normality of logarithm of the MCF estimates.
##'
##' @aliases mcf,formula-method
##'
##' @importFrom stats na.exclude na.fail na.omit na.pass qnorm quantile sd
##' @export
setMethod(
    f = "mcf", signature = "formula",
    definition = function(object, data, subset, na.action,
                          variance = c("LawlessNadeau", "Poisson", "bootstrap"),
                          logConfInt = FALSE,
                          level = 0.95,
                          control = list(), ...)
    {
        ## specify the type of variance
        variance <- match.arg(variance)
        ## basic check on inputs
        if (! isLogicalOne(logConfInt, error_na = TRUE))
            stop(wrapMessages(
                "The argument 'logConfint' has to be",
                "either 'TRUE' or 'FALSE'."
            ))
        if (! isNumOne(level, error_na = TRUE) || level <= 0 || level >= 1)
            stop("Confidence level must be between 0 and 1.")

        ## get the data and take care of the possible subset
        if (missing(data))
            data <- environment(object)
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
        ## mcall$data <- substitute(data)
        mcall$data <- data
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
        if (! is.Survr(resp))
            stop("Response in formula must be a survival recurrent object.")
        ## check covariate in formula
        if (! NCOL(mm))
            stop("'1' or covariates must be specified in formula.")
        Terms <- stats::terms(object)
        ord <- attr(Terms, "order")

        if (length(ord) & any(ord != 1))
            stop("Interaction term is not supported for this function.")

        ## update control list
        control <- do.call(mcf_formula_control, control)

        ## for possible missing values in covaraites
        if (length(na.action <- attr(mf, "na.action"))) {
            ## update if there is missing value removed
            attr(resp, "ord") <- order(resp[, "ID"], resp[, "time"])
            attr(resp, "ID") <- attr(resp, "ID")[- na.action]
            ## check data for possible error caused by removal of missing values
            if (control$verbose)
                message("Observations with missing value in covariates ",
                        "are removed.\nChecking the new dataset again...\n",
                        appendLF = FALSE)
            check_Survr(resp, check = TRUE)
            if (control$verbose)
                message("Done!")
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
        colnames(dat) <- c("ID", "time", "event", "origin", covar_names)

        ## if no covariates specified
        if (! nBeta) {
            outDat <- sMcf(dat,
                           variance = variance,
                           logConfInt = logConfInt,
                           level = level,
                           control = control,
                           groupLabel = NULL)
            ## remove all censoring rows? probably no for the plot method
            ## outDat <- base::subset(outDat, event == 1)
            rownames(outDat) <- NULL
            ## revert subject ID
            dat$ID <- attr(resp, "ID")

            return(
                new("mcf.formula",
                    formula = object,
                    data = if (control$keep.data) dat else data.frame(),
                    MCF = outDat,
                    origin = min(dat$origin),
                    multiGroup = FALSE,
                    variance = variance,
                    logConfInt = logConfInt,
                    level = level)
            )
        }

        ## else at least one covariate are specified
        ## get the levels for each covaraite in form of grid
        xGrid <- unique(dat1[do.call(order, as.list(dat1)), , drop = FALSE])
        levs <- apply(xGrid, 1L, paste, collapse = ":")
        datLevs <- apply(dat1, 1L, paste, collapse = ":")
        ## number of levels
        num_levels <- NROW(xGrid)
        if (num_levels == 1L)
            warning("There is only one level in the covariates.")

        ## initialize origin for each level...
        originVec <- rep(NA, num_levels)
        ## loop part
        outList <- lapply(seq(num_levels), function(i) {
            subDat <- dat[datLevs %in% levs[i], ]
            ## ...compute origin for each level
            originVec[i] <<- min(subDat$origin)
            oneSmcfDat <- sMcf(dat = subDat,
                               variance = variance,
                               logConfInt = logConfInt,
                               level = level,
                               control = control,
                               groupLabel = levs[i])
            do.call(data.frame, c(as.list(oneSmcfDat),
                                  as.list(xGrid[i, , drop = FALSE])))
        })
        outDat <- do.call(rbind, outList)
        ## name origin vector
        names(originVec) <- levs
        rownames(outDat) <- NULL

        ## whether to keep data in output
        if (control$keep.data) {
            ## revert subject ID
            dat$ID <- attr(resp, "ID")
        } else {
            dat <- data.frame()
        }

        ## return
        methods::new("mcf.formula",
                     formula = object,
                     data = dat,
                     MCF = outDat,
                     origin = originVec,
                     multiGroup = TRUE,
                     variance = variance,
                     logConfInt = logConfInt,
                     level = level)
    }
)


### internal function ==========================================================
## simple wrapper function for sMcf_point and addVar_sMcf
sMcf <- function(dat, variance, logConfInt, level, control, groupLabel)
{
    sMcfDat <- sMcf_point(dat, groupLabel = groupLabel)
    addVar_sMcf(dat, sMcfDat, variance, logConfInt, level, control)
}

## sample MCF point estimates
sMcf_point <- function(inpDat, groupLabel = NULL)
{
    ## throw out warning if no any event
    if (all(! inpDat$event)) {
        warning(wrapMessages(
            "No event found",
            if (is.null(groupLabel))
                "overall."
            else
                sprintf("in the %s group.", groupLabel)
        ), call. = FALSE)
    }

    ## different origins?
    origin_idx <- length(unique(inpDat$origin)) > 1
    if (origin_idx) {
        ## already sorted by ID, time, and (1 - event) by Survr
        dup_idx <- duplicated(inpDat$ID)
        originVec <- inpDat$origin[! dup_idx]
        rightVec <- inpDat[! inpDat$event, "time"]
    }
    ## if time ties, put event time before censoring time
    inpDat <- inpDat[with(inpDat, order(time, - event)), seq_len(4L)]
    num_at_risk <- if (origin_idx) {
                       ## for non-zero origin
                       sapply(inpDat$time, function(xTime) {
                           sum(xTime >= originVec & xTime <= rightVec)
                       })
                   } else {
                       ## for all equal origin
                       length(unique(inpDat$ID)) - cumsum(inpDat$event <= 0)
                   }
    increment <- with(inpDat, ifelse(event > 0, event / num_at_risk, 0))

    ## index of unique event and censoring time
    uni_first_idx <- ! duplicated(inpDat$time, fromLast = FALSE)
    uni_last_idx <- ! duplicated(inpDat$time, fromLast = TRUE)

    ## sample mcf at each unique time point
    smcf <- cumsum(increment)[uni_last_idx]

    ## return
    data.frame(time = inpDat$time[uni_first_idx],
               numRisk = num_at_risk[uni_first_idx],
               instRate = diff(c(0, smcf)),
               MCF = smcf,
               se = NA,
               lower = NA,
               upper = NA)
}

## add variance estimates for sample MCF point estimates
addVar_sMcf <- function(dat, sMcfDat, variance, logConfInt, level, control)
{
    ## point estimates of mcf
    smcf <- sMcfDat$MCF

    ## sample variance estimator
    if (variance == "LawlessNadeau") {
        var_smcf <- var_lawlessNadeau(dat, sMcfDat)
        se_smcf <- sqrt(var_smcf)
    } else if (variance == "Poisson") {
        varVec <- with(sMcfDat, ifelse(numRisk > 0, instRate / numRisk, 0))
        var_smcf <- cumsum(varVec)
        se_smcf <- sqrt(var_smcf)
    } else if (variance == "bootstrap") {
        sMcfBootMat <- sMcf_boot(dat, sMcfDat, B = control$B)
        se_smcf <- switch(
            control$se.method,
            "sample.se" = seBoot_sampleSE(sMcfBootMat),
            "normality" = seBoot_normality(sMcfBootMat, upperQuan = 0.75)
        )
    }

    ## confidence interval
    if (variance == "bootstrap" && control$ci.method != "normality") {
        ciMat <- switch(
            control$ci.method,
            percentile = ciBoot_percentile(sMcfBootMat, level, logConfInt)
        )
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
            lower <- smcf - criVal
        }
    }
    ## update mcf data
    sMcfDat$se <- se_smcf
    sMcfDat$lower <- lower
    sMcfDat$upper <- upper
    ## return
    sMcfDat
}

## function computing variance estimates by Lawless and Nadeau (1995)
var_lawlessNadeau <- function(dat, sMcfDat)
{
    ## following notations in Lawless and Nadeau (1995)
    ## at unique time points...
    tj <- sMcfDat$time
    ## ... number of processes under risk
    delta_tj <- sMcfDat$numRisk
    ## ... and events or costs
    m_tj <- sMcfDat$instRate
    ## function for each process
    var_comp <- function(subDat, tj, delta_tj, m_tj) {
        ## compute delta_i(tj) as delta_i_tj
        origin <- subDat$origin[1L]
        cenTime <- with(subDat, time[event <= 0])
        delta_i_tj <- as.numeric(tj >= origin & tj <= cenTime)
        ## compute n_i(t_j) as n_i_tj
        n_i_tj <- rep(0, length(tj))
        names(n_i_tj) <- as.character(tj)
        ## if multiple events exists at the same time
        tj_i <- unique(subDat$time)
        n_i_tj[as.character(tj_i)] <- if (any(duplicated(subDat$time)))
                                          with(subDat, tapply(event, time, sum))
                                      else
                                          subDat$event
        ## apply formula (2.3)--(2.5)
        res_ij <- ifelse(delta_tj > 0,
                         delta_i_tj / delta_tj * (n_i_tj - m_tj),
                         0)
        ## return
        cumsum(res_ij) ^ 2
    }
    ## apply var_comp to each process
    varList <- by(dat[, c("ID", "time", "event", "origin")], dat$ID,
                  var_comp, tj = tj, delta_tj = delta_tj, m_tj = m_tj,
                  simplify = FALSE)
    ## just in case dimension was dropped
    if (nrow(sMcfDat) > 1)
        colSums(do.call(rbind, varList))
    else
        sum(do.call(c, varList))
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
    bootList <- rowIndList[subInd]
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
    bootDat$ID <- as.numeric(factor(bootDat$ID))
    ## apply sMcf_point to bootDat
    mcfDat <- sMcf_point(bootDat)
    stepMcf <- with(mcfDat, stepfun(time, c(0, MCF)))
    ## variance estimate for inst rate?
    stepMcf(grid)
}

## estimate SE based on normality assumption
seBoot_normality <- function(sMcfBootMat, upperQuan = 0.75)
{
    if (! isNumOne(upperQuan, error_na = TRUE) ||
        upperQuan >= 1 || upperQuan <= 0.5)
        stop(wrapMessages(
            "The upper quantile ('upperQuan') has to be between 0.5 and 1."
        ), call. = FALSE)
    quanDiff <- qnorm(upperQuan) - qnorm(1 - upperQuan)
    tmp <- apply(sMcfBootMat, 1L, quan_fun, level = 2 * upperQuan - 1)
    (tmp[2L, ] - tmp[1L, ]) / quanDiff
}

## estimate SE by sample SE
seBoot_sampleSE <- function(sMcfBootMat)
{
    apply(sMcfBootMat, 1L, stats::sd)
}

## estimate CI by normality assumption
ciBoot_percentile <- function(sMcfBootMat, level, logConfInt)
{
    if (logConfInt) {
        exp(apply(log(sMcfBootMat), 1L, quan_fun))
    } else {
        apply(sMcfBootMat, 1L, quan_fun)
    }
}

## estimate CI by BCa method (accelation and bias-correction)
## ciBoot_BCa <- function(sMcfBootMat, level) {
##     ## TODO
## }

## extract required quantiles
quan_fun <- function(x, level = 0.95)
{
    half_level <- level / 2 * c(- 1, 1)
    quan2 <- 0.5 + half_level
    stats::quantile(x, probs = quan2)
}


## function for preparing control list
mcf_formula_control <- function(B = 2e2,
                                se.method = c("sample.se", "normality"),
                                ci.method = c("normality", "percentile"),
                                keep.data = TRUE,
                                verbose = TRUE,
                                ...)
{
    if (! isNumOne(B, error_na = TRUE) || B <= 1)
        stop(wrapMessages(
            "The number of bootstarp sample, 'B' should be a numeric number."
        ), call. = FALSE)
    if (B < 30)
        warning(wrapMessages(
            "A larger number of bootstarp samples is suggested."
        ), call. = FALSE)
    se.method <- match.arg(se.method)
    ci.method <- match.arg(ci.method)
    if (! isLogicalOne(keep.data, error_na = TRUE))
        stop("The option 'keep.data' should be either 'TRUE' or 'FALSE'.",
             call. = FALSE)
    ## return
    list(
        B = B,
        se.method = se.method,
        ci.method = ci.method,
        keep.data = keep.data,
        verbose = verbose
    )
}
