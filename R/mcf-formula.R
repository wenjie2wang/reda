##
## R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
## Copyright (C) 2015-2020
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
##'     method, \code{"bootstrap"} for bootstrap method, \code{"CSV"} for
##'     variance estimates of the corresponding cumulative sample mean function
##'     (CSM) by the cumulative sample variance method (Cook and Lawless, 2007),
##'     and \code{"none"} for no variance estimates.  Partial matching on the
##'     names is allowed.
##' @param logConfInt A logical value. If \code{FALSE} (the default), the
##'     confidence interval are constructed based on the normality of the MCF
##'     estimates. Otherwise, the confidence intervals of given level are
##'     constucted based on the normality of logarithm of the MCF estimates.
##' @param adjustRiskset A logical value indicating whether to adjust the size
##'     of risk-set.  If \code{TRUE} by default, the size of risk-set will be
##'     adjusted based on at-risk indicators and Nelson-Aalen estimator will be
##'     returned.  Otherwise, the cumulative sample mean (CSM) function given by
##'     Cook and Lawless (2007) will be returned without adjustment on size of
##'     risk-set.
##'
##' @aliases mcf,formula-method
##'
##' @importFrom stats na.exclude na.fail na.omit na.pass qnorm quantile sd
##' @export
setMethod(
    f = "mcf", signature = "formula",
    definition = function(object, data, subset, na.action,
                          variance = c("LawlessNadeau", "Poisson",
                                       "bootstrap", "CSV", "none"),
                          logConfInt = FALSE, adjustRiskset = TRUE,
                          level = 0.95, control = list(), ...)
    {
        ## specify the type of variance
        variance <- match.arg(variance)
        ## basic check on inputs
        if (! isLogicalOne(logConfInt, error_na = TRUE))
            stop(wrapMessages(
                "The argument 'logConfint' has to be",
                "either 'TRUE' or 'FALSE'."
            ), call. = FALSE)
        if (! isNumOne(level, error_na = TRUE) || level <= 0 || level >= 1)
            stop("Confidence level must be between 0 and 1.", call. = FALSE)

        ## get the data and take care of the possible subset
        if (missing(data))
            data <- environment(object)
        if (! missing(subset)) {
            sSubset <- substitute(subset)
            subIdx <- eval(sSubset, data, parent.frame())
            if (! is.logical(subIdx))
                stop("'subset' must be logical.", call. = FALSE)
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

        ## check response constructed from Recur
        resp <- stats::model.extract(mf, "response")
        if (! (is.Recur(resp) || is.Survr(resp)) )
            stop("Response in formula must be an 'Recur' object.",
                 call. = FALSE)
        ## check covariate in formula
        if (! NCOL(mm))
            stop("'1' or covariates must be specified in formula.",
                 call. = FALSE)
        Terms <- stats::terms(object)
        ord <- attr(Terms, "order")

        if (length(ord) & any(ord != 1))
            stop("Interaction term is not supported for this function.",
                 call. = FALSE)

        ## update control list
        control <- do.call(mcf_formula_control, control)
        point_method <- as.integer(adjustRiskset)
        var_method <- match(
            variance,
            c("LawlessNadeau", "Poisson", "bootstrap", "CSV")
        )
        ## no variance estimates otherwise
        if (is.na(var_method)) {
            var_method <- 0L
        }
        ci_method <- if (variance == "bootstrap" &&
                         control$ci.method == "percentile") {
                         3
                     } else if (logConfInt) {
                         2
                     } else {
                         1
                     }
        ci_level <- level
        var_bootstrap_method <- match(control$se.method,
                                      c("sample.se", "normality"))
        var_bootstrap_B <- control$B

        ## for possible missing values in covaraites
        if (length(na.action <- attr(mf, "na.action"))) {
            ## update if there is missing value removed
            attr(resp, "ID") <- attr(resp, "ID")[- na.action]
            ## check data for possible error caused by removal of missing values
            if (control$verbose)
                message("Observations with missing value in covariates ",
                        "are removed.\nChecking the new dataset again...\n",
                        appendLF = FALSE)
            if (is.Recur(resp)) {
                resp <- check_Recur(resp)
            } else {
                resp <- check_Survr(resp, check = TRUE)
            }
            if (control$verbose)
                message("Done!")
        }

        ## number of covariates
        nBeta <- length(mf) - 1L

        ## data processed
        dat <- as.data.frame(mf[[1L]])
        ## for compatibility of Survr
        if (! is.Recur(resp)) {
            ## create time1 for Survr object
            dat <- with(dat, Recur(time = time, id = ID,
                                   event = event, origin = origin))
            dat <- as.data.frame(dat)
        }
        dat <- dat[, c("time1", "time2", "id", "event", "terminal")]

        ## add covariates if any covariate is specified
        if (nBeta) {
            dat1 <- data.frame(lapply(mf[- 1L], factor))
            ## add covariates' names
            colnames(dat1) <- names(mf)[- 1L]
            dat <- cbind(dat, dat1)
        } else {
            ## no covariates
            outDat <- sMcf(dat,
                           point_method = point_method,
                           var_method = var_method,
                           ci_method = ci_method,
                           ci_level = ci_level,
                           var_bootstrap_method = var_bootstrap_method,
                           var_bootstrap_B = var_bootstrap_B,
                           groupLabel = NULL)
            ## remove all censoring rows? probably no for the plot method
            ## outDat <- base::subset(outDat, event == 1)
            rownames(outDat) <- NULL
            ## revert subject ID
            dat$id <- attr(resp, "ID")

            return(
                new("mcf.formula",
                    formula = object,
                    data = if (control$keep.data) dat else data.frame(),
                    MCF = outDat,
                    origin = min(dat$time1),
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
            warning("There is only one level in the covariates.",
                    call. = FALSE)

        ## initialize origin for each level...
        originVec <- rep(NA, num_levels)
        ## loop part
        outList <- lapply(seq(num_levels), function(i) {
            subDat <- dat[datLevs %in% levs[i], ]
            ## ...compute origin for each level
            originVec[i] <<- min(subDat$time1)
            oneSmcfDat <- sMcf(subDat,
                               point_method = point_method,
                               var_method = var_method,
                               ci_method = ci_method,
                               ci_level = ci_level,
                               var_bootstrap_method = var_bootstrap_method,
                               var_bootstrap_B = var_bootstrap_B,
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
            dat$id <- attr(resp, "ID")
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
## a warpper function for the underlying Rcpp routine
sMcf <- function(inpDat, point_method, var_method, ci_method, ci_level,
                 var_bootstrap_method, var_bootstrap_B, groupLabel = NULL)
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
    ## call cpp routine
    res_list <- cpp_np_mcf(
        time1 = inpDat$time1,
        time2 = inpDat$time2,
        id = inpDat$id,
        event = inpDat$event,
        point_method = point_method,
        var_method = var_method,
        ci_method = ci_method,
        ci_level = ci_level,
        var_bootstrap_method = var_bootstrap_method,
        var_bootstrap_B = var_bootstrap_B
    )
    ## if no variance estimate, set ci and se to NA
    if (var_method == 0L) {
        nout <- length(res_list$jump_time)
        res_list$se_cum_rate <- res_list$lower_cum_rate <-
            res_list$upper_cum_rate <- rep(NA_real_, nout)
    }
    ## return
    data.frame(
        time = res_list$jump_time,
        numRisk = res_list$riskset_size,
        instRate = res_list$inst_rate,
        MCF = res_list$cum_rate,
        se = res_list$se_cum_rate,
        lower = res_list$lower_cum_rate,
        upper = res_list$upper_cum_rate
    )
}

## function for preparing control list
mcf_formula_control <- function(B = 2e2,
                                se.method = c("sample.se", "normality"),
                                ci.method = c("normality", "percentile"),
                                keep.data = TRUE,
                                verbose = TRUE,
                                ...)
{
    if (! isNumOne(B, error_na = TRUE) || B < 2)
        stop(wrapMessages(
            "The number of bootstarp samples 'B'",
            "should be a postive integer > 1."
        ), call. = FALSE)
    if (B < 30) {
        warning(wrapMessages(
            "A larger number of bootstarp samples is suggested."
        ), call. = FALSE)
    }
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
