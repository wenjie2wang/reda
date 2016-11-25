################################################################################
##
##   R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
##   Copyright (C) 2015-2016
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
##' @param ... Other arguments for future usage.
##' @param na.action A function that indicates what should the procedure do
##' if the data contains \code{NA}s.  The default is set by the
##' na.action setting of \code{\link[base]{options}}.
##' The "factory-fresh" default is \code{\link[stats]{na.omit}}.
##' Other possible values inlcude \code{\link{na.fail}},
##' \code{\link{na.exclude}}, and \code{\link{na.pass}}.
##' \code{help(na.fail)} for details.
##' @param level An optional numeric value
##' indicating the confidence level required. The default value is 0.95.
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
##' Nelson, W. B. (2003). \emph{Recurrent events data analysis for product
##' repairs, disease recurrences, and other applications} (Vol. 10). SIAM.
##'
##' Lawless, J. F. and Nadeau, C. (1995). Some Simple Robust Methods for the
##' Analysis of Recurrent Events. \emph{Technometrics}, 37, 158--168.
##'
##' ReliaWiki. (2012, March 19). Recurrent Event Data Analysis.
##' Retrieved November 23, 2015, from
##' \url{http://reliawiki.org/index.php/Recurrent_Event_Data_Analysis}
##' @seealso
##' \code{\link{rateReg}} for model fitting;
##' \code{\link{plot}} for plotting MCF.
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
##' ## create customized levels in legend
##' levs <- with(simuDat, expand.grid(levels(group), levels(gender)))
##' levs <- do.call(paste, c(as.list(levs), sep = " & "))
##'
##' ## plot sample MCF
##' plot(simuMcf, conf.int = TRUE, lty = 1 : 4,
##'      legendName = "Treatment & Gender", legendLevels = levs)
##'
##' ## For estimated MCF from a fitted model,
##' ## see examples given in function rateReg.
##' @export
setGeneric(name = "mcf",
           def = function(object, ...) {
               standardGeneric("mcf")
           })


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
##' @importFrom stats na.fail na.omit na.exclude na.pass
##' @export
setMethod(
    f = "mcf", signature = "formula",
    definition = function(object, data, subset, na.action,
                          variance = c("Poisson", "LawlessNadeau"),
                          logConfInt = TRUE, level = 0.95, ...) {
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
                stop("'subset' must be logical")
            subIdx <- subIdx & ! is.na(subIdx)
            data <- data[subIdx, ]
        }

        ## Prepare data: ID, time, event ~ 1 or X
        mcall <- match.call(expand.dots = FALSE)
        mcall[[1]] <- as.name("mcf")
        names(mcall) <- sub(pattern = "object",
                            replacement = "formula", names(mcall))
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
            stop("Interaction terms are not valid for this function.")

        ## for possible missing values in covaraites
        if (length(na.action <- attr(mf, "na.action"))) {
            ## update if there is missing value removed
            attr(resp, "ord") <- order(resp[, "ID"], resp[, "time"])
            attr(resp, "ID_") <- attr(resp, "ID_")[- na.action]
            check_Survr(resp, check = TRUE)
        }

        ## Sorted data by ID, time, and event
        ord <- attr(resp, "ord")
        ## number of covariates excluding intercept
        nBeta <- ncol(mm) - 1L

        ## data matrix processed
        if (nBeta) {
            xMat <- mm[ord, - 1L, drop = FALSE]
            ## covariates' names
            covar_names <- colnames(mm)[- 1L]
        } else {
            xMat <- covar_names <- NULL
        }
        dat <- as.data.frame(cbind(resp[ord, ], xMat))
        colnames(dat) <- c("ID", "time", "event", covar_names)
        ## revert subject ID
        dat$ID <- attr(resp, "ID_")[ord]

        ## ouput: variable names in formula
        varNames <- as.character(object[[2L]])[- 1L]

        ## output: na.action
        na.action <- if (is.null(attr(mm, "na.action"))) {
                         options("na.action")[[1]]
                     } else {
                         paste0("na.", class(attr(mm, "na.action")))
                     }

        ## if no covariates specified
        if (! nBeta) {
            outDat <- sMcf(inpDat = dat, variance = variance,
                           logConfInt = logConfInt, level = level)
            colnames(outDat)[seq_len(3)] <- varNames
            ## remove all censoring rows
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
        xGrid <- unique(xMat)
        levs <- apply(xGrid, 1, paste, collapse = "_")
        datLevs <- apply(xMat, 1, paste, collapse = "_")
        ## number of levels
        num_levels <- NROW(xGrid)
        if (num_levels == 1L)
            warning("The covariate has only one level.")

        outDat <- data.frame(matrix(NA, nrow = nrow(mm), ncol = 7L + nBeta))
        for (i in seq(num_levels)) {
            subDat <- dat[datLevs %in% levs[i], ]
            rowLen <- nrow(subDat)
            rowInd <- if(i == 1L) {
                          seq_len(rowLen)
                      } else {
                          seq(from = rowInd[length(rowInd)] + 1L, by = 1L,
                              length.out = rowLen)
                      }
            outDat[rowInd, seq_len(7)] <- sMcf(subDat, variance,
                                               logConfInt, level)
            outDat[rowInd, - seq_len(7)] <- rep(xGrid[i, ], each = rowLen)
        }
        colnames(outDat) <- c(varNames, "MCF", "se",
                              "lower", "upper", covar_names)
        ## remove all censoring rows? not now for plot
        ## outDat <- base::subset(outDat, event == 1)
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


##' @describeIn mcf Estimated MCF from a fitted model.
##'
##' @param newdata An optional data frame. If specified, the data frame should
##' have the same column names as the covariate names appearing in the formula
##' of original fitting.
##' @param groupName An optional length-one charactor vector to specify the
##' name for grouping each unique row in \code{newdata}, such as "gender"
##' for "male" and "female". The default value is "group".
##' @param groupLevels An optional charactor vector to specify the levels for
##' each unique row in \code{newdata}, such as "treatment" and "control".
##' The default values are capital letters starting from "A".
##' @param control An optional list to specify the time grid
##' where the MCF are estimated.
##' The availble elements of the control list include
##' \code{grid}, \code{length.out}, \code{from} and \code{to}.
##' The time grid can be directly specified via element \code{grid}.
##' A dense time grid is suggested.
##' Element \code{length.out} represents the length of grid points.
##' The dafault value is 1,000.
##' Element \code{from} means the starting point of grid with default 0.
##' Element \code{to} represnts the endpoint of grid
##' with the right boundary knot as default.
##' When \code{grid} is missing, the grid will be generated
##' by \code{seq} (from package \pkg{base})
##' with arguments \code{from}, \code{to} and \code{length.out}.
##' @aliases mcf,rateReg-method
##' @importFrom stats na.fail na.omit na.exclude na.pass
##' @importFrom splines2 ibs iSpline
##' @export
setMethod(
    f = "mcf", signature = "rateReg",
    definition = function(object, newdata, groupName, groupLevels,
                          level = 0.95, na.action, control = list(), ...) {

        ## extract model fitting information from object
        beta <- object@estimates$beta[, 1L]
        alpha <- object@estimates$alpha[, 1L]
        fcovnames <- as.character(object@call[[2L]][[3L]])
        covnames <- fcovnames[fcovnames != "+"]
        nBeta <- length(beta)
        knots <- object@spline$knots
        degree <- object@spline$degree
        Boundary.knots <- object@spline$Boundary.knots

        ## control
        controlist <- c(control, list(Boundary.knots_ = Boundary.knots))
        control <- do.call("rateReg_mcf_control", controlist)
        gridTime <- control$grid
        n_xx <- control$length.out

        ## baseline rate basis matrix
        iMat <- splines2::ibs(x = gridTime, knots = knots,
                              degree = degree, intercept = TRUE,
                              Boundary.knots = Boundary.knots)
        bMat <- attr(iMat, "bsMat")

        ## estimated MCF
        estMcf <- as.vector(iMat %*% alpha)

        ## covariance matrix of beta and alpha
        n_par <- nrow(object@fisher)
        covInd <- seq_len(n_par)[- (nBeta + 1L)]
        covMat <- solve(object@fisher)[covInd, covInd]

        ## nonsense, just to suppress Note from R CMD check --as-cran
        `(Intercept)` <- NULL

        ## about newdata
        tt <- stats::terms(object@formula)
        Terms <- stats::delete.response(tt)
        if (missing("na.action"))
            na.action <- options("na.action")[[1L]]

        if (missing(newdata)) {
            X <- matrix(rep(0, nBeta), nrow = 1)
            colnames(X) <- covnames
            rownames(X) <- "1"
        } else {
            mf <- stats::model.frame(Terms, newdata,
                                     na.action = na.action,
                                     xlev = object@xlevels)
            if (is.null(attr(mf, "na.action"))) {
                na.action <- options("na.action")[[1L]]
            } else {
                na.action <- paste0("na.", class(attr(mf, "na.action")))
            }
            X <- stats::model.matrix(Terms, mf, contrasts.arg =
                                                    object@contrasts$constracts)
            ## remove intercept and deplicated rows
            X <- unique(base::subset(X, select = - `(Intercept)`))
            if (ncol(X) != nBeta)
                stop(paste("The number of input covariates does not",
                           "match with 'rateReg' object."))
        }

        ## mcf for possible multigroups
        ndesign <- nrow(X)
        multiGroup <- ndesign > 1

        coveff <- as.numeric(exp(X %*% beta))
        outDat <- matrix(NA, ncol = 4, nrow = ndesign * n_xx)

        for (i in seq(ndesign)) {
            ## Delta-method
            gradMat <- gradDelta(Xi = X[i, ], iMat, estMcf, coveffi = coveff[i])
            varTime <- apply(gradMat, 1, function (gradVec, covMat) {
                crossprod(gradVec, covMat) %*% gradVec
            }, covMat = covMat)
            confBand <- stats::qnorm((1 + level) / 2) * sqrt(varTime)
            mcf_i <- estMcf * coveff[i]
            lower <- pmax(0, mcf_i - confBand)
            upper <- mcf_i + confBand
            ind <- seq(from = n_xx * (i - 1) + 1, to = n_xx * i, by = 1)
            outDat[ind, 1L] <- gridTime
            outDat[ind, 2L] <- mcf_i
            outDat[ind, 3L] <- lower
            outDat[ind, 4L] <- upper
        }
        outDat <- as.data.frame(outDat)
        colnames(outDat) <- c("time", "MCF", "lower", "upper")
        if (multiGroup) {
            if (missing(groupLevels))
                groupLevels <- LETTERS[seq(ndesign)]
            if (missing(groupName))
                groupName <- "group"
            tempcol <- factor(rep(groupLevels[seq(ndesign)], each = n_xx))
            outDat <- cbind(outDat, tempcol)
            colnames(outDat)[5L] <- groupName
        }
        ## output
        out <- new("rateRegMcf",
                   call = object@call,
                   formula = object@formula,
                   spline = object@spline$spline,
                   knots = knots, degree = degree,
                   Boundary.knots = Boundary.knots,
                   newdata = X, MCF = outDat, level = level,
                   na.action = na.action, control = control,
                   multiGroup = multiGroup)
        ## return
        out
    })


### internal function ==========================================================
## compute sample MCF
sMcf <- function(inpDat, variance, logConfInt, level) {
    ## if time ties, put event time before censoring time
    inpDat <- inpDat[with(inpDat, base::order(time, 1 - event)), seq_len(3)]

    num_pat <- length(unique(inpDat$ID))
    num_at_risk <- num_pat - cumsum(inpDat$event == 0)
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
        incre_var <- ifelse(inpDat$event == 0, 0, incre)
        var_smcf <- cumsum(incre_var)
        se_smcf <- sqrt(var_smcf)
    } else {
        stop("Invalid 'variance' specified.")
    }

    ## Confidence interval for log(MCF) or MCF
    if (logConfInt) {
        criVal <- stats::qnorm(0.5 + level / 2) * se_smcf / smcf
        wtonexp <- exp(criVal)
        upper <- smcf * wtonexp
        lower <- smcf / wtonexp
    } else {
        criVal <- stats::qnorm(0.5 + level / 2) * se_smcf
        upper <- smcf + criVal
        lower <- smcf - criVal
    }

    ## return
    data.frame(inpDat, MCF = smcf, se = se_smcf,
               lower = lower, upper = upper)
}


## Delta-method, compute the gradient on beta and alpha
gradDelta <- function(Xi, iMat, estMcf, coveffi) {
    gradBeta <- estMcf %o% Xi * coveffi
    gradAlpha <- iMat * coveffi
    ## retrun
    cbind(gradBeta, gradAlpha)
}


## control function
rateReg_mcf_control <- function(grid, length.out = 1e3, from, to, ...,
                                Boundary.knots_) {
    ## controls for function MCF with signiture rateReg
    from <- if (missing(from)) Boundary.knots_[1]
    to  <- if (missing(to)) Boundary.knots_[2]
    if (! missing(grid)) {
        if (! is.numeric(grid) || is.unsorted(grid))
            stop("'grid' specified must be an increasing numeric vector.")
        length.out <- length(grid)
        from <- min(grid)
        to <- max(grid)
    } else {
        grid <- seq(from = from, to = to, length.out = length.out)
    }
    if (min(grid) < Boundary.knots_[1] || max(grid) > Boundary.knots_[2])
        stop("'grid' must be within the coverage of boundary knots.")

    ## return
    list(grid = grid, length.out = length.out, from = from, to = to)
}
