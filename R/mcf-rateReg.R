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


## collation after class.R and mcf-generic.R
##' @include class.R
##' @include mcf-generic.R
NULL


##' @describeIn mcf Estimated MCF from a fitted model.
##'
##' @param newdata An optional data frame. If specified, the data frame should
##'     have the same column names as the covariate names appearing in the
##'     formula of original fitting.
##' @param groupName An optional length-one charactor vector to specify the name
##'     for grouping each unique row in \code{newdata}, such as "gender" for
##'     "male" and "female". The default value is "Group".
##' @param groupLevels An optional charactor vector to specify the levels for
##'     each unique row in \code{newdata}, such as "treatment" and "control".
##'     The default values are \code{"Level"} followed by a numeric sequence
##'     with length of number of levels.
##'
##' @aliases mcf,rateReg-method
##'
##' @importFrom stats na.fail na.omit na.exclude na.pass
##' @importFrom splines2 ibs iSpline
##' @export
setMethod(
    f = "mcf", signature = "rateReg",
    definition = function(object, newdata, groupName, groupLevels,
                          level = 0.95, na.action, control = list(), ...)
    {
        ## extract model fitting information from object
        beta <- object@estimates$beta[, "coef"]
        alpha <- object@estimates$alpha[, "coef"]
        covnames <- names(beta)
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
        spline <- object@spline$spline
        if (identical(spline, "bSplines")) {
            iMat <- splines2::ibs(x = gridTime, knots = knots,
                                  degree = degree, intercept = TRUE,
                                  Boundary.knots = Boundary.knots)
            bMat <- attr(iMat, "bsMat")
        } else if (identical(spline, "mSplines")) {
            iMat <- splines2::iSpline(x = gridTime, knots = knots,
                                      degree = degree, intercept = TRUE,
                                      Boundary.knots = Boundary.knots)
            bMat <- attr(iMat, "msMat")
        } else
            stop("Unknown splines.")

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
        if (missing(na.action) || is.null(na.action))
            na.action <- options("na.action")[[1L]]

        if (missing(newdata)) {
            X <- matrix(rep(0, nBeta), nrow = 1)
            colnames(X) <- covnames
            rownames(X) <- "1"
        } else {
            mf <- stats::model.frame(Terms, newdata,
                                     na.action = na.action,
                                     xlev = object@xlevels)
            na.action <- paste0("na.", class(attr(mf, "na.action")))
            X <- stats::model.matrix(
                            Terms, mf,
                            contrasts.arg = object@contrasts$constracts
                        )
            ## remove intercept and deplicated rows
            X <- unique(base::subset(X, select = - `(Intercept)`))
            if (ncol(X) != nBeta)
                stop(wrapMessages(
                    "The number of input covariates does not",
                    "match with 'rateReg' object."
                ))
        }

        ## mcf for possible multigroups
        ndesign <- nrow(X)
        designSeq <- seq_len(ndesign)
        multiGroup <- ndesign > 1

        coveff <- as.numeric(exp(X %*% beta))
        outDat <- matrix(NA_real_, ncol = 4, nrow = ndesign * n_xx)

        for (i in designSeq) {
            ## Delta-method
            gradMat <- gradDelta(Xi = X[i, ], iMat, estMcf, coveffi = coveff[i])
            varTime <- apply(gradMat, 1, function (gradVec, covMat) {
                crossprod(gradVec, covMat) %*% gradVec
            }, covMat = covMat)
            confBand <- stats::qnorm((1 + level) / 2) * sqrt(varTime)
            mcf_i <- estMcf * coveff[i]
            lower <- pmax(0, mcf_i - confBand)
            upper <- mcf_i + confBand
            ind <- seq.int(from = n_xx * (i - 1) + 1, to = n_xx * i, by = 1)
            outDat[ind, 1L] <- gridTime
            outDat[ind, 2L] <- mcf_i
            outDat[ind, 3L] <- lower
            outDat[ind, 4L] <- upper
        }
        outDat <- as.data.frame(outDat)
        colnames(outDat) <- c("time", "MCF", "lower", "upper")
        if (multiGroup) {
            if (missing(groupLevels))
                groupLevels <- paste("Level", designSeq)
            if (missing(groupName))
                groupName <- "Group"
            tempcol <- factor(rep(groupLevels[designSeq], each = n_xx))
            outDat <- cbind(outDat, tempcol)
            colnames(outDat)[5L] <- groupName
        }
        ## return
        methods::new("mcf.rateReg",
                     call = object@call,
                     formula = object@formula,
                     spline = object@spline$spline,
                     knots = knots,
                     degree = degree,
                     Boundary.knots = Boundary.knots,
                     newdata = X,
                     MCF = outDat,
                     level = level,
                     na.action = na.action,
                     control = control,
                     multiGroup = multiGroup)
    })


### internal function ==========================================================
## Delta-method, compute the gradient on beta and alpha
gradDelta <- function(Xi, iMat, estMcf, coveffi) {
    gradBeta <- estMcf %o% Xi * coveffi
    gradAlpha <- iMat * coveffi
    ## retrun
    cbind(gradBeta, gradAlpha)
}


## control function
rateReg_mcf_control <- function(grid, length.out = 1e3, from, to, ...,
                                Boundary.knots_)
{
    ## controls for function MCF with signiture rateReg
    from <- if (missing(from)) Boundary.knots_[1L] else from
    to  <- if (missing(to)) Boundary.knots_[2L] else to
    if (! missing(grid)) {
        if (! isNumVector(grid))
            stop("The 'grid' specified must be numeric vector.",
                 call. = FALSE)
        grid <- sort(as.numeric(na.omit(grid)))
        length.out <- length(grid)
        from <- min(grid)
        to <- max(grid)
    } else {
        grid <- seq.int(from = from, to = to, length.out = length.out)
    }
    if (min(grid) < Boundary.knots_[1] || max(grid) > Boundary.knots_[2])
        stop(wrapMessages(
            "The 'grid' specified must be within",
            "the coverage of the boundary knots."
        ), call. = FALSE)

    ## return
    list(grid = grid, length.out = length.out, from = from, to = to)
}
