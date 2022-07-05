##
## R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
## Copyright (C) 2015-2022
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


##' Plot Baseline Rate or Mean Cumulative Function (MCF)
##'
##' S4 class methods plotting sample MCF from data, estimated MCF, or estimated
##' baseline hazard rate function from a fitted model by using \code{ggplot2}
##' plotting system.  The plots generated are thus able to be further customized
##' properly.
##'
##' @name plot-method
##'
##' @param x An object used to dispatch a method.
##' @param y An argument that should be missing and ignored now.  Its existence
##'     is just for satisfying the definition of generaic function \code{plot}
##'     in package \code{graphics} for methods' dispatching.
##' @param conf.int A logical value indicating whether to plot confidence
##'     interval.  The default value is \code{FALSE}.
##' @param lty An optional numeric vector indicating line types specified to
##'     different groups: 0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 4 =
##'     dotdash, 5 = longdash, 6 = twodash.
##' @param col An optional character vector indicating line colors specified to
##'     different groups.
##' @param ... Other arguments for further usage.
##'
##' @return A \code{ggplot} object.
##'
##' @examples
##' ## See examples given in function mcf and rateReg.
##'
##' @seealso
##' \code{\link{mcf}} for estimation of MCF;
##' \code{\link{rateReg}} for model fitting.
##'
##' @importFrom ggplot2 aes aes_string element_text geom_line geom_step
##'     geom_text ggplot ggtitle scale_color_manual scale_linetype_manual theme
##'     ylab
##'
##' @importFrom graphics plot
NULL


##' @rdname plot-method
##'
##' @aliases plot,mcf.formula-method
##'
##' @param legendName An optional length-one charactor vector to specify the
##'     name for grouping each unique row in \code{newdata}, such as "gender"
##'     for "male" and "female". The default value is generated from the
##'     \code{object}.
##' @param legendLevels An optional charactor vector to specify the levels for
##'     each unique row in \code{newdata}, such as "treatment" and "control".
##'     The default values are generated from the \code{object}.
##' @param mark.time A logical value with default value \code{FALSE}.  If
##'     \code{TRUE}, each censoring time is marked by "+" on the MCF curves.
##'     Otherwise, the censoring time would not be marked.
##' @param addOrigin A logical value indicating whether the MCF curves start
##'     from origin time. The default value is \code{FALSE}.
##'
##' @export
setMethod(
    f = "plot", signature = c("mcf.formula", "missing"),
    definition = function(x, y,
                          lty, col,
                          legendName, legendLevels,
                          conf.int = FALSE,
                          mark.time = FALSE,
                          addOrigin = FALSE,
                          ...)
    {
        ## nonsense, just to suppress Note from R CMD check --as-cran
        MCF <- instRate <- lower <- upper <- design <- time <- NULL

        ## mcf data
        MCFdat <- x@MCF
        ## if it is overall sample MCF
        if (! x@multiGroup) {
            if (addOrigin) {
                ## add starting point at origin time
                originDat <- MCFdat[1L, ]
                originDat[, ] <- 0
                originDat[, "time"] <- x@origin
                originDat[, "instRate"] <- NA
                MCFdat <- rbind(originDat, MCFdat)
            }
            ## set default line type and color
            if (missing(lty)) lty <- 1
            if (missing(col)) col <- "black"
            p <- ggplot(data = MCFdat, aes_string(x = "Time")) +
                geom_step(mapping = aes(x = time, y = MCF),
                          linetype = lty, color = col)
            ## mark censoring time
            if (mark.time) {
                cenDat <- base::subset(MCFdat, ! is.na(instRate) &
                                               instRate <= 0)
                p <- p + geom_text(data = cenDat,
                                   aes(label = "+", x = time, y = MCF),
                                   vjust = 0.3, hjust = 0.5,
                                   color = col, show.legend = FALSE)
            }
            ## confidence interval
            if (conf.int) {
                p <- p + geom_step(mapping = aes(x = time, y = lower),
                                   linetype = "3313", color = col) +
                    geom_step(mapping = aes(x = time, y = upper),
                              linetype = "3313", color = col)
            }
        } else {
            ## if it is MCF for multiple groups
            ## caution: this value may vary if mcf-formula method changes
            numColMcf <- 7L
            desDat <- MCFdat[, - seq_len(numColMcf), drop = FALSE]
            groupName <- paste(colnames(desDat), collapse = " & ")
            ## keep order for factor variables
            desDat <- desDat[do.call(order, as.list(desDat)), , drop = FALSE]
            desVec <- do.call(paste, c(as.list(desDat), sep = " & "))
            desLevs <- unique(desVec)
            MCFdat$design <- factor(desVec, levels = desLevs)
            nDesign <- length(desLevs)
            desInd <- seq_len(nDesign)
            ## set possibly customized group name and levels
            legendName <- if (missing(legendName)) {
                              groupName
                          } else {
                              as.character(legendName)[1L]
                          }
            if (! missing(legendLevels)) {
                if (length(legendLevels) != nDesign)
                    stop(wrapMessages(
                        "The length of 'legendLevels' must ",
                        "match the number of designs."
                    ), call. = FALSE)
                desLevs <- as.character(legendLevels)
                MCFdat$design <- factor(desVec,
                                        levels = levels(MCFdat$design),
                                        labels = desLevs)
            }
            if (addOrigin) {
                ## add starting point at origin time for each group
                originDat <- MCFdat[rep(1L, nDesign),]
                originDat[, 2 : numColMcf] <- 0
                originDat[, "time"] <- x@origin
                originDat[, "design"] <- desLevs
                originDat[, "instRate"] <- NA
                MCFdat <- rbind(originDat, MCFdat)
            }
            ## about lty
            ## 0 = blank, 1 = solid, 2 = dashed, 3 = dotted,
            ## 4 = dotdash, 5 = longdash, 6 = twodash
            ## set line types and colors
            lts <- if (missing(lty)) {
                       stats::setNames(rep(1, nDesign), desLevs)
                   } else {
                       stats::setNames(lty[desInd], desLevs)
                   }
            lcs <- if (missing(col)) {
                       stats::setNames(gg_color_hue(nDesign), desLevs)
                   } else {
                       stats::setNames(col[desInd], desLevs)
                   }
            p <- ggplot(data = MCFdat, aes_string(x = "Time")) +
                geom_step(mapping = aes(x = time, y = MCF, color = design,
                                        linetype = design)) +
                scale_color_manual(values = lcs, name = legendName) +
                scale_linetype_manual(values= lts, name = legendName)
            if (mark.time) {
                cenDat <- base::subset(MCFdat, ! is.na(instRate) &
                                               instRate <= 0)
                p <- p + geom_text(data = cenDat,
                                   aes(label = "+", x = time,
                                       y = MCF, color = design),
                                   vjust = 0.3, hjust = 0.5,
                                   show.legend = FALSE)
            }
            if (conf.int) {
                p <- p + geom_step(mapping = aes(x = time, y = lower,
                                                 color = design),
                                   linetype = "3313") +
                    geom_step(mapping = aes(x = time, y = upper,
                                            color = design),
                              linetype = "3313")
            }
        }
        p + ylab("MCF Estimates")
    })


##' @rdname plot-method
##'
##' @aliases plot,mcf.rateReg-method
##'
##' @export
setMethod(
    f = "plot", signature = c("mcf.rateReg", "missing"),
    definition = function(x, y, conf.int = FALSE, lty, col, ...)
    {
        ## nonsense, just to suppress Note from R CMD check --as-cran
        MCF <- lower <- upper <- time <- NULL
        ## mcf data
        MCFdat <- x@MCF
        ## if MCF is just for one certain group
        if (! x@multiGroup) {
            if (missing(lty)) lty <- 1
            if (missing(col)) col <- "black"
            p <- ggplot(data = MCFdat, aes_string(x = "Time")) +
                geom_line(mapping = aes(x = time, y = MCF),
                          linetype = lty, color = col)
            if (conf.int) {
                p <- p + geom_line(mapping = aes(x = time, y = lower),
                                   linetype = "3313", color = col) +
                    geom_line(mapping = aes(x = time, y = upper),
                              linetype = "3313", color = col)
            }
        } else {
            legendName <- colnames(MCFdat)[ncol(MCFdat)]
            MCFdat$Design <- MCFdat[, legendName]
            Design <- factor(MCFdat$Design)
            nDesign = length(desLevs <- levels(Design))
            desInd <- seq_len(nDesign)

            ## about lty
            ## 0 = blank, 1 = solid, 2 = dashed, 3 = dotted,
            ## 4 = dotdash, 5 = longdash, 6 = twodash
            ## set line types and colors
            lts <- if (missing(lty)) {
                       stats::setNames(rep(1, nDesign), desLevs)
                   } else {
                       stats::setNames(lty[desInd], desLevs)
                   }
            lcs <- if (missing(col)) {
                       stats::setNames(gg_color_hue(nDesign), desLevs)
                   } else {
                       stats::setNames(col[desInd], desLevs)
                   }
            p <- ggplot(data = MCFdat, aes_string(x = "Time")) +
                geom_line(mapping = aes(x = time, y = MCF, color = Design,
                                        linetype = Design)) +
                scale_color_manual(values = lcs, name = legendName) +
                scale_linetype_manual(values= lts, name = legendName)
            if (conf.int) {
                p <- p + geom_line(mapping = aes(x = time, y = lower,
                                                 color = Design),
                                   linetype = "3313") +
                    geom_line(mapping = aes(x = time, y = upper,
                                            color = Design),
                              linetype = "3313")
            }
        }
        p + ylab("MCF Estimates")
    })


##' @rdname plot-method
##'
##' @aliases plot,baseRate.rateReg-method
##'
##' @export
setMethod(
    f = "plot", signature = c("baseRate.rateReg", "missing"),
    definition = function(x, y, conf.int = FALSE, lty, col, ...)
    {
        ## nonsense, just to suppress Note from R CMD check --as-cran
        lower <- upper <- time <- NULL

        dat <- x@baseRate
        if (missing(lty)) lty <- 1
        if (missing(col)) col <- "black"
        p <- ggplot(data = dat, aes_string(x = "Time")) +
            geom_line(mapping = aes(x = time, y = baseRate),
                      linetype = lty, color = col)
        if (conf.int) {
            p <- p + geom_line(mapping = aes(x = time, y = lower),
                               linetype = "3313", color = col) +
                geom_line(mapping = aes(x = time, y = upper),
                          linetype = "3313", color = col)
        }
        p + ylab("Baseline Hazard Rate Estimates")
    })


##' @rdname plot-method
##' @aliases plot,mcfDiff-method
##' @export
setMethod(
    f = "plot", signature = c("mcfDiff", "missing"),
    definition = function(x, y,
                          lty, col,
                          legendName, legendLevels,
                          conf.int = TRUE,
                          addOrigin = FALSE,
                          ...)
    {
        ## nonsense, just to suppress Note from R CMD check --as-cran
        MCF <- lower <- upper <- time <- NULL

        ## mcf data
        MCFdat <- x@MCF
        if (addOrigin) {
            ## add starting point at origin time
            originDat <- MCFdat[1L, ]
            originDat[1L, ] <- 0
            originDat[, "time"] <- min(x@origin)
            MCFdat <- rbind(originDat, MCFdat)
        }
        ## set default line type and color
        if (missing(lty)) lty <- 1
        if (missing(col)) col <- "black"
        p <- ggplot(data = MCFdat, aes_string(x = "Time")) +
            geom_step(mapping = aes(x = time, y = MCF),
                      linetype = lty, color = col)
        ## confidence interval
        if (conf.int) {
            p <- p + geom_step(mapping = aes(x = time, y = lower),
                               linetype = "3313", color = col) +
                geom_step(mapping = aes(x = time, y = upper),
                          linetype = "3313", color = col)
        }
        p <- p + ylab("MCF Difference")
        p
    })


### internal function ==========================================================
## function to emulate the default colors used in ggplot2
##' @importFrom grDevices hcl
gg_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100)[seq_len(n)]
}
