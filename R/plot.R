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


##' Plot Mean Cumulative Function (MCF)
##'
##'
##' S4 class methods plotting sample MCF from data or
##' estimated MCF from a fitted model
##' by using \code{ggplot2} plotting system.
##' The plots generated are able to be further customized properly.
##'
##' @name plotMcf-method
##' @param x An object used to dispatch a method.
##' @param y An argument that should be missing and ignored now.
##' Its existence is just for satisfying the definition of generaic function
##' \code{plot} in package \code{graphics} for methods' dispatching.
##' @param conf.int A logical value indicating
##' whether to plot confidence interval.
##' The default value is \code{FALSE}.
##' @param ... Other arguments for further usage.
##' @param mark.time A logical value with default \code{FALSE}.
##' If \code{TRUE}, each censoring time is marked by "+" on the MCF curves.
##' Otherwise, the censoring time would not be marked.
##' @param lty An optional numeric vector indicating
##' line types specified to different groups:
##' 0 = blank, 1 = solid, 2 = dashed, 3 = dotted,
##' 4 = dotdash, 5 = longdash, 6 = twodash.
##' @param col An optional character vector indicating
##' line colors specified to different groups.
##' @return A \code{ggplot} object.
##' @examples
##' ## See examples given in function mcf and rateReg.
##' @seealso
##' \code{\link{mcf}} for estimation of MCF;
##' \code{\link{rateReg}} for model fitting.
##' @importFrom graphics plot
NULL

##' @rdname plotMcf-method
##' @aliases plot,sampleMcf-method
##' @param legendName An optional length-one charactor vector to specify the
##' name for grouping each unique row in \code{newdata}, such as "gender"
##' for "male" and "female". The default value is generated from the
##' \code{object}.
##' @param legendLevels An optional charactor vector to specify the levels for
##' each unique row in \code{newdata}, such as "treatment" and "control".
##' The default values are generated from the \code{object}.
##' @importFrom ggplot2 ggplot geom_step aes aes_string scale_color_manual
##' scale_linetype_manual ylab ggtitle geom_text theme element_text
##' @export
setMethod(
    f = "plot", signature = c("sampleMcf", "missing"),
    definition = function(x, y, conf.int = FALSE, mark.time = FALSE,
                          lty, col, legendName, legendLevels, ...) {

        ## nonsense, just to suppress Note from R CMD check --as-cran
        MCF <- event <- lower <- upper <- design <- time <- NULL

        ## rename first three columns
        MCFdat <- x@MCF
        colnames(MCFdat)[seq_len(3)] <- c("ID", "time", "event")

        ## if MCF is just for one certain group
        if (! x@multiGroup) {
            ## add starting point at time 0
            MCFdat <- MCFdat[c(1L, seq_len(nrow(MCFdat))), ]
            MCFdat[1L, 2 : 7] <- 0
            if (missing(lty)) lty <- 1
            if (missing(col)) col <- "black"
            p <- ggplot(data = MCFdat, aes_string(x = "Time")) +
                geom_step(mapping = aes(x = time, y = MCF),
                          linetype = lty, color = col)
            ## mark censoring time
            if (mark.time) {
                cenDat <- base::subset(MCFdat, time > 0 & event == 0)
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
            desDat <- MCFdat[, - seq_len(7), drop = FALSE]
            groupName <- paste(colnames(desDat), collapse = "&")
            desList <- as.list(desDat)
            MCFdat$design <- factor(do.call(paste, c(desList, sep = "&")))
            nDesign = length(desLevs <- levels(MCFdat$design))

            ## set possibly customized group name and levels
            legendName <- if (missing(legendName)) {
                              groupName
                          } else {
                              as.character(legendName)[1L]
                          }
            if (! missing(legendLevels)) {
                if (length(legendLevels) != nDesign)
                    stop(paste("The length of 'legendLevels' must",
                               "match the number of designs."))
                desLevs <- levels(MCFdat$design) <-
                    as.character(legendLevels)
            }

            ## add starting point at time 0 before each level
            idx <- which(! duplicated(MCFdat$design))
            sortIdx <- sort(c(idx, seq_len(nrow(MCFdat))))
            MCFdat <- MCFdat[sortIdx, ]
            desInd <- seq_len(nDesign)
            MCFdat[idx + desInd - 1, 2 : 7] <- 0

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
                cenDat <- base::subset(MCFdat, time > 0 & event == 0)
                p <- p + geom_text(data = cenDat,
                                   aes(label = "+", x = time, y = MCF,
                                       linetype = design, color = design),
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
        p <- p + ylab("MCF") + ggtitle("Sample Mean Cumulative Function")
        ## centered alignment
        p <- p + theme(plot.title = element_text(hjust = 0.5))
        p
    })


##' @rdname plotMcf-method
##' @aliases plot,rateRegMcf-method
##' @importFrom ggplot2 ggplot geom_line aes aes_string scale_color_manual
##' scale_linetype_manual ylab ggtitle theme element_text
##' @export
setMethod(
    f = "plot", signature = c("rateRegMcf", "missing"),
    definition = function(x, y, conf.int = FALSE, lty, col, ...) {

        ## nonsense, just to suppress Note from R CMD check --as-cran
        MCF <- lower <- upper <- time <- NULL

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
        p <- p + ylab("MCF") + ggtitle("Estimated Mean Cumulative Function")
        ## centered alignment
        p <- p + theme(plot.title = element_text(hjust = 0.5))
        p
    })



### internal function ==========================================================
## function to emulate the default colors used in ggplot2
##' @importFrom grDevices hcl
gg_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100)[seq_len(n)]
}
