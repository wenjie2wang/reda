################################################################################
##
##   R package heart by Haoda Fu, Jun Yan, and Wenjie Wang
##   Copyright (C) 2015
##
##   This file is part of the R package heart.
##
##   The R package heart is free software: you can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package heart is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package heart. If not, see <http://www.gnu.org/licenses/>.
##
################################################################################


#' Plot Mean Cumulative Function (MCF)
#' 
#' An S4 class generic function dispatched to a certain method 
#' to plot mean cumulative function by using ggplot2 plotting system. 
#' So the plots generated are able to be further customized properly.
#' 
#' @param object \code{\link{empirMCF-class}} or 
#' \code{\link{heartMCF-class}} object.
#' @param conf.int logical indicating whether to plot confidence interval.
#' The default value is FALSE.
#' @param ... other arguments for further usage.
#' @param linetypes an optional numeric vector indicating
#' line types specified to different groups with 
#' 0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 
#' 4 = dotdash, 5 = longdash, 6 = twodash.
#' @param linecolors an optional character vector indicating
#' line colors specified to different groups. 
#' @return ggplot object.
#' @seealso \code{\link{MCF}} 
#' @examples 
#' library(heart)
#' data(simuDat)
#' 
#' ## empirical MCF
#' sampleMCF <- MCF(Survr(ID, time, event) ~ group, data = simuDat)
#' plotMCF(sampleMCF, linetypes = c(1, 3), linecolors = c("orange", "navy"))
#' 
#' ## estimated MCF for baseline rate function from HEART model
#' heartfit <- heart(formula = Survr(ID, time, event) ~ X1 + group, 
#'                   data = simuDat, subset = ID %in% 75:125,
#'                   baselinepieces = seq(28, 168, length = 6))
#' baselineMCF <- MCF(heartfit)
#' plotMCF(baselineMCF, conf.int = TRUE, linecolors = "blueviolet") + 
#'   ggplot2::theme_bw()
#' @export
setGeneric(name = "plotMCF",
           def = function(object, conf.int = FALSE, ...) {
             standardGeneric("plotMCF")
           })


#' @describeIn plotMCF Plot empirical mean cumulative function (MCF)
#' @aliases plotMCF,empirMCF-method
#' @importFrom utils tail
#' @importFrom stats setNames
#' @importFrom ggplot2 ggplot geom_step aes aes_string scale_color_manual
#' scale_linetype_manual ylab ggtitle
#' @export
setMethod(f = "plotMCF", signature = "empirMCF", 
          definition = function(object, conf.int = FALSE, 
                                linetypes, linecolors, ...) {
            MCFdat <- object@MCF
            ## if MCF is just for one certain group
            if (! object@multigroup) {
              p <- ggplot2::ggplot(data = MCFdat, 
                                   ggplot2::aes_string(x = "Time")) + 
                ggplot2::geom_step(mapping = ggplot2::aes(x = time, y = MCF))
            } else {
              ## function to emulate the default colors used in ggplot2
              gg_color_hue <- function(n){
                hues = seq(15, 375, length=n+1)
                return(hcl(h=hues, l=65, c=100)[1:n])
              }
              
              legendname <- utils::tail(colnames(MCFdat), n = 1)
              MCFdat$design <- MCFdat[, legendname]
              Design <- factor(MCFdat$design)
              ndesign = length(levels(Design))
              
              ## about linetypes
              # 0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 
              # 4 = dotdash, 5 = longdash, 6 = twodash
              ## set line types and colors
              if(missing(linetypes)){
                lts <- stats::setNames(rep(1, ndesign), levels(Design)) 
              }else{
                lts <- stats::setNames(linetypes[seq(ndesign)], levels(Design))
              }
              if(missing(linecolors)){
                lcs <- stats::setNames(gg_color_hue(ndesign), levels(Design))
              }else{
                lcs <- stats::setNames(linecolors[seq(ndesign)], levels(Design))
              }
              p <- ggplot2::ggplot(data = MCFdat, 
                                   ggplot2::aes_string(x = "Time")) +
                ggplot2::geom_step(mapping = ggplot2::aes(x = time, y = MCF, 
                                                          color = design, 
                                                          linetype = design)) +
                ggplot2::scale_color_manual(values = lcs, name = legendname) +
                ggplot2::scale_linetype_manual(values= lts, name = legendname)
            }
            p <- p + ggplot2::ylab("MCF") + 
              ggplot2::ggtitle("Empirical Mean Cumulative Function")
            return(p)
          })


#' @describeIn plotMCF Estimated mean cumulative function (MCF) 
#' for baseline rate function.
#' @aliases plotMCF,heartMCF-method
#' @importFrom stats setNames
#' @importFrom ggplot2 ggplot geom_line aes aes_string scale_color_manual
#' scale_linetype_manual ylab ggtitle 
#' @export
setMethod(f = "plotMCF", signature = "heartMCF", 
          definition = function(object, conf.int = FALSE, 
                                linetypes, linecolors, ...) {
            MCFdat <- object@MCF
            ## if MCF is just for one certain group
            if (! object@multigroup) {
              p <- ggplot2::ggplot(data = MCFdat, 
                                   ggplot2::aes_string(x = "Time")) + 
                ggplot2::geom_line(mapping = ggplot2::aes(x = time, y = MCF))
              if (conf.int) {
                p <- p + 
                  ggplot2::geom_line(
                    mapping = ggplot2::aes(x = time, y = lower), 
                    linetype = "3313") +
                  ggplot2::geom_line(
                    mapping = ggplot2::aes(x = time, y = upper), 
                    linetype = "3313")
              }
            } else {
              ## function to emulate the default colors used in ggplot2
              gg_color_hue <- function(n){
                hues = seq(15, 375, length=n+1)
                return(hcl(h=hues, l=65, c=100)[1:n])
              }
              
              legendname <- tail(colnames(MCFdat), n = 1)
              MCFdat$Design <- MCFdat[, legendname]
              Design <- factor(MCFdat$Design)
              ndesign = length(levels(Design))
              
              ## about linetypes
              # 0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 
              # 4 = dotdash, 5 = longdash, 6 = twodash
              ## set line types and colors
              if(missing(linetypes)){
                lts <- stats::setNames(rep(1, ndesign), levels(Design)) 
              }else{
                lts <- stats::setNames(linetypes[seq(ndesign)], levels(Design))
              }
              if(missing(linecolors)){
                lcs <- stats::setNames(gg_color_hue(ndesign), levels(Design))
              }else{
                lcs <- stats::setNames(linecolors[seq(ndesign)], levels(Design))
              }
              p <- ggplot2::ggplot(data = MCFdat, 
                                   ggplot2::aes_string(x = "Time")) +
                ggplot2::geom_line(mapping = ggplot2::aes(x = time, y = MCF, 
                                                          color = Design, 
                                                          linetype = Design)) +
                ggplot2::scale_color_manual(values = lcs, name = legendname) +
                ggplot2::scale_linetype_manual(values= lts, name = legendname)
              if (conf.int) {
                p <- p + 
                  ggplot2::geom_line(
                    mapping = ggplot2::aes(x = time, y = lower, color = Design), 
                    linetype = "3313") +
                  ggplot2::geom_line(
                    mapping = ggplot2::aes(x = time, y = upper, color = Design), 
                    linetype = "3313")
              }
            }
            p <- p + ggplot2::ylab("MCF") + 
              ggplot2::ggtitle("Estimated Mean Cumulative Function")
            return(p)
          })


