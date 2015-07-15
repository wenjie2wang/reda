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

#' Plot Empirical Mean Cumulative Function (MCF) from Sample Data.
#' 
#' \code{plot_sampleMCF} plots empirical mean cumulative function (MCF).  
#' 
#' The function plots empirical mean cumulative function 
#' by using ggplot2 plotting system.  
#' So the plots generated are able to be further customized properly.
#' @usage plot_sampleMCF(MCF, linetypes, linecolors, ...)
#' @param MCF data.frame generated from function 
#' \code{\link[heart]{sample_MCF}}.
#' @param linetypes line types specified to different groups.
#' @param linecolors line colors specified to different groups.
#' @param ... further arguments.
#' @return ggplot object.
#' @export
plot_sampleMCF <- function(MCF, linetypes, linecolors, ...) {
  ## if it is overall MCF
  if (ncol(MCF) == 5) {
    p <- ggplot2::ggplot(data = MCF, ggplot2::aes_string(x = "Time")) + 
      ggplot2::geom_step(mapping = ggplot2::aes(x = time, y = MCF))
   } else {

  ## function to emulate the default colors used in ggplot2
  gg_color_hue <- function(n){
    hues = seq(15, 375, length=n+1)
    return(hcl(h=hues, l=65, c=100)[1:n])
  }
  
  legendname <- tail(colnames(MCF), n = 1)
  MCF$group <- MCF[, legendname]
  Design <- factor(MCF$group)
  ndesign = length(levels(Design))
  
  ## about linetypes
  # 0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 
  # 4 = dotdash, 5 = longdash, 6 = twodash
  ## set line types and colors
  if(missing(linetypes)){
    lts <- setNames(rep(1, ndesign), levels(Design)) 
  }else{
    lts <- setNames(linetypes[seq(ndesign)], levels(Design))
  }
  if(missing(linecolors)){
    lcs <- gg_color_hue(ndesign)
  }else{
    lcs <- linecolors[seq(ndesign)]
  }
  p <- ggplot2::ggplot(data = MCF, ggplot2::aes_string(x = "Time")) +
    ggplot2::geom_step(mapping = ggplot2::aes(x = time, y = MCF, 
                                              color = group, 
                                              linetype = group)) +
    ggplot2::scale_color_manual(values = lcs, name = legendname) +
    ggplot2::scale_linetype_manual(values= lts, name = legendname)
  }
  p <- p + ggplot2::ylab("MCF") + 
    ggplot2::ggtitle("Empirical Mean Cumulative Function") +
    ggplot2::xlim(0, max(MCF$time)) + 
    ggplot2::ylim(0, max(MCF$MCF))
  return(p)
}

#' Plot Estimated Mean Cumulative Function (MCF) of Baseline Rate Function.
#' 
#' \code{plot_heartMCF} produces the plot of baseline function.
#' 
#' The function estimates mean cumulative function of baseline rate function
#' by using ggplot2 plotting system.  
#' So the plots generated are able to be further customized properly.
#' @usage plot_MCF(object, CI = TRUE, level = 0.95, ...)
#' @param object heart object.
#' @param CI logic value, TRUE or FALSE to specify whether to include confidence
#' interval in the plot. 
#' @param level numeric value \eqn{\in (0, 1)} confidence level.
#' @param ... further arguments.
#' @return ggplot object
#' @export
plot_heartMCF <- function(object, CI = TRUE, level = 0.95, ...) {
  baselinepieces <- object@baselinepieces
  n_xx <- 1000
  n_pieces <- length(baselinepieces)
  BL_segments <- c(baselinepieces[1], diff(baselinepieces))
  xx <- seq(0, max(baselinepieces), length = n_xx)
  indx <- sapply(xx, whereT, baselinepieces = baselinepieces)
  LinCom_M <- t(apply(array(indx), 1, 
                      function(ind_indx) {
                        c(BL_segments[1:ind_indx], rep(0, n_pieces - ind_indx))
                      }))
  CMF_B4_indx <- c(0, baselinepieces)[indx]
  LinCom_M[(indx - 1) * n_xx + 1:n_xx] <- xx - CMF_B4_indx
  n_par <- dim(object@hessian)[1]
  Cov_M <- solve(object@hessian)[c((n_par - n_pieces + 1):n_par), 
                                 c((n_par - n_pieces + 1):n_par)]
  CI_band <- qnorm((1 + level)/2) * 
    sqrt(diag(LinCom_M %*% Cov_M %*% t(LinCom_M)))
  baseline_mean <- LinCom_M %*% object@estimates$alpha[, 1]
  ymax <- max(baseline_mean + CI_band)
  lower <- baseline_mean - CI_band
  upper <- baseline_mean + CI_band
  plotdat <- data.frame(time = xx, MCF = baseline_mean, 
                        lower = lower, upper = upper)
  p <- ggplot2::ggplot(data = plotdat, ggplot2::aes_string(x = "Time"))
  p <-  p + ggplot2::geom_line(mapping = ggplot2::aes(x = time, y = MCF))
  if (CI) {
    p <- p + 
      ggplot2::geom_line(mapping = ggplot2::aes(x = time, y = lower), 
                         linetype = "3313") +
      ggplot2::geom_line(mapping = ggplot2::aes(x = time, y = upper), 
                         linetype = "3313")
  }
  p <- p + ggplot2::ylab("MCF") + ggplot2::ggtitle("Mean Cumulative Function") +
    ggplot2::xlim(0, max(baselinepieces)) + 
    ggplot2::ylim(0, ymax)
  ## return
  p
}

