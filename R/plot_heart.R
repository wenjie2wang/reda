##############################################################################
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
##############################################################################


#' plot baseline.
#' 
#' \code{plot} produces the plot of baseline function.
#' 
#' This is a generic function using R base plotting system, 
#' which probably can be rewritten using ggplot2 later.
#' (This is a test Roxygen comments)
plot.baseline <- function(est.results, BaselinePieces, 
                          CI = TRUE, CI.level = 0.95) {
  n.xx <- 1000
  n.pieces <- length(BaselinePieces)
  BL.segments <- c(BaselinePieces[1], diff(BaselinePieces))
  xx <- seq(0, max(BaselinePieces), length = n.xx)
  indx <- sapply(xx, whereT, BaselinePieces = BaselinePieces)
  LinCom.M <- t(apply(array(indx), 1, function(ind.indx) {
    c(BL.segments[1:ind.indx], rep(0, n.pieces - ind.indx))}))
  CMF.B4.indx <- c(0, BaselinePieces)[indx]
  LinCom.M[(indx - 1) * n.xx + 1:n.xx] <- xx - CMF.B4.indx
  n.par <- dim(est.results$"Fisher Info Matrix")[1]
  Cov.M <- solve(est.results$"Fisher Info Matrix")[c((n.par - n.pieces + 1):n.par), 
                                                   c((n.par - n.pieces + 1):n.par)]
  CI.band <- qnorm((1 + CI.level)/2) * 
    sqrt(diag(LinCom.M %*% Cov.M %*% t(LinCom.M)))
  baseline.mean <- LinCom.M %*% est.results$alpha[, 1]
  ymax <- max(baseline.mean + CI.band)
  plot(xx, baseline.mean, type = "l", lwd = 2, ylim = c(0, ymax), 
       xlim = c(0, max(BaselinePieces)), xlab = "Time", 
       ylab = "MCF", main = "Mean Cumulative Function")
  lines(xx, baseline.mean + CI.band, lty = 3)
  lines(xx, baseline.mean - CI.band, lty = 3)
}


