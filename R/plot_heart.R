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


#' plot baseline.
#' 
#' \code{plot} produces the plot of baseline function.
#' 
#' This is a generic function using R base plotting system, 
#' which probably can be rewritten using ggplot2 later.
#' (This is a test Roxygen comments)
plot.baseline <- function(object, CI = TRUE, level = 0.95) {
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
  plot(xx, baseline_mean, type = "l", lwd = 2, ylim = c(0, ymax), 
       xlim = c(0, max(baselinepieces)), xlab = "Time", 
       ylab = "MCF", main = "Mean Cumulative Function")
  lines(xx, baseline_mean + CI_band, lty = 3)
  lines(xx, baseline_mean - CI_band, lty = 3)
}


