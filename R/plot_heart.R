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


