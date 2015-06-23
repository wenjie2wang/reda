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


whereT <- function(tt, BaselinePieces) {
  min(which(tt <= BaselinePieces))
}


rho_0 <- function(par.BaselinePW, BaselinePieces, Tvec) {
  indx <- apply(as.array(Tvec), 1, whereT, BaselinePieces)
  return(par.BaselinePW[indx])
}


mu0 <- function(par.BaselinePW, BaselinePieces, Tvec) {
  indx <- apply(as.array(Tvec), 1, whereT, BaselinePieces)
  BL.segments <- c(BaselinePieces[1], diff(BaselinePieces))
  #The CMF at each time point  
  CumMean.Pieces <- diffinv(BL.segments * par.BaselinePW)[-1]  
  CumMean.Pieces[indx] - (BaselinePieces[indx] - Tvec) * par.BaselinePW[indx]
}


dmu0_alpha <- function(tt, BaselinePieces) {
  BL.segments <- c(BaselinePieces[1], diff(BaselinePieces))
  indx <- min(which(tt <= BaselinePieces))
  value <- BL.segments
  n.pieces <- length(BaselinePieces)
  if (indx == n.pieces) {
    value[n.pieces] <- tt - BaselinePieces[n.pieces - 1]
  } else if (indx > 1) {
    value[(indx + 1):n.pieces] <- 0
    value[indx] <- tt - BaselinePieces[indx - 1]
  } else {
    value[(indx + 1):n.pieces] <- 0
    value[indx] <- tt
  }
  return(value)
}


logL.HEART <- function(par, dataM, BaselinePieces) {
  npieces <- length(BaselinePieces)
  if (BaselinePieces[npieces] < max(dataM$Time)) {
    BaselinePieces[npieces] <- max(dataM$Time) + 1e-08
    warning("Extend the Baseline Pieces to adjust the data")
  }
  nbeta <- dim(dataM)[2] - 3
  par.beta <- par[1:nbeta]
  par.theta <- par[nbeta + 1]
  par.alpha <- par[(nbeta + 2):length(par)]
  m <- length(unique(dataM$ID))
  expXBeta <- exp(as.matrix(dataM[, 4:(3 + nbeta)]) %*% (as.matrix(par.beta)))
  rho_0_ij <- rho_0(par.BaselinePW = par.alpha, BaselinePieces = BaselinePieces, 
                    dataM$Time[dataM$Event == 1])
  rho_i <- expXBeta[dataM$Event == 1] * rho_0_ij
  rho_i[rho_i < 1e-100] <- 1e-100
  sum.log.rho_i <- sum(log(rho_i))
  # these codes to make sure that the order will not change 
  # if the patient ID is not ordered
  n_ij <- table(dataM$ID)[order(unique(dataM$ID))] - 1  
  # if there is a subject with 0 event, 
  # the sequence will not be generated for this subject
  theta_j_1 <- par.theta + sequence(n_ij) - 1  
  theta_j_1[theta_j_1 < 1e-100] <- 1e-100
  sum.log_theta_j_1 <- sum(log(theta_j_1))
  mu0i <- mu0(par.BaselinePW = par.alpha, BaselinePieces = BaselinePieces, 
              dataM$Time[dataM$Event == 0])
  mui <- mu0i * expXBeta[dataM$Event == 0]
  mui_theta <- par.theta + mui
  mui_theta[mui_theta < 1e-100] <- 1e-100
  sum.log_theta_mui <- sum((n_ij + par.theta) * log(mui_theta))
  if (par.theta < 1e-100) {
    par.theta <- 1e-100
  }
  logLH <- m * par.theta * log(par.theta) + sum.log.rho_i + 
    sum.log_theta_j_1 - sum.log_theta_mui
  penal <- ifelse(par.theta < 0 | min(par.alpha) < 0, 1e+50, 0)
  negLH <- -logLH + penal
  ### Calculate the gradient
  dl_dbeta <- apply(diag((n_ij - mui)/(par.theta + mui) * par.theta) %*% 
                      as.matrix(dataM[dataM$Event == 0, 4:(3 + nbeta)]), 2, sum)
  dl_dtheta <- m + m * log(par.theta) + 
    sum(1/(par.theta + sequence(n_ij) - 1)) - 
    sum((n_ij + par.theta)/(par.theta + mui)) - sum(log(mui_theta))
  indx <- apply(as.array(dataM$Time[dataM$Event == 1]), 1, 
                whereT, BaselinePieces)
  if (length(unique(indx)) < length(par.alpha)) {
    stop("Some segements have zero events!")
  }
  indx_taui <- apply(as.array(dataM$Time[dataM$Event == 0]), 1, 
                     whereT, BaselinePieces)
  dl_dalpha.part2 <- diag((n_ij + par.theta) / (par.theta + mui) * 
                            expXBeta[dataM$Event == 0]) %*% 
    t(apply(array(dataM$Time[dataM$Event == 0]), 1, dmu0_alpha, BaselinePieces))
  dl_dalpha <- 1 / par.alpha * table(indx) - apply(dl_dalpha.part2, 2, sum)
  attr(negLH, "gradient") <- -c(dl_dbeta, dl_dtheta, dl_dalpha)
  return(negLH)
}




#' Fitting Heart Model: Piece-wise Gamma Frailty Model for Recurrent Events. 
#' The model is named after the paper title of \emph{Fu et al. (2014)},  
#' Hypoglycemic Events Analysis via Recurrent Time-to-Event (HEART) Models
#'
#' \code{functionname} returns fitted model results.
#'
#' This is a test Roxygen comments
#'
#' @param ... Numeric, complex, or logical vectors.
#' @param na.rm A logical scalar. 
#' @return If all inputs are integer and logical, then the output
#'   will be an integer. If integer overflow
#'   \url{http://en.wikipedia.org/wiki/Integer_overflow} occurs, the output
#'   will be NA with a warning. Otherwise it will be a length-one numeric or
#'   complex vector.
#' @examples
#' sum(1:10)
#' sum(1:5, 6:10)
#' sum(F, F, F, T, T)
heart <- function(dataM, BaselinePieces, ini) {
  fit <- nlm(logL.HEART, ini, dataM = dataM, 
             BaselinePieces = BaselinePieces, hessian = T)
  length.par <- length(ini)
  length.alpha <- length(BaselinePieces)
  length.beta <- length.par - length.alpha - 1
  est.beta <- matrix(NA, nrow = length.beta, ncol = 3)
  colnames(est.beta) <- c("beta", "se(beta)", "two sided p-value")
  
  se.vec <- sqrt(diag(solve(fit$hessian)))
  
  est.beta[, 1] <- fit$estimate[1:length.beta]
  est.beta[, 2] <- se.vec[1:length.beta]
  est.beta[, 3] <- 2 * pnorm(-abs(est.beta[, 1]/est.beta[, 2]))
  
  est.theta <- matrix(NA, nrow = 1, ncol = 2)
  colnames(est.theta) <- c("theta", "se(theta)")
  est.theta[1, ] <- c(fit$estimate[length.beta + 1], se.vec[length.beta + 1])
  
  est.alpha <- matrix(NA, nrow = length.alpha, ncol = 2)
  colnames(est.alpha) <- c("alpha", "se(alpha)")
  
  est.alpha[, 1] <- fit$estimate[(length.beta + 2):length.par]
  est.alpha[, 2] <- se.vec[(length.beta + 2):length.par]
  results <- list(beta = est.beta, theta = est.theta, alpha = est.alpha, 
                  Convergence = fit$code, `Fisher Info Matrix` = fit$hessian)
  return(results)
}

