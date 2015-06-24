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
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOse_  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package heart. If not, see <http://www.gnu.org/licenses/>.
##
##############################################################################


## internal functions ========================================================
whereT <- function(tt, BaselinePieces) {
  ## return
  min(which(tt <= BaselinePieces))
}

rho_0 <- function(par_BaselinePW, BaselinePieces, Tvec) {
  indx <- apply(as.array(Tvec), 1, whereT, BaselinePieces)
  ## return
  par_BaselinePW[indx]
}

mu0 <- function(par_BaselinePW, BaselinePieces, Tvec) {
  indx <- apply(as.array(Tvec), 1, whereT, BaselinePieces)
  BL_segments <- c(BaselinePieces[1], diff(BaselinePieces))
  ## The CMF at each time point  
  CumMean_Pieces <- diffinv(BL_segments * par_BaselinePW)[-1]  
  ## return
  CumMean_Pieces[indx] - (BaselinePieces[indx] - Tvec) * par_BaselinePW[indx]
}

dmu0_alpha <- function(tt, BaselinePieces) {
  BL_segments <- c(BaselinePieces[1], diff(BaselinePieces))
  indx <- min(which(tt <= BaselinePieces))
  value <- BL_segments
  n_pieces <- length(BaselinePieces)
  if (indx == n_pieces) {
    value[n_pieces] <- tt - BaselinePieces[n_pieces - 1]
  } else if (indx > 1) {
    value[(indx + 1):n_pieces] <- 0
    value[indx] <- tt - BaselinePieces[indx - 1]
  } else {
    value[(indx + 1):n_pieces] <- 0
    value[indx] <- tt
  }
  ## return
  value
}

logL_heart <- function(par, data, BaselinePieces) {
  npieces <- length(BaselinePieces)
  if (BaselinePieces[npieces] < max(data$Time)) {
    BaselinePieces[npieces] <- max(data$Time) + 1e-08
    warning("Extend the Baseline Pieces to adjust the data")
  }
  nbeta <- ncol(data) - 3
  ## par = \THETA in the paper
  par_beta <- par[1:nbeta]
  par_theta <- par[nbeta + 1]
  par_alpha <- par[(nbeta + 2):length(par)]
  m <- length(unique(data$ID))
  expXBeta <- exp(as.matrix(data[, 4:(3 + nbeta)]) %*% (as.matrix(par_beta)))
  rho_0_ij <- rho_0(par_BaselinePW = par_alpha, BaselinePieces = BaselinePieces, 
                    data$Time[data$Event == 1])
  rho_i <- expXBeta[data$Event == 1] * rho_0_ij
  rho_i[rho_i < 1e-100] <- 1e-100
  sum.log.rho_i <- sum(log(rho_i))
  ## these codes to make sure that the order will not change 
  ## if the patient ID is not ordered
  n_ij <- table(data$ID)[order(unique(data$ID))] - 1  
  ## if there is a subject with 0 event, 
  ## the sequence will not be generated for this subject
  theta_j_1 <- par_theta + sequence(n_ij) - 1  
  theta_j_1[theta_j_1 < 1e-100] <- 1e-100
  sum.log_theta_j_1 <- sum(log(theta_j_1))
  mu0i <- mu0(par_BaselinePW = par_alpha, BaselinePieces = BaselinePieces, 
              data$Time[data$Event == 0])
  mui <- mu0i * expXBeta[data$Event == 0]
  mui_theta <- par_theta + mui
  mui_theta[mui_theta < 1e-100] <- 1e-100
  sum.log_theta_mui <- sum((n_ij + par_theta) * log(mui_theta))
  if (par_theta < 1e-100) {
    par_theta <- 1e-100
  }
  logLH <- m * par_theta * log(par_theta) + sum.log.rho_i + 
    sum.log_theta_j_1 - sum.log_theta_mui
  penal <- ifelse(par_theta < 0 | min(par_alpha) < 0, 1e+50, 0)
  negLH <- -logLH + penal
  ## Calculate the gradient
  dl_dbeta <- apply(diag((n_ij - mui)/(par_theta + mui) * par_theta) %*% 
                      as.matrix(data[data$Event == 0, 4:(3 + nbeta)]), 2, sum)
  dl_dtheta <- m + m * log(par_theta) + 
    sum(1/(par_theta + sequence(n_ij) - 1)) - 
    sum((n_ij + par_theta)/(par_theta + mui)) - sum(log(mui_theta))
  indx <- apply(as.array(data$Time[data$Event == 1]), 1, 
                whereT, BaselinePieces)
  if (length(unique(indx)) < length(par_alpha)) {
    stop("Some segements have zero events!")
  }
  indx_taui <- apply(as.array(data$Time[data$Event == 0]), 1, 
                     whereT, BaselinePieces)
  dl_dalpha_part2 <- diag((n_ij + par_theta) / (par_theta + mui) * 
                            expXBeta[data$Event == 0]) %*% 
    t(apply(array(data$Time[data$Event == 0]), 1, dmu0_alpha, BaselinePieces))
  dl_dalpha <- 1 / par_alpha * table(indx) - apply(dl_dalpha_part2, 2, sum)
  attr(negLH, "gradient") <- -c(dl_dbeta, dl_dtheta, dl_dalpha)
  ## return
  negLH
}


## functions to export ========================================================

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
heart <- function(formula, data, BaselinePieces, ini) {
  ## record the function call 
  call <- match.call()
  # Prepare data matrix LRX
  mf <- model.frame(formula, data)
  mm <- model.matrix(formula, data)
  
  fit <- nlm(logL_heart, ini, data = data, 
             BaselinePieces = BaselinePieces, hessian = TRUE)
  length_par <- length(ini)
  length_alpha <- length(BaselinePieces)
  length_beta <- length_par - length_alpha - 1
  est_beta <- matrix(NA, nrow = length_beta, ncol = 3)
  colnames(est_beta) <- c("beta", "se(beta)", "two sided p-value")
  
  se_vec <- sqrt(diag(solve(fit$hessian)))
  
  est_beta[, 1] <- fit$estimate[1:length_beta]
  est_beta[, 2] <- se_vec[1:length_beta]
  est_beta[, 3] <- 2 * pnorm(-abs(est_beta[, 1]/est_beta[, 2]))
  
  est_theta <- matrix(NA, nrow = 1, ncol = 2)
  colnames(est_theta) <- c("theta", "se(theta)")
  est_theta[1, ] <- c(fit$estimate[length_beta + 1], se_vec[length_beta + 1])
  
  est_alpha <- matrix(NA, nrow = length_alpha, ncol = 2)
  colnames(est_alpha) <- c("alpha", "se(alpha)")
  
  est_alpha[, 1] <- fit$estimate[(length_beta + 2):length_par]
  est_alpha[, 2] <- se_vec[(length_beta + 2):length_par]
  results <- list(beta = est_beta, theta = est_theta, alpha = est_alpha, 
                  Convergence = fit$code, `Fisher Info Matrix` = fit$hessian)
  return(results)
}

