# ## simulation study #1
# ## 200 patients with follow-up period: 24 * 7 = 168 days.
# ## covariate: treatment group, factor with level: treatment and control
# ## fit 6 baselinepieces heart model due to 6 visits
# beta <- 0.3
# theta <- 0.5
# alpha <- 0.06
# 
# ## generate event times for each process 
# ## which depends on random seed by myseed and covariate by x
simu1_fun <- function(ID = "1", beta = 0.3, theta = 0.5, alpha = 0.06, 
                      baselinepieces = 168, x = 0, tau = 168) {
  whereT <- function(tt, baselinepieces) {
    ## return
    min(which(tt <= baselinepieces))
  }
  ## generate rate functions for a single process
  r <- rgamma(n = 1, shape = theta, scale = 1/theta)
  rho <- r * alpha * exp(crossprod(beta, x))
  ## step 1: calculate the supremum value of rate function for each process
  rho_m <- max(rho)
  ## step 2: simulate W_i every time
  W <- NULL
  timeinter <- NULL
  while (is.null(W) || sum(W) < tau) {
    W <- c(W, rexp(n = 1, rate = rho_m))
    ## step 3
    timeinter <- c(timeinter, sum(W))
  }
  timeinter <- unique(round(timeinter, digits = 0))
  timeinter <- timeinter[timeinter < tau]
  ## step 4
  tempn <- length(timeinter)
  U <- runif(n = tempn, min = 0, max = 1)
  rho_t <- apply(as.array(timeinter), 1, 
                 whereT, baselinepieces = baselinepieces)
  ind <- U <= rho_t / rho_m
  timeout <- c(timeinter[ind], tau)
  eventout <- c(rep(1, length(timeout) - 1), 0)
  x <- t(x)
  res <- data.frame(ID = ID, time = timeout, event = eventout, x)
  return(res)
}

## generate simulation data to export as example data 
set.seed(1216)
simuDat <- NULL
for (i in 1:200) {
  simuDat <- rbind(simuDat, 
                   simu1_fun(ID = i, beta = c(0.5, 0.3), alpha = 0.06,
                             x = rbind(ifelse(i <= 100, 0, 1), rnorm(1))))
}
simuDat$X2 <- round(simuDat$X2, digits = 2)
simuDat$X1 <- factor(simuDat$X1, levels = c(0, 1), labels = c("Treat", "Contr"))
colnames(simuDat)[4:5] <- c("group", "X1")
## save(simuDat, file = "data/simuDat.RData")

## simulation verification for variance of mcf
### function part
## beta should be length 2
verf1 <- function(npat = 200, beta = c(0.5, 0.3), ...) {
  nbeta <- length(beta)
  datlist <- vector(mode = "list", length = npat)
  for (i in seq(npat)) {
      datlist[[i]] <- simu1_fun(ID = i, beta = beta, alpha = 0.06,
                           x = rbind(ifelse(i <= npat/2, 0, 1), rnorm(1)))
  }
  dat <- data.frame(do.call(rbind, datlist))
  colnames(dat) <- c("ID", "time", "event", "group", "X1")
  dat$group <- factor(dat$group,
                      levels = c(0, 1), labels = c("Treat", "Contr"))
  fit <- heart(formula = Survr(ID, time, event) ~ group + X1,
               data = dat, baselinepieces = seq(28, 168, length = 6))
  est_beta <- coef(fit)
  est_se <- summary(fit)@coefficients[, 2]
  mcf_base <- mcf(fit)@MCF
  mcf_newdat <- mcf(fit, newdata = data.frame(X1 = rep(0.1, 2), 
                             group = gl(2, 1, labels = c("Treat", "Contr"))),
                    grouplevels = c("Treat", "Contr"))@MCF
  ## return
  list(est_beta = est_beta, est_se = est_se,
       mcf_base = mcf_base, mcf_newdat = mcf_newdat)
}
## nrep > 1
summerz <- function(nrep = 200, beta = c(0.5, 0.3), ...) {
  est <- matrix(NA, nrow = nrep, ncol = 4)
  MCF_base <- MCF_newdat <- NULL
  for (i in seq(nrep)) {
    temp <- verf1(beta = beta, ...)
    est[i, ] <- c(temp$est_beta, temp$est_se)
    MCF_base <- rbind(MCF_base, cbind(temp$mcf_base, simu = i))
    MCF_newdat <- rbind(MCF_newdat, cbind(temp$mcf_newdat, simu = i))
  }
  ## beta, bar(beta), se(hat(beta)), bar(hat(se(beta)))
  bars <- colMeans(est)
  bar_beta <- bars[1:2]
  bar_se_beta <- bars[3:4]
  se_hat_beta <- sqrt(diag(var(est[, 1:2])))
  out1 <- c(nrep, beta, bar_beta, se_hat_beta, bar_se_beta)
  names(out1) <- c("rep.", "group", "X1", "est_group", "est_X1",
                  "se_est_group", "se_est_X1", "bar_se_group", "bar_se_X1")
  ## return
  list(betab = out1, MCF_base = MCF_base, MCF_newdat = MCF_newdat)
}

### for estimation of beta
library(heart)
set.seed(1216)
simures50 <- summerz(50)
simures50$betab
##       rep.        group           X1    est_group       est_X1 se_est_group 
## 50.0000000    0.5000000    0.3000000    0.4856200    0.2749800    0.2000803 
##  se_est_X1 bar_se_group    bar_se_X1 
##  0.1064778    0.2055400    0.1017600 

set.seed(1216)
simures100 <- summerz(100)
simures100$betab
##       rep.        group           X1    est_group       est_X1 se_est_group 
## 100.0000000    0.5000000    0.3000000    0.4613900    0.2680900    0.1845496 
##   se_est_X1 bar_se_group    bar_se_X1 
##   0.0976210    0.2051700    0.1008600

set.seed(1216)
simures200 <- summerz(200)
simures200$betab
##         rep.        group           X1    est_group       est_X1 se_est_group 
## 200.00000000   0.50000000   0.30000000   0.45154000   0.26878000   0.18728614 
##    se_est_X1 bar_se_group    bar_se_X1 
##   0.09362135   0.20497000   0.10125000 


#### needs updating ============================================================
### baseline
 pdf("./test_mcf1.pdf", width = 7, height = 5)
plot(MCF ~ time, data = MCF_null1, type = "l", col = "red",
     xlim = c(0, 200), ylim = c(0, 12))
lines(MCF_null1$upper ~ MCF_null1$time, col = "red", lty = 3)
lines(MCF_null1$lower ~ MCF_null1$time, col = "red", lty = 3)
for (j in 1:100) {
    with(subset(res_mcf1, simu == j),
         lines(MCF ~ time, type = "l", col = "grey", lty = 2))
}
dev.off()

### MCF for newdata
library(ggplot2)
pdf("./test_mcf2.pdf", width = 7, height = 5)
ggplot(data = res_mcf2) +
    geom_line(aes(x = time, y = MCF, group = group:factor(simu), color = group), linetype = 3) +
        geom_line(data = MCF_null2, aes(x = time, y = MCF, color = group)) +
            geom_line(data = MCF_null2, aes(x = time, y = lower, color = group), linetype = 2) +
                geom_line(data = MCF_null2, aes(x = time, y = upper, color = group), linetype = 2)
dev.off()


    
