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


    
