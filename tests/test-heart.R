testdat <- data.frame(ID = c(rep(1, 3), rep(2, 2), rep(3, 2), rep(4, 3)),
                      Time = c(20, 70, 100, 60, 100, 50, 70, 30, 50, 75),
                      Event = c(1, 1, 0, 1, 0, 1, 0, 1, 1, 0), 
                      X1 = c(rep(5, 3), rep(8, 2), rep(7, 2), rep(NA, 3)),
                      X2 = gl(n = 2, k = 5, length = 10))
testdat
# #    ID Time Event X1 X2
# # 1   1   20     1  5  1
# # 2   1   70     1  5  1
# # 3   1  100     0  5  1
# # 4   2   60     1  8  1
# # 5   2  100     0  8  1
# # 6   3   50     1  7  2
# # 7   3   70     0  7  2
# # 8   4   30     1 NA  2
# # 9   4   50     1 NA  2
# # 10  4   75     0 NA  2
# 
# ## test
heartfit <- heart(formula = survrec::Survr(ID, Time, Event) ~ X1 + X2, 
                  data = testdat, subset = ID %in% 1:4, 
                  na.action = na.omit, 
                  contrasts = list("contr.sum"),
                  baselinepieces = c(50, 100), 
                  start = list(beta = c(0.3, 1), 
                               theta = 0.5, 
                               alpha = c(0.15, 0.16)),
                  control = list(gradtol = 1e-6, stepmax = 1e5, 
                                 steptol = 1e-6, iterlim = 1e2))
show(heartfit) # or simply call: heartfit
str(heartfit)
summary(heartfit)
coef(heartfit)
confint(heartfit)
baseline(heartfit)
plot.baseline(heartfit)
# 
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
simu1_fun <- function(beta = 0.3, theta = 0.5, alpha = 0.06, baselinepieces = 168,
                  x = 0, myseed = 1216, tau = 168) {
  ## generate rate functions for a single process
  set.seed(myseed)
  r <- rgamma(n = 1, shape = theta, scale = 1/theta)
  rho <- r * alpha * exp(crossprod(beta, x))
  ## step 1: calculate the supremum value of rate function for each process
  rho_m <- max(rho)
  ## step 2: simulate W_i every time
  set.seed(myseed)
  W <- NULL
  timeinter <- NULL
  while (is.null(W) || sum(W) <= tau) {
  W <- c(W, rexp(n = 1, rate = rho_m))
  ## step 3
  timeinter <- c(timeinter, sum(W))
  }
  timeinter <- head(timeinter, length(timeinter) - 1)
  timeinter <- round(timeinter, digits = 0)
  ## step 4
  set.seed(myseed)
  tempn <- length(timeinter)
  U <- runif(n = tempn, min = 0, max = 1)
  rho_t <- apply(as.array(timeinter), 1, 
                 whereT, baselinepieces = baselinepieces)
  ind <- U <= rho_t / rho_m
  timeout <- c(timeinter[ind], tau)
  eventout <- c(rep(1, tempn), 0)
  res <- data.frame(ID = myseed, time = timeout, event = eventout, x)
  return(res)
}
# 
simu1_dat <- NULL
for (i in 1:200) {
  simu1_dat <- rbind(simu1_dat, simu1_fun(myseed = i, x = ifelse(i<=100, 0, 1)))
}

simu1 <- heart(formula = survrec::Survr(ID, time, event) ~ x, 
               data = simu1_dat, baselinepieces = seq(28, 168, length = 5))
show(simu1) 
summary(simu1)
coef(simu1)
confint(simu1)
baseline(simu1)
plot.baseline(simu1)
# 
# 
# 
# 
# 
# ## dataset: colon from package survrec
# # require(survrec)
# # data(colon)
# # testdat <- colon
# # BaselinePieces <- quantile(testdat$time, probs = c(0.25, 0.5, 0.75, 1))
# # attr(BaselinePieces, "names") <- NULL
# # BaselinePieces
# # [1]  21  216  880 2175
# 
# 
# 
# 
