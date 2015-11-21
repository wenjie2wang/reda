################################################################################
## simulation study on performance of spline rate function
## 200 patients with follow-up period: 24 * 7 = 168 days.
## covariate: treatment group, factor with level: treatment and control
## set five internal knots due to 6 visits
## version controlled by git
################################################################################

### computation part ===========================================================
## setups
source("simuFun.R")
source("../R/fit.R")
source("../R/class.R")
source("../R/dataCheck")

### fit rateReg model to recover the truth =====================================

## setting 1: 6 pieces' piecewise constant rate function =======================
nRepeat <- 1e3
lenPara <- 9
nMC <- detectCores()
set.seed(1216)

## nSubject = 200
cl <- makeCluster(nMC)
registerDoParallel(cl)
const_2e2 <- foreach(j = seq(nRepeat),
                     .packages = c("foreach", "splines"),
                     .export = c("Survr", "rateReg_control", "rateReg_start"),
                     .combine = "rbind") %dorng% {
                         temp <- exportClass()                        
                         simuFit(nSubject = 200)
                     }
stopCluster(cl)

## nSubject = 500
cl <- makeCluster(nMC)
registerDoParallel(cl)
const_5e2 <- foreach(j = seq(nRepeat),
                     .packages = c("foreach", "splines"),
                     .export = c("Survr", "rateReg_control", "rateReg_start"),
                     .combine = "rbind") %dorng% {
                         temp <- exportClass()                        
                         simuFit(nSubject = 500)
                     }
stopCluster(cl)

## nSubject = 1000
cl <- makeCluster(nMC)
registerDoParallel(cl)
const_1e3 <- foreach(j = seq(nRepeat),
                     .packages = c("foreach", "splines"),
                     .export = c("Survr", "rateReg_control", "rateReg_start"),
                     .combine = "rbind") %dorng% {
                         temp <- exportClass()                        
                         simuFit(nSubject = 1000)
                     }
stopCluster(cl)
## save results for constant rate function
save(const_2e2, const_5e2, const_1e3, file = "simuConst.RData")


## setting 2: cubic spline with 2 knots (df = 6) ===============================
nRepeat <- 1e3
lenPara <- 9
nMC <- detectCores()
set.seed(1216)

## nSubject = 200
cl <- makeCluster(nMC)
registerDoParallel(cl)
spline_2e2 <- foreach(j = seq(nRepeat),
                      .packages = c("foreach", "splines"),
                      .export = c("Survr", "rateReg_control", "rateReg_start"),
                      .combine = "rbind") %dorng% {
                          temp <- exportClass()                        
                          simuFit(nSubject = 200, knots0 = c(56, 112),
                                  degree = 3, boundaryKnots = c(0 ,168))
                      }
stopCluster(cl)

## nSubject = 500
cl <- makeCluster(nMC)
registerDoParallel(cl)
spline_5e2 <- foreach(j = seq(nRepeat),
                      .packages = c("foreach", "splines"),
                      .export = c("Survr", "rateReg_control", "rateReg_start"),
                      .combine = "rbind") %dorng% {
                          temp <- exportClass()                        
                          simuFit(nSubject = 500, knots0 = c(56, 112),
                                  degree = 3, boundaryKnots = c(0 ,168))
                      }
stopCluster(cl)

## nSubject = 1000
cl <- makeCluster(nMC)
registerDoParallel(cl)
spline_1e3 <- foreach(j = seq(nRepeat),
                      .packages = c("foreach", "splines"),
                      .export = c("Survr", "rateReg_control", "rateReg_start"),
                      .combine = "rbind") %dorng% {
                          temp <- exportClass()                        
                          simuFit(nSubject = 1000, knots0 = c(56, 112),
                                  degree = 3, boundaryKnots = c(0 ,168))
                      }
stopCluster(cl)
## save results for spline rate function
save(spline_2e2, spline_5e2, spline_1e3, file = "simuSpline.RData")


## setting 3: general rate function fitted by cubic spline with 5 knots, df = 9
nRepeat <- 1e3
lenPara <- 12
nMC <- detectCores()
set.seed(1216)

## nSubject = 200
cl <- makeCluster(nMC)
registerDoParallel(cl)
rho0_2e2 <- foreach(j = seq(nRepeat),
                      .packages = c("foreach", "splines"),
                      .export = c("Survr", "rateReg_control", "rateReg_start"),
                      .combine = "rbind") %dorng% {
                          temp <- exportClass()                        
                          simuFit(nSubject = 200,
                                  knots0 = seq(28, 140, by = 28),
                                  degree = 3, boundaryKnots = c(0 ,168),
                                  rho0 = rho0)
                      }
stopCluster(cl)

## nSubject = 500
cl <- makeCluster(nMC)
registerDoParallel(cl)
rho0_5e2 <- foreach(j = seq(nRepeat),
                      .packages = c("foreach", "splines"),
                      .export = c("Survr", "rateReg_control", "rateReg_start"),
                      .combine = "rbind") %dorng% {
                          temp <- exportClass()                        
                          simuFit(nSubject = 500,
                                  knots0 = seq(28, 140, by = 28),
                                  degree = 3, boundaryKnots = c(0 ,168),
                                  rho0 = rho0)
                      }
stopCluster(cl)

## nSubject = 1000
cl <- makeCluster(nMC)
registerDoParallel(cl)
rho0_1e3 <- foreach(j = seq(nRepeat),
                      .packages = c("foreach", "splines"),
                      .export = c("Survr", "rateReg_control", "rateReg_start"),
                      .combine = "rbind") %dorng% {
                          temp <- exportClass()                        
                          simuFit(nSubject = 1000,
                                  knots0 = seq(28, 140, by = 28),
                                  degree = 3, boundaryKnots = c(0 ,168),
                                  rho0 = rho0)
                      }
stopCluster(cl)

## save results for spline rate function
save(rho0_2e2, rho0_5e2, rho0_1e3, file = "simuRho0.RData")


## summary simulation results ==================================================

## 6 pieces' pieceswise constant rate function =================================
load("simuConst.RData")

## number of subjects: 200
simuSummary(const_2e2)
##        est0     barEst       seEst       barSe
## beta1  0.50 0.49220661 0.197414091 0.216742836
## beta2  0.30 0.31031885 0.104295079 0.104834380
## theta  0.50 0.51282320 0.056370478 0.055600520
## alpha1 0.06 0.06020693 0.009090975 0.009948663
## alpha2 0.04 0.04002695 0.006102372 0.006753476
## alpha3 0.05 0.05014509 0.007506700 0.008357467
## alpha4 0.03 0.03013336 0.004623664 0.005188329
## alpha5 0.04 0.04021215 0.006164543 0.006782509
## alpha6 0.05 0.05017386 0.007623767 0.008360235

## number of subjects: 500
simuSummary(const_5e2)
##        est0     barEst       seEst       barSe
## beta1  0.50 0.50109893 0.129873250 0.136348847
## beta2  0.30 0.29765965 0.063317661 0.065498842
## theta  0.50 0.50719455 0.035997125 0.034650161
## alpha1 0.06 0.06006513 0.005940667 0.006271868
## alpha2 0.04 0.04010771 0.003972226 0.004276209
## alpha3 0.05 0.05001616 0.004954101 0.005267118
## alpha4 0.03 0.03001310 0.002980693 0.003265476
## alpha5 0.04 0.04005678 0.004034287 0.004271137
## alpha6 0.05 0.04997881 0.004873124 0.005263798

## number of subjects: 1000
simuSummary(const_1e3)
##        est0     barEst       seEst       barSe
## beta1  0.50 0.49990637 0.091057826 0.096457650
## beta2  0.30 0.30278210 0.046281404 0.046132662
## theta  0.50 0.50333405 0.023828806 0.024295578
## alpha1 0.06 0.06000667 0.004039262 0.004432715
## alpha2 0.04 0.04002552 0.002700854 0.003019274
## alpha3 0.05 0.05006079 0.003395143 0.003729328
## alpha4 0.03 0.03001408 0.002091581 0.002310097
## alpha5 0.04 0.03996712 0.002770022 0.003015148
## alpha6 0.05 0.05001748 0.003365394 0.003726254


## spline rate function with 2 internal knots of degree 3 ======================
load("simuSpline.RData")

## number of subjects: 200
simuSummary(spline_2e2)
##        est0     barEst       seEst       barSe
## beta1  0.50 0.49160553 0.198331247 0.209499802
## beta2  0.30 0.31038242 0.104564136 0.104943798
## theta  0.50 0.51332248 0.057516896 0.056122484
## alpha1 0.06 0.06056272 0.010703184 0.010934975
## alpha2 0.04 0.03969611 0.008619126 0.008700001
## alpha3 0.05 0.05037548 0.009791963 0.009891794
## alpha4 0.03 0.02998315 0.007494648 0.007663052
## alpha5 0.04 0.04054036 0.008502376 0.008522231
## alpha6 0.05 0.04968507 0.009037034 0.009333009

## number of subjects: 500
simuSummary(spline_5e2)
##        est0     barEst       seEst       barSe
## beta1  0.50 0.50072403 0.130246591 0.131980551
## beta2  0.30 0.29764326 0.063811856 0.065593249
## theta  0.50 0.50715403 0.036312710 0.034934291
## alpha1 0.06 0.06014956 0.007027485 0.006873441
## alpha2 0.04 0.04003227 0.005271337 0.005509937
## alpha3 0.05 0.05020904 0.006229036 0.006230321
## alpha4 0.03 0.02982337 0.004628202 0.004816345
## alpha5 0.04 0.04024599 0.005252251 0.005354947
## alpha6 0.05 0.04969597 0.005748535 0.005901190

## number of subjects: 1000
simuSummary(spline_1e3)
##        est0     barEst       seEst       barSe
## beta1  0.50 0.49985661 0.090953517 0.093393886
## beta2  0.30 0.30267620 0.046378888 0.046201531
## theta  0.50 0.50328502 0.023885374 0.024494199
## alpha1 0.06 0.06011441 0.004697691 0.004859330
## alpha2 0.04 0.03994988 0.003763963 0.003890165
## alpha3 0.05 0.05014693 0.004254719 0.004400801
## alpha4 0.03 0.02986683 0.003399430 0.003403897
## alpha5 0.04 0.04011860 0.003679048 0.003779682
## alpha6 0.05 0.04987170 0.004036874 0.004186134

## underlying true spline rate function
bsVec <- bs(1:167, knots = c(56, 112), degree = 3, Boundary.knots = c(0, 168),
            intercept = TRUE) %*% c(0.06, 0.04, 0.05, 0.03, 0.04, 0.05)
plot(1:167, bsVec, type = "l")
dev.off()


## generate plots of fitted rate function by cubic splines =====================
load("simuRho0.RData")
pdf(file = "simuPlot.pdf", width = 7, height = 5)
plotRate(rho0_2e2) + ggtitle("200 Subjects") + ylim(c(0.02, 0.12))
plotRate(rho0_5e2) + ggtitle("500 Subjects") + ylim(c(0.02, 0.12))
plotRate(rho0_1e3) + ggtitle("1000 Subjects") + ylim(c(0.02, 0.12))
dev.off()

## undelying general rate function used
f <- function (x) { 0.03 * exp(x / 168) + 0.01 * sin(10 * x / 168)}
x <- seq(1, 168, length.out = 1000)
y <- f(x)
plot(x, y, type = "l", col = "red")

