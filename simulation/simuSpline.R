################################################################################
## simulation study on performance of spline rate function
## 200 patients with follow-up period: 24 * 7 = 168 days.
## covariate: treatment group, factor with level: treatment and control
## set five internal knots due to 6 visits
################################################################################

### computation part ===========================================================
## setups
source("simuFun.R")
source("../R/fit.R")
source("../R/class.R")
## source("../R/show.R")
if (! require(foreach)) {install.packages("foreach"); library(foreach)}
if (! require(doParallel)) {install.packages("doParallel"); library(doParallel)}
if (! require(doRNG)) {install.packages("doRNG"); library(doRNG)}

## fit rateReg model to recover the truth ======================================
## setting 1: 6 pieces' piecewise constant rate function
nRepeat <- 1e3
lenPara <- 9
nMC <- detectCores()
set.seed(1216)

## nSubject = 200
nMC <- detectCores()
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

## setting 2: spline rate function with 2 internal knots and degree 3 (df = 6)
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


## underlying true spline rate function
require(splines)
bsVec <- bs(1:167, knots = c(56, 112), degree = 3, Boundary.knots = c(0, 168),
            intercept = TRUE) %*% c(0.06, 0.04, 0.05, 0.03, 0.04, 0.05)
plot(1:167, bsVec, type = "l")
dev.off()
