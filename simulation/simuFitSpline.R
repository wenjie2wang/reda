################################################################################
## simulation: Fitting models where true model is cubic spline rate
## version controlled by git
################################################################################

rm(list = ls())
## apply package snow for parallel computing and
## package doRNG to ensure the reproducibility
source("simuFun.R")
source("../R/fit.R")
source("../R/class.R")
source("../R/dataCheck.R")
source("../R/mcf.R")
source("../R/plot.R")

### setting 2: True model is cubic spline rate =================================

## setups
nRepeat <- 1e3
nMC <- detectCores()
set.seed(1216)

## nSubject = 200
cl <- makeCluster(nMC, type = "SOCK")
clusterExport(cl, as.list(ls()))
clusterEvalQ(cl, library(splines))
clusterEvalQ(cl, library(foreach))
clusterSetupRNG(cl)
fitListSpline_2e2 <-
    parLapply(cl, seq(nRepeat),
              function(a) {
                  temp <- exportClass()
                  simuDat <- simuDataset(nSubject = 200,
                                         beta = c(0.5, 0.2),
                                         alpha = c(0.06, 0.04, 0.05,
                                                   0.03, 0.04, 0.05),
                                         knots = c(56, 112),
                                         degree = 3L,
                                         theta = 0.5,
                                         boundaryKnots0 = c(0, 168))
                  ## fitting
                  constFit <- rateReg(Survr(ID, time, event) ~ x1 + x2,
                                       data = simuDat,
                                       knots = seq(28, 140, by = 28))
                  splineFit <- rateReg(Survr(ID, time, event) ~ x1 + x2,
                                       data = simuDat,
                                       knots = c(56, 112), degree = 3)
                  ## return
                  list(constFit = constFit, splineFit = splineFit)
              })
## Cluster ends
stopCluster(cl)

## nSubject = 500
cl <- makeCluster(nMC, type = "SOCK")
clusterExport(cl, as.list(ls()))
clusterEvalQ(cl, library(splines))
clusterEvalQ(cl, library(foreach))
clusterSetupRNG(cl)
fitListSpline_5e2 <-
    parLapply(cl, seq(nRepeat),
              function(a) {
                  temp <- exportClass()
                  simuDat <- simuDataset(nSubject = 500,
                                         beta = c(0.5, 0.2),
                                         alpha = c(0.06, 0.04, 0.05,
                                                   0.03, 0.04, 0.05),
                                         knots = c(56, 112),
                                         degree = 3L,
                                         theta = 0.5,
                                         boundaryKnots0 = c(0, 168))
                  ## fitting
                  constFit <- rateReg(Survr(ID, time, event) ~ x1 + x2,
                                       data = simuDat,
                                       knots = seq(28, 140, by = 28))
                  splineFit <- rateReg(Survr(ID, time, event) ~ x1 + x2,
                                       data = simuDat,
                                       knots = c(56, 112), degree = 3)
                  ## return
                  list(constFit = constFit, splineFit = splineFit)
              })
## Cluster ends
stopCluster(cl)

## nSubject = 1000
cl <- makeCluster(nMC, type = "SOCK")
clusterExport(cl, as.list(ls()))
clusterEvalQ(cl, library(splines))
clusterEvalQ(cl, library(foreach))
clusterSetupRNG(cl)
fitListSpline_1e3 <-
    parLapply(cl, seq(nRepeat),
              function(a) {
                  temp <- exportClass()
                  simuDat <- simuDataset(nSubject = 1000,
                                         beta = c(0.5, 0.2),
                                         alpha = c(0.06, 0.04, 0.05,
                                                   0.03, 0.04, 0.05),
                                         knots = c(56, 112),
                                         degree = 3L,
                                         theta = 0.5,
                                         boundaryKnots0 = c(0, 168))
                  ## fitting
                  constFit <- rateReg(Survr(ID, time, event) ~ x1 + x2,
                                       data = simuDat,
                                       knots = seq(28, 140, by = 28))
                  splineFit <- rateReg(Survr(ID, time, event) ~ x1 + x2,
                                       data = simuDat,
                                       knots = c(56, 112), degree = 3)
                  ## return
                  list(constFit = constFit, splineFit = splineFit)
              })
## Cluster ends
stopCluster(cl)

## save results
save(fitListSpline_2e2, fitListSpline_5e2, fitListSpline_1e3,
     file = "fitListSpline.RData")
