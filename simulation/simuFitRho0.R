################################################################################
## simulation: Fitting models by spline rate function with degree from 0 to 3
## true rate function is rho0: 0.03 * exp(t / 168) + 0.01 * sin(10 * t / 168)
## true mcf is mu0t: 0.03*168*(exp(t/168)-1)-0.01*(168/10)*(cos(10*t/168)-1)
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
fitListRho0_2e2 <-
    parLapply(cl, seq(nRepeat),
              function(a) {
                  temp <- exportClass()
                  simuDat <- simuDataset(nSubject = 200,
                                         beta = c(0.5, 0.2),
                                         theta = 0.5,
                                         boundaryKnots0 = c(0, 168),
                                         rho0 = rho0)
                  ## fitting
                  tempFitList <- lapply(seq(0, 3), function (b) {
                      rateReg(Survr(ID, time, event) ~ x1 + x2,
                                       data = simuDat, degree = b,
                                       knots = seq(28, 140, by = 28))
                  })
                  names(tempFitList) <- paste("splineFit",
                                              seq(0, 3), sep = "")
                  ## return
                  tempFitList
              })
## Cluster ends
stopCluster(cl)

## nSubject = 500
cl <- makeCluster(nMC, type = "SOCK")
clusterExport(cl, as.list(ls()))
clusterEvalQ(cl, library(splines))
clusterEvalQ(cl, library(foreach))
clusterSetupRNG(cl)
fitListRho0_5e2 <-
    parLapply(cl, seq(nRepeat),
              function(a) {
                  temp <- exportClass()
                  simuDat <- simuDataset(nSubject = 500,
                                         beta = c(0.5, 0.2),
                                         theta = 0.5,
                                         boundaryKnots0 = c(0, 168),
                                         rho0 = rho0)
                  ## fitting
                  tempFitList <- lapply(seq(0, 3), function (b) {
                      rateReg(Survr(ID, time, event) ~ x1 + x2,
                                       data = simuDat, degree = b,
                                       knots = seq(28, 140, by = 28))
                  })
                  names(tempFitList) <- paste("splineFit",
                                              seq(0, 3), sep = "")
                  ## return
                  tempFitList
              })
## Cluster ends
stopCluster(cl)


## nSubject = 1000
cl <- makeCluster(nMC, type = "SOCK")
clusterExport(cl, as.list(ls()))
clusterEvalQ(cl, library(splines))
clusterEvalQ(cl, library(foreach))
clusterSetupRNG(cl)
fitListRho0_1e3 <-
    parLapply(cl, seq(nRepeat),
              function(a) {
                  temp <- exportClass()
                  simuDat <- simuDataset(nSubject = 1000,
                                         beta = c(0.5, 0.2),
                                         theta = 0.5,
                                         boundaryKnots0 = c(0, 168),
                                         rho0 = rho0)
                  ## fitting
                  tempFitList <- lapply(seq(0, 3), function (b) {
                      rateReg(Survr(ID, time, event) ~ x1 + x2,
                                       data = simuDat, degree = b,
                                       knots = seq(28, 140, by = 28))
                  })
                  names(tempFitList) <- paste("splineFit",
                                              seq(0, 3), sep = "")
                  ## return
                  tempFitList
              })
## Cluster ends
stopCluster(cl)

## save results
save(fitListRho0_2e2, fitListRho0_5e2, fitListRho0_1e3,
     file = "fitListRho0.RData")
