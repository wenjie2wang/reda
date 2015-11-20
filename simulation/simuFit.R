################################################################################
## simulation study to check the implementation of delta-method
## on computing mean cumulative function
## version controlled by git
################################################################################

rm(list = ls())
### computation part ===========================================================
## test on baseline mcf for general rate function
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

## nSubject = 100
cl <- makeCluster(nMC, type = "SOCK")
clusterExport(cl, as.list(ls()))
clusterEvalQ(cl, library(splines))
clusterEvalQ(cl, library(foreach))
clusterSetupRNG(cl)
fitList_1e2 <- parLapply(cl, seq(nRepeat),
                         function(a) {
                             temp <- exportClass()
                             simuDat <- simuDataset(nSubject = 100,
                                                    beta = c(0.05, 0.1),
                                                    theta = 0.5,
                                                    boundaryKnots0 = c(0, 168),
                                                    rho0 = rho0)
                             simuFit(data = simuDat)
                         })
## Cluster ends
stopCluster(cl)

## nSubject = 200
cl <- makeCluster(nMC, type = "SOCK")
clusterExport(cl, as.list(ls()))
clusterEvalQ(cl, library(splines))
clusterEvalQ(cl, library(foreach))
clusterSetupRNG(cl)
fitList_2e2 <- parLapply(cl, seq(nRepeat),
                         function(a) {
                             temp <- exportClass()
                             simuDat <- simuDataset(nSubject = 200,
                                                    beta = c(0.05, 0.1),
                                                    theta = 0.5,
                                                    boundaryKnots0 = c(0, 168),
                                                    rho0 = rho0)
                             simuFit(data = simuDat)
                         })
## Cluster ends
stopCluster(cl)

## nSubject = 500
cl <- makeCluster(nMC, type = "SOCK")
clusterExport(cl, as.list(ls()))
clusterEvalQ(cl, library(splines))
clusterEvalQ(cl, library(foreach))
clusterSetupRNG(cl)
fitList_5e2 <- parLapply(cl, seq(nRepeat),
                         function(a) {
                             temp <- exportClass()
                             simuDat <- simuDataset(nSubject = 500,
                                                    beta = c(0.05, 0.1),
                                                    theta = 0.5,
                                                    boundaryKnots0 = c(0, 168),
                                                    rho0 = rho0)
                             simuFit(data = simuDat)
                         })
## Cluster ends
stopCluster(cl)

## nSubject = 1000
cl <- makeCluster(nMC, type = "SOCK")
clusterExport(cl, as.list(ls()))
clusterEvalQ(cl, library(splines))
clusterEvalQ(cl, library(foreach))
clusterSetupRNG(cl)
fitList_1e3 <- parLapply(cl, seq(nRepeat),
                         function(a) {
                             temp <- exportClass()
                             simuDat <- simuDataset(nSubject = 1000,
                                                    beta = c(0.05, 0.1),
                                                    theta = 0.5,
                                                    boundaryKnots0 = c(0, 168),
                                                    rho0 = rho0)
                             simuFit(data = simuDat)
                         })
## Cluster ends
stopCluster(cl)

## save results
save(fitList_1e2, fitList_2e2, fitList_5e2, fitList_1e3,
     file = "Fits.RData")
