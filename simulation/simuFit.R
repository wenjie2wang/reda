################################################################################
## simulation: Fitting models by spline rate function with degree 0 and 2.
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

## load results
load("Fits.RData")

## test on mcf =================================================================
## 500 subjects ----------------------------------------------------------------
simuMcfList_5e2 <- lapply(fitList_5e2, function (subList) {
    simuMcf(piecesFit = subList[[1]], splineFit = subList[[2]],
            control = list(length.out = 1e3))
})

## piecewise constant
piecesMcfDat_5e2 <- sumrzMcf(mcfList = simuMcfList_5e2, pieces = TRUE)
knitr::kable(piecesMcfDat_5e2[c(150, 300, 450, 600, 750), ],
             format = "markdown", row.names = FALSE)
## |      time|     mcf0| meanEstMcf|  seEstMcf| meanEstSe|
## |---------:|--------:|----------:|---------:|---------:|
## |  25.22523| 1.673609|   1.701109| 0.1634993| 0.1728879|
## |  50.45045| 3.610998|   3.569245| 0.3324687| 0.3531763|
## |  75.67568| 5.185105|   5.178153| 0.4779824| 0.5086337|
## | 100.90090| 6.927818|   7.068303| 0.6487428| 0.6907854|
## | 126.12613| 9.618056|   9.709862| 0.8861581| 0.9451675|

## spline
splineMcfDat_5e2 <- sumrzMcf(mcfList = simuMcfList_5e2, pieces = FALSE)
knitr::kable(splineMcfDat_5e2[c(150, 300, 450, 600, 750), ],
             format = "markdown", row.names = FALSE)
## |      time|     mcf0| meanEstMcf|  seEstMcf| meanEstSe|
## |---------:|--------:|----------:|---------:|---------:|
## |  25.22523| 1.673609|   1.733466| 0.1658168| 0.1679909|
## |  50.45045| 3.610998|   3.606346| 0.3366631| 0.3422896|
## |  75.67568| 5.185105|   5.169873| 0.4775300| 0.4864452|
## | 100.90090| 6.927818|   6.998241| 0.6425708| 0.6548391|
## | 126.12613| 9.618056|   9.581469| 0.8748747| 0.8930209|


## 1000 subjects ---------------------------------------------------------------
simuMcfList_1e3 <- lapply(fitList_1e3, function (subList) {
    simuMcf(piecesFit = subList[[1]], splineFit = subList[[2]],
            control = list(length.out = 1e3))
})

## piecewise constant
piecesMcfDat_1e3 <- sumrzMcf(mcfList = simuMcfList_1e3, pieces = TRUE)
knitr::kable(piecesMcfDat_1e3[c(150, 300, 450, 600, 750), ],
            format = "markdown", row.names = FALSE)
## |      time|     mcf0| meanEstMcf|  seEstMcf| meanEstSe|
## |---------:|--------:|----------:|---------:|---------:|
## |  25.22523| 1.673609|   1.701010| 0.1167240| 0.1223022|
## |  50.45045| 3.610998|   3.568453| 0.2383975| 0.2498061|
## |  75.67568| 5.185105|   5.176904| 0.3425099| 0.3597620|
## | 100.90090| 6.927818|   7.066408| 0.4652128| 0.4885885|
## | 126.12613| 9.618056|   9.706862| 0.6364174| 0.6684798|

## spline
splineMcfDat_1e3 <- sumrzMcf(mcfList = simuMcfList_1e3, pieces = FALSE)
knitr::kable(splineMcfDat_1e3[c(150, 300, 450, 600, 750), ],
             format = "markdown", row.names = FALSE)


### mcf given newdata ==========================================================
newdata0 <- data.frame(x1 = 1, x2 = 0.2)

## 500 subjects ----------------------------------------------------------------
simuMcfList_5e2 <- lapply(fitList_5e2, function (subList) {
    simuMcf(piecesFit = subList[[1]], splineFit = subList[[2]],
            control = list(length.out = 1e3), beta0 = c(0.05, 0.1),
            newdata = newdata0)
})

## piecewise constant
piecesMcfDat_5e2 <- sumrzMcf(mcfList = simuMcfList_5e2, pieces = TRUE)
knitr::kable(piecesMcfDat_5e2[c(150, 300, 450, 600, 750), ],
             format = "markdown", row.names = FALSE)
## |      time|      mcf0| meanEstMcf|  seEstMcf| meanEstSe|
## |---------:|---------:|----------:|---------:|---------:|
## |  25.22523|  1.794959|   1.821238| 0.1652493| 0.1759147|
## |  50.45045|  3.872825|   3.821582| 0.3388247| 0.3583586|
## |  75.67568|  5.561068|   5.544625| 0.4905164| 0.5156378|
## | 100.90090|  7.430142|   7.568677| 0.6670038| 0.6999467|
## | 126.12613| 10.315443|  10.397333| 0.9120321| 0.9575289|

## spline
splineMcfDat_5e2 <- sumrzMcf(mcfList = simuMcfList_5e2, pieces = FALSE)
knitr::kable(splineMcfDat_5e2[c(150, 300, 450, 600, 750), ],
             format = "markdown", row.names = FALSE)
## |      time|      mcf0| meanEstMcf|  seEstMcf| meanEstSe|
## |---------:|---------:|----------:|---------:|---------:|
## |  25.22523|  1.794959|   1.855828| 0.1670828| 0.1779083|
## |  50.45045|  3.872825|   3.861266| 0.3426203| 0.3624675|
## |  75.67568|  5.561068|   5.535719| 0.4897158| 0.5150458|
## | 100.90090|  7.430142|   7.493690| 0.6611899| 0.6933169|
## | 126.12613| 10.315443|  10.259952| 0.9017246| 0.9455661|


## 1000 subjects ---------------------------------------------------------------
simuMcfList_1e3 <- lapply(fitList_1e3, function (subList) {
    simuMcf(piecesFit = subList[[1]], splineFit = subList[[2]],
            control = list(length.out = 1e3), beta0 = c(0.05, 0.1),
            newdata = newdata0)
})

## piecewise constant
piecesMcfDat_1e3 <- sumrzMcf(mcfList = simuMcfList_1e3, pieces = TRUE)
knitr::kable(piecesMcfDat_1e3[c(150, 300, 450, 600, 750), ],
             format = "markdown", row.names = FALSE)
## |      time|      mcf0| meanEstMcf|  seEstMcf| meanEstSe|
## |---------:|---------:|----------:|---------:|---------:|
## |  25.22523|  1.794959|   1.823861| 0.1167318| 0.1245864|
## |  50.45045|  3.872825|   3.826176| 0.2374220| 0.2537448|
## |  75.67568|  5.561068|   5.550947| 0.3430728| 0.3650906|
## | 100.90090|  7.430142|   7.576988| 0.4660087| 0.4955682|
## | 126.12613| 10.315443|  10.408228| 0.6371197| 0.6779052|

## spline
splineMcfDat_1e3 <- sumrzMcf(mcfList = simuMcfList_1e3, pieces = FALSE)
knitr::kable(splineMcfDat_1e3[c(150, 300, 450, 600, 750), ],
             format = "markdown", row.names = FALSE)
## |      time|      mcf0| meanEstMcf|  seEstMcf| meanEstSe|
## |---------:|---------:|----------:|---------:|---------:|
## |  25.22523|  1.794959|   1.857228| 0.1182266| 0.1259107|
## |  50.45045|  3.872825|   3.865476| 0.2400480| 0.2566265|
## |  75.67568|  5.561068|   5.541308| 0.3425140| 0.3646246|
## | 100.90090|  7.430142|   7.500691| 0.4616755| 0.4907947|
## | 126.12613| 10.315443|  10.269740| 0.6293189| 0.6693737|



