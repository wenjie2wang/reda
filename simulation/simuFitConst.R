################################################################################
## simulation: Fitting models where true model is piecewise constant rate.
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

### setting 1: True model is piecewise constant rate ===========================

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
fitListConst_2e2 <-
    parLapply(cl, seq(nRepeat),
              function(a) {
                  temp <- exportClass()
                  simuDat <- simuDataset(nSubject = 200,
                                         beta = c(0.5, 0.2),
                                         alpha = c(0.06, 0.04, 0.05,
                                                   0.03, 0.04, 0.05),
                                         knots = seq(28, 140, by =28),
                                         degree = 0L,
                                         theta = 0.5,
                                         boundaryKnots0 = c(0, 168))
                  ## fitting
                  constFit <- rateReg(Survr(ID, time, event) ~ x1 + x2,
                                       data = simuDat,
                                       knots = seq(28, 140, by = 28))
                  splineFit <- rateReg(Survr(ID, time, event) ~ x1 + x2,
                                       data = simuDat,
                                       knots = c(56, 112), degree = 3L)
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
fitListConst_5e2 <-
    parLapply(cl, seq(nRepeat),
              function(a) {
                  temp <- exportClass()
                  simuDat <- simuDataset(nSubject = 500,
                                         beta = c(0.5, 0.2),
                                         alpha = c(0.06, 0.04, 0.05,
                                                   0.03, 0.04, 0.05),
                                         knots = seq(28, 140, by =28),
                                         degree = 0L,
                                         theta = 0.5,
                                         boundaryKnots0 = c(0, 168))
                  ## fitting
                  constFit <- rateReg(Survr(ID, time, event) ~ x1 + x2,
                                       data = simuDat,
                                       knots = seq(28, 140, by = 28))
                  splineFit <- rateReg(Survr(ID, time, event) ~ x1 + x2,
                                       data = simuDat,
                                       knots = c(56, 112), degree = 3L)
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
fitListConst_1e3 <-
    parLapply(cl, seq(nRepeat),
              function(a) {
                  temp <- exportClass()
                  simuDat <- simuDataset(nSubject = 1000,
                                         beta = c(0.5, 0.2),
                                         alpha = c(0.06, 0.04, 0.05,
                                                   0.03, 0.04, 0.05),
                                         knots = seq(28, 140, by =28),
                                         degree = 0L,
                                         theta = 0.5,
                                         boundaryKnots0 = c(0, 168))
                  ## fitting
                  constFit <- rateReg(Survr(ID, time, event) ~ x1 + x2,
                                       data = simuDat,
                                       knots = seq(28, 140, by = 28))
                  splineFit <- rateReg(Survr(ID, time, event) ~ x1 + x2,
                                       data = simuDat,
                                       knots = c(56, 112), degree = 3L)
                  ## return
                  list(constFit = constFit, splineFit = splineFit)
              })
## Cluster ends
stopCluster(cl)

## save results
save(fitListConst_2e2, fitListConst_5e2, fitListConst_1e3,
     file = "fitListConst.RData")

## load results
load("fitListConst.RData")

simuMcfList_2e2 <- lapply(fitListConst_2e2, function (subList) {
    simuMcf(piecesFit = subList[[1]], splineFit = subList[[2]],
            control = list(length.out = 1e3))
})

## piecewise constant
constMcfDat_2e2 <- sumrzMcf(mcfList = simuMcfList_2e2,
                             pieces = TRUE, interval = FALSE)
knitr::kable(constMcfDat_2e2[c(150, 300, 450, 600, 750), ],
             format = "markdown", row.names = FALSE)
## |      time|     mcf0| meanEstMcf|  seEstMcf| meanEstSe|
## |---------:|--------:|----------:|---------:|---------:|
## |  25.22523| 1.673609|   1.503743| 0.2275003| 0.2491166|
## |  50.45045| 3.610998|   2.559814| 0.3815526| 0.4174935|
## |  75.67568| 5.185105|   3.760182| 0.5575898| 0.6081592|
## | 100.90090| 6.927818|   4.678453| 0.6915673| 0.7546819|
## | 126.12613| 9.618056|   5.570025| 0.8188145| 0.8965075|

constMcfDat_2e2 <- sumrzMcf(mcfList = simuMcfList_2e2,
                            pieces = TRUE, interval = TRUE)
knitr::kable(constMcfDat_2e2[c(150, 300, 450, 600, 750), ],
             format = "markdown", row.names = FALSE)


## spline
splineMcfDat_5e2 <- sumrzMcf(mcfList = simuMcfList_5e2,
                             pieces = FALSE, interval = FALSE)
knitr::kable(splineMcfDat_5e2[c(150, 300, 450, 600, 750), ],
             format = "markdown", row.names = FALSE)
## |      time|     mcf0| meanEstMcf|  seEstMcf| meanEstSe|
## |---------:|--------:|----------:|---------:|---------:|
## |  25.22523| 1.673609|   1.733466| 0.1658168| 0.1679909|
## |  50.45045| 3.610998|   3.606346| 0.3366631| 0.3422896|
## |  75.67568| 5.185105|   5.169873| 0.4775300| 0.4864452|
## | 100.90090| 6.927818|   6.998241| 0.6425708| 0.6548391|
## | 126.12613| 9.618056|   9.581469| 0.8748747| 0.8930209|

splineMcfDat_5e2 <- sumrzMcf(mcfList = simuMcfList_5e2,
                             pieces = FALSE, interval = TRUE)
knitr::kable(splineMcfDat_5e2[c(150, 300, 450, 600, 750), ],
             format = "markdown", row.names = FALSE)
## |      time|     mcf0| meanEstMcf| empirLower| empirUpper| meanLower| meanUpper|
## |---------:|--------:|----------:|----------:|----------:|---------:|---------:|
## |  25.22523| 1.673609|   1.733466|   1.425918|   2.062419|  1.404209|  2.062722|
## |  50.45045| 3.610998|   3.606346|   2.972171|   4.292111|  2.935470|  4.277221|
## |  75.67568| 5.185105|   5.169873|   4.270245|   6.133883|  4.216458|  6.123288|
## | 100.90090| 6.927818|   6.998241|   5.809591|   8.292400|  5.714780|  8.281702|
## | 126.12613| 9.618056|   9.581469|   7.952174|  11.386450|  7.831180| 11.331758|


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

piecesMcfDat_1e3 <- sumrzMcf(mcfList = simuMcfList_1e3,
                             pieces = TRUE, interval = TRUE)
knitr::kable(piecesMcfDat_1e3[c(150, 300, 450, 600, 750), ],
             format = "markdown", row.names = FALSE)
## |      time|     mcf0| meanEstMcf| empirLower| empirUpper| meanLower| meanUpper|
## |---------:|--------:|----------:|----------:|----------:|---------:|---------:|
## |  25.22523| 1.673609|   1.701010|   1.488561|   1.936994|  1.461302|  1.940717|
## |  50.45045| 3.610998|   3.568453|   3.139228|   4.049163|  3.078842|  4.058064|
## |  75.67568| 5.185105|   5.176904|   4.570538|   5.849626|  4.471783|  5.882024|
## | 100.90090| 6.927818|   7.066408|   6.246502|   7.999583|  6.108793|  8.024024|
## | 126.12613| 9.618056|   9.706862|   8.578829|  10.998980|  8.396666| 11.017059|


## spline
splineMcfDat_1e3 <- sumrzMcf(mcfList = simuMcfList_1e3, pieces = FALSE)
knitr::kable(splineMcfDat_1e3[c(150, 300, 450, 600, 750), ],
             format = "markdown", row.names = FALSE)
## |      time|     mcf0| meanEstMcf|  seEstMcf| meanEstSe|
## |---------:|--------:|----------:|---------:|---------:|
## |  25.22523| 1.673609|   1.732144| 0.1184477| 0.1187636|
## |  50.45045| 3.610998|   3.605116| 0.2411114| 0.2421041|
## |  75.67568| 5.185105|   5.167942| 0.3422455| 0.3440593|
## | 100.90090| 6.927818|   6.995259| 0.4608451| 0.4631410|
## | 126.12613| 9.618056|   9.577678| 0.6279406| 0.6316149|

splineMcfDat_1e3 <- sumrzMcf(mcfList = simuMcfList_1e3,
                             pieces = FALSE, interval = TRUE)
knitr::kable(splineMcfDat_1e3[c(150, 300, 450, 600, 750), ],
             format = "markdown", row.names = FALSE)
## |      time|     mcf0| meanEstMcf| empirLower| empirUpper| meanLower| meanUpper|
## |---------:|--------:|----------:|----------:|----------:|---------:|---------:|
## |  25.22523| 1.673609|   1.732144|   1.518767|   1.964190|  1.499371|  1.964916|
## |  50.45045| 3.610998|   3.605116|   3.165087|   4.079934|  3.130601|  4.079632|
## |  75.67568| 5.185105|   5.167942|   4.565302|   5.846254|  4.493598|  5.842285|
## | 100.90090| 6.927818|   6.995259|   6.190990|   7.912490|  6.087519|  7.902998|
## | 126.12613| 9.618056|   9.577678|   8.452229|  10.848511|  8.339736| 10.815621|


### mcf given newdata ==========================================================
newdata0 <- data.frame(x1 = 1, x2 = 0.2)

## 500 subjects ----------------------------------------------------------------
simuMcfList_5e2 <- lapply(fitList_5e2, function (subList) {
    simuMcf(piecesFit = subList[[1]], splineFit = subList[[2]],
            control = list(length.out = 1e3), beta0 = c(0.05, 0.1),
            newdata = newdata0)
})

## piecewise constant
piecesMcfDat_5e2 <- sumrzMcf(mcfList = simuMcfList_5e2,
                             pieces = TRUE, beta0 = c(0.05, 0.1))
knitr::kable(piecesMcfDat_5e2[c(150, 300, 450, 600, 750), ],
             format = "markdown", row.names = FALSE)
## |      time|      mcf0| meanEstMcf|  seEstMcf| meanEstSe|
## |---------:|---------:|----------:|---------:|---------:|
## |  25.22523|  1.794959|   1.821238| 0.1652493| 0.1759147|
## |  50.45045|  3.872825|   3.821582| 0.3388247| 0.3583586|
## |  75.67568|  5.561068|   5.544625| 0.4905164| 0.5156378|
## | 100.90090|  7.430142|   7.568677| 0.6670038| 0.6999467|
## | 126.12613| 10.315443|  10.397333| 0.9120321| 0.9575289|

piecesMcfDat_5e2 <- sumrzMcf(mcfList = simuMcfList_5e2, pieces = TRUE,
                             interval = TRUE, beta0 = c(0.05, 0.1))
knitr::kable(piecesMcfDat_5e2[c(150, 300, 450, 600, 750), ],
             format = "markdown", row.names = FALSE)
## |      time|      mcf0| meanEstMcf| empirLower| empirUpper| meanLower| meanUpper|
## |---------:|---------:|----------:|----------:|----------:|---------:|---------:|
## |  25.22523|  1.794959|   1.821238|   1.497600|   2.146683|  1.476452|  2.166024|
## |  50.45045|  3.872825|   3.821582|   3.174934|   4.519040|  3.119212|  4.523951|
## |  75.67568|  5.561068|   5.544625|   4.604660|   6.521795|  4.533994|  6.555257|
## | 100.90090|  7.430142|   7.568677|   6.305213|   8.880944|  6.196807|  8.940548|
## | 126.12613| 10.315443|  10.397333|   8.657027|  12.220107|  8.520611| 12.274055|


## spline
splineMcfDat_5e2 <- sumrzMcf(mcfList = simuMcfList_5e2, pieces = FALSE,
                             beta0 = c(0.05, 0.1))
knitr::kable(splineMcfDat_5e2[c(150, 300, 450, 600, 750), ],
             format = "markdown", row.names = FALSE)
## |      time|      mcf0| meanEstMcf|  seEstMcf| meanEstSe|
## |---------:|---------:|----------:|---------:|---------:|
## |  25.22523|  1.794959|   1.855828| 0.1670828| 0.1779083|
## |  50.45045|  3.872825|   3.861266| 0.3426203| 0.3624675|
## |  75.67568|  5.561068|   5.535719| 0.4897158| 0.5150458|
## | 100.90090|  7.430142|   7.493690| 0.6611899| 0.6933169|
## | 126.12613| 10.315443|  10.259952| 0.9017246| 0.9455661|

splineMcfDat_5e2 <- sumrzMcf(mcfList = simuMcfList_5e2, pieces = FALSE,
                             interval = TRUE, beta0 = c(0.05, 0.1))
knitr::kable(splineMcfDat_5e2[c(150, 300, 450, 600, 750), ],
             format = "markdown", row.names = FALSE)
## |      time|      mcf0| meanEstMcf| empirLower| empirUpper| meanLower| meanUpper|
## |---------:|---------:|----------:|----------:|----------:|---------:|---------:|
## |  25.22523|  1.794959|   1.855828|   1.542901|   2.191780|  1.507134|  2.204522|
## |  50.45045|  3.872825|   3.861266|   3.211315|   4.562294|  3.150843|  4.571689|
## |  75.67568|  5.561068|   5.535719|   4.592000|   6.531754|  4.526248|  6.545191|
## | 100.90090|  7.430142|   7.493690|   6.246326|   8.807911|  6.134814|  8.852567|
## | 126.12613| 10.315443|  10.259952|   8.541514|  12.079639|  8.406677| 12.113228|


## 1000 subjects ---------------------------------------------------------------
simuMcfList_1e3 <- lapply(fitList_1e3, function (subList) {
    simuMcf(piecesFit = subList[[1]], splineFit = subList[[2]],
            control = list(length.out = 1e3), beta0 = c(0.05, 0.1),
            newdata = newdata0)
})

## piecewise constant
piecesMcfDat_1e3 <- sumrzMcf(mcfList = simuMcfList_1e3, pieces = TRUE,
                             beta0 = c(0.05, 0.1))
knitr::kable(piecesMcfDat_1e3[c(150, 300, 450, 600, 750), ],
             format = "markdown", row.names = FALSE)
## |      time|      mcf0| meanEstMcf|  seEstMcf| meanEstSe|
## |---------:|---------:|----------:|---------:|---------:|
## |  25.22523|  1.794959|   1.823861| 0.1167318| 0.1245864|
## |  50.45045|  3.872825|   3.826176| 0.2374220| 0.2537448|
## |  75.67568|  5.561068|   5.550947| 0.3430728| 0.3650906|
## | 100.90090|  7.430142|   7.576988| 0.4660087| 0.4955682|
## | 126.12613| 10.315443|  10.408228| 0.6371197| 0.6779052|

piecesMcfDat_1e3 <- sumrzMcf(mcfList = simuMcfList_1e3, pieces = TRUE,
                             interval = TRUE, beta0 = c(0.05, 0.1))
knitr::kable(piecesMcfDat_1e3[c(150, 300, 450, 600, 750), ],
             format = "markdown", row.names = FALSE)
## |      time|      mcf0| meanEstMcf| empirLower| empirUpper| meanLower| meanUpper|
## |---------:|---------:|----------:|----------:|----------:|---------:|---------:|
## |  25.22523|  1.794959|   1.823861|   1.590093|   2.036326|  1.579676|  2.068046|
## |  50.45045|  3.872825|   3.826176|   3.346015|   4.280056|  3.328845|  4.323507|
## |  75.67568|  5.561068|   5.550947|   4.844747|   6.222542|  4.835382|  6.266511|
## | 100.90090|  7.430142|   7.576988|   6.638429|   8.490706|  6.605692|  8.548284|
## | 126.12613| 10.315443|  10.408228|   9.149667|  11.649883|  9.079558| 11.736898|

## spline
splineMcfDat_1e3 <- sumrzMcf(mcfList = simuMcfList_1e3, pieces = FALSE,
                             beta0 = c(0.05, 0.1))
knitr::kable(splineMcfDat_1e3[c(150, 300, 450, 600, 750), ],
             format = "markdown", row.names = FALSE)
## |      time|      mcf0| meanEstMcf|  seEstMcf| meanEstSe|
## |---------:|---------:|----------:|---------:|---------:|
## |  25.22523|  1.794959|   1.857228| 0.1182266| 0.1259107|
## |  50.45045|  3.872825|   3.865476| 0.2400480| 0.2566265|
## |  75.67568|  5.561068|   5.541308| 0.3425140| 0.3646246|
## | 100.90090|  7.430142|   7.500691| 0.4616755| 0.4907947|
## | 126.12613| 10.315443|  10.269740| 0.6293189| 0.6693737|

splineMcfDat_1e3 <- sumrzMcf(mcfList = simuMcfList_1e3, pieces = FALSE,
                             interval = TRUE, beta0 = c(0.05, 0.1))
knitr::kable(splineMcfDat_1e3[c(150, 300, 450, 600, 750), ],
             format = "markdown", row.names = FALSE)
## |      time|      mcf0| meanEstMcf| empirLower| empirUpper| meanLower| meanUpper|
## |---------:|---------:|----------:|----------:|----------:|---------:|---------:|
## |  25.22523|  1.794959|   1.857228|   1.613600|   2.084441|  1.610447|  2.104008|
## |  50.45045|  3.872825|   3.865476|   3.382275|   4.332422|  3.362498|  4.368455|
## |  75.67568|  5.561068|   5.541308|   4.827370|   6.202830|  4.826657|  6.255959|
## | 100.90090|  7.430142|   7.500691|   6.569486|   8.400238|  6.538751|  8.462630|
## | 126.12613| 10.315443|  10.269740|   9.032528|  11.499710|  8.957791| 11.581688|


