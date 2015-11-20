################################################################################
## simulation study to check the implementation of delta-method
## on computing mean cumulative function
## version controlled by git
################################################################################

### computation part ===========================================================
## setups
source("simuFun.R")
source("../R/fit.R")
source("../R/class.R")
source("../R/dataCheck.R")
source("../R/mcf.R")
source("../R/plot.R")

## test on baseline mcf
nRepeat <- 4
nMC <- detectCores()
set.seed(1216)

## nSubject = 200
cl <- makeCluster(nMC)
registerDoParallel(cl)
baseMcfConst_2e2 <- foreach(j = seq(nRepeat),
                     .packages = c("foreach", "splines"),
                     .export = c('Survr', 'rateReg_control', 'rateReg_start',
                                 'rateReg_mcf_control', 'gradDelta')) %dorng% {
                         temp <- exportClass()                        
                         simuDat <- simuDataset(nSubject = 200,
                                                beta = c(0.5, 0.3),
                                                theta = 0.5,
                                                alpha = c(0.06, 0.04, 0.05,
                                                          0.03, 0.04, 0.05),
                                                knots = seq(28, 140, by = 28),
                                                degree = 0L,
                                                boundaryKnots0 = c(0, 168))
                         simuMcf(data = simuDat)

                     }
stopCluster(cl)

## nSubject = 500


## nSubject = 1000


## save results for constant rate function


## setting 2: Underlying spline rate ===========================================


## setting 3: Underlying general rate ==========================================


## ### tryout ==================================================================
## ## setting 1: Underlying constant rate ======================================
## ## generate sample data
## simuDat <- simuDataset(nSubject = 200, beta = c(0.5, 0.3), theta = 0.5,
##                        alpha = c(0.06, 0.04, 0.05, 0.03, 0.04, 0.05),
##                        knots = seq(28, 140, by = 28),
##                        degree = 0L, boundaryKnots0 = c(0, 168))
## simuDat$x1 <- factor(simuDat$x1, levels = c(0, 1),
##                      labels = c("Treat", "Contr"))
## colnames(simuDat)[4:5] <- c("group", "x1")
## ## fitting
## piecesFit <- rateReg(Survr(ID, time, event) ~ group + x1, data = simuDat,
##                      knots = seq(28, 140, by = 28))
## splineFit <- rateReg(Survr(ID, time, event) ~ group + x1, data = simuDat,
##                      knots = c(56, 112), degree = 3)
## ## compute mcf
## gridTime <- head(seq(0, 168, length.out = 1e3), - 1)[- 1]
## constMcf <- mu0(piecesFit@estimates$alpha[, 1], Tvec = gridTime,
##                 bKnots = c(piecesFit@knots, 168), degree = 0L,
##                 boundaryKnots = c(0, 168))
## bsMat <- bs(gridTime, knots = splineFit@knots, degree = splineFit@degree,
##             intercept = TRUE, Boundary.knots = c(0, 168))
## splineMcf <- mu0(splineFit@estimates$alpha[, 1], Tvec = gridTime,
##                  bKnots = c(splineFit@knots, 168), degree = 3L,
##                  boundaryKnots = c(0, 168), bsMat_est = bsMat)
## ggDat <- data.frame("time" = rep(gridTime, times = 2),
##                     "mcf" = c(constMcf, splineMcf),
##                     "type" = gl(2, length(gridTime),
##                                 label = c("const.", "spline")))
## pdf(file = "mcfTrueConst.pdf", width = 7, height = 5)
## ggplot(data = ggDat, aes(x = time, y = mcf, color = type)) + geom_line() +
##     ggtitle("Estimated baseline MCF")
## plotMcf(mcf(piecesFit), conf.int = TRUE) + ylim(c(0, 15)) +
##     ggtitle("Estimated baseline MCF from HEART model")
## plotMcf(mcf(splineFit), conf.int = TRUE) + ylim(c(0, 15)) +
##     ggtitle("Estimated baseline MCF from splines rate function")
## dev.off()

## ## setting 2: Underlying spline rate ========================================
## ## generate sample data
## simuDat <- simuDataset(nSubject = 200, beta = c(0.5, 0.3), theta = 0.5,
##                        alpha = c(0.06, 0.04, 0.05, 0.03, 0.04, 0.05),
##                        knots = c(56, 112), degree = 3,
##                        boundaryKnots0 = c(0, 168))
## simuDat$x1 <- factor(simuDat$x1, levels = c(0, 1),
##                      labels = c("Treat", "Contr"))
## colnames(simuDat)[4:5] <- c("group", "x1")
## ## fitting
## piecesFit <- rateReg(Survr(ID, time, event) ~ group + x1, data = simuDat,
##                      knots = seq(28, 140, by = 28))
## splineFit <- rateReg(Survr(ID, time, event) ~ group + x1, data = simuDat,
##                      knots = c(56, 112), degree = 3)
## ## compute mcf
## gridTime <- head(seq(0, 168, length.out = 1e3), - 1)[- 1]
## constMcf <- mu0(piecesFit@estimates$alpha[, 1], Tvec = gridTime,
##                 bKnots = c(piecesFit@knots, 168), degree = 0L,
##                 boundaryKnots = c(0, 168))
## bsMat <- bs(gridTime, knots = splineFit@knots, degree = splineFit@degree,
##             intercept = TRUE, Boundary.knots = c(0, 168))
## splineMcf <- mu0(splineFit@estimates$alpha[, 1], Tvec = gridTime,
##                  bKnots = c(splineFit@knots, 168), degree = 3L,
##                  boundaryKnots = c(0, 168), bsMat_est = bsMat)
## ggDat <- data.frame("time" = rep(gridTime, times = 2),
##                     "mcf" = c(constMcf, splineMcf),
##                     "type" = gl(2, length(gridTime),
##                                 label = c("const.", "spline")))
## pdf(file = "mcfTrueSpline.pdf", width = 7, height = 5)
## ggplot(data = ggDat, aes(x = time, y = mcf, color = type)) + geom_line() +
##     ggtitle("Estimated baseline MCF")
## plotMcf(mcf(piecesFit), conf.int = TRUE) + ylim(c(0, 15)) +
##     ggtitle("Estimated baseline MCF from HEART model")
## plotMcf(mcf(splineFit), conf.int = TRUE) + ylim(c(0, 15)) +
##     ggtitle("Estimated baseline MCF from splines rate function")
## dev.off()

## ## setting 3: Underlying general rate =======================================
## ## generate sample data
## simuDat <- simuDataset(nSubject = 200, beta = c(0.2, 0.3),
##                        boundaryKnots0 = c(0, 168), rho0 = rho0)
## simuDat$x1 <- factor(simuDat$x1, levels = c(0, 1),
##                      labels = c("Treat", "Contr"))
## colnames(simuDat)[4:5] <- c("group", "x1")
## ## fitting
## piecesFit <- rateReg(Survr(ID, time, event) ~ group + x1, data = simuDat,
##                      knots = seq(28, 140, by = 28))
## splineFit <- rateReg(Survr(ID, time, event) ~ group + x1, data = simuDat,
##                      knots = c(56, 112), degree = 3)
## ## compute mcf
## gridTime <- head(seq(0, 168, length.out = 1e3), - 1)[- 1]
## mcf0 <- mu0t(gridTime)
## constMcf <- mu0(piecesFit@estimates$alpha[, 1], Tvec = gridTime,
##                 bKnots = c(piecesFit@knots, 168), degree = 0L,
##                 boundaryKnots = c(0, 168))
## bsMat <- bs(gridTime, knots = splineFit@knots, degree = splineFit@degree,
##             intercept = TRUE, Boundary.knots = c(0, 168))
## splineMcf <- mu0(splineFit@estimates$alpha[, 1], Tvec = gridTime,
##                  bKnots = c(splineFit@knots, 168), degree = 3L,
##                  boundaryKnots = c(0, 168), bsMat_est = bsMat)
## ggDat <- data.frame("time" = rep(gridTime, times = 3),
##                     "mcf" = c(mcf0, constMcf, splineMcf),
##                     "type" = gl(3, length(gridTime),
##                                 label = c("true", "const.", "spline")))
## pdf(file = "mcfGeneral.pdf", width = 7, height = 5)
## ggplot(data = ggDat, aes(x = time, y = mcf, color = type)) + geom_line() +
##     ggtitle("Estimated baseline MCF")
## plotMcf(mcf(piecesFit), conf.int = TRUE) +
##     ggtitle("Estimated baseline MCF from HEART model")
## plotMcf(mcf(splineFit), conf.int = TRUE) + 
##     ggtitle("Estimated baseline MCF from splines rate function")
## dev.off()
