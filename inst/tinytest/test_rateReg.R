## Testing exception handling of rateReg

library(reda)

isNumMatrix <- reda:::isNumMatrix
isNumOne <- reda:::isNumOne

data(simuDat)

## error if formula is not specified
expect_error(rateReg(data = simuDat), "formula")

## compatibility of Survr
expect_warning(rateReg(Survr(ID, time, event) ~ group, simuDat),
               "deprecated")

## error if subset is not logical
expect_error(rateReg(Recur(time, ID, event) ~ group, simuDat, subset = 1),
             "subset")

## error if formula response is not of class 'Recur'
expect_error(rateReg(ID ~ group, simuDat), "formula")

## warning if some spline basis does cover any event
expect_error(
    rateReg(Recur(time, ID, event) ~ group, simuDat,
            knots = c(50, 100, 150, 170),
            control = list(Boundary.knots = c(0, 180))),
    "does not capture any event"
)

## error if verbose is not logical vector of length one
expect_error(
    rateReg(Recur(time, ID, event) ~ group, simuDat,
            control = list(verbose = 1)),
    "logical value"
)

## error if something is wrong with the starting values
expect_error(
    rateReg(Recur(time, ID, event) ~ group, simuDat,
            start = list(beta = c(0.1, 1))),
    "coefficients"
)
expect_error(
    rateReg(Recur(time, ID, event) ~ group, simuDat,
            start = list(theta = 0)),
    "frailty"
)




## Quick tests for normal usages
## try the case without any covariates
expect_equal(coef(rateReg(Recur(time, ID, event) ~ 1, simuDat)),
             numeric(0))
expect_equal(coef(rateReg(Recur(time, ID, event) ~ 1, simuDat,
                          spline = "mSplines")),
             numeric(0))

## test on subsetting
expect_equal(length(coef(
    rateReg(Recur(time, ID, event) ~ group + gender,
            simuDat, subset = ID %in% seq_len(50))
)), 2)

## test on na.action on missing values in covariates
tmpDat <- subset(simuDat, ID %in% seq_len(50))
tmpDat[6 : 8, "x1"] <- NA
expect_equal(attr(
    rateReg(Recur(time, ID, event) ~ group + x1,
            tmpDat, na.action = na.exclude), "na.action"
), "na.exclude")
expect_error(
    rateReg(Recur(time, ID, event) ~ group + x1,
            tmpDat, na.action = "na.fail"),
    "missing values"
)

## test on contrasts
expect_equal(names(coef(
    rateReg(Recur(time, ID, event) ~ x1 + group + gender,
            simuDat, ID %in% seq_len(50),
            contrasts = list(group = "contr.sum",
                             gender = "contr.poly"))
)), c("x1", "group1", "gender.L"))

## test related methods
## set up three fitted objects
testDat <- base::subset(simuDat, ID %in% seq_len(50))
constFit <- rateReg(Recur(time, ID, event) ~ x1 + group + gender,
                    testDat)
piecesFit <- rateReg(Recur(time, ID, event) ~ x1 + group + gender,
                     testDat, knots = seq.int(28, 140, 28))
splineFit <- rateReg(Recur(time, ID, event) ~ x1 + group + gender,
                     testDat, knots = c(60, 90, 120), degree = 3,
                     spline = "mSplines")
## test summary
expect_equivalent(class(summary(constFit)), "summary.rateReg")

## skip testing coef since it has been covered in previous tests

## test confint
expect_equal(isNumMatrix(confint(constFit), 3L, 2L), TRUE)
expect_equal(isNumMatrix(confint(piecesFit, parm = 1:2), 2L, 2L), TRUE)
expect_equal(isNumMatrix(confint(splineFit, parm = "x1"), 1L, 2L), TRUE)
expect_error(confint(splineFit, factor(1)), "parm")

## test AIC and BIC
expect_equal(isNumOne(AIC(constFit)), TRUE)
expect_equal(isNumOne(BIC(constFit)), TRUE)
expect_equal(dim(AIC(constFit, piecesFit, splineFit)), c(3L, 2L))
expect_equal(dim(BIC(constFit, piecesFit, splineFit)), c(3L, 2L))

## test baseRate
br_constFit <- baseRate(constFit)
expect_equivalent(class(br_constFit), "baseRate.rateReg")
expect_equivalent(class(baseRate(splineFit)), "baseRate.rateReg")
tmpFit <- splineFit
tmpFit@spline$spline <- "iSplines"
expect_error(baseRate(tmpFit), "Unknown splines type.")
## test plot,baseRate.rateReg-method
expect_equivalent(class(plot(br_constFit, conf.int = TRUE)),
                  c("gg", "ggplot"))
## trigger warnings
## set.seed(123)
## sinDat <- simEventData(100, rho = function(tVec) 1 - sin(tVec))
## sinFit <- rateReg(Recur(time, ID, event) ~ 1, sinDat,
##                   knots = c(2), degree = 3)
## expect_error(baseRate(sinFit), "variance-covariance")

## test mcf,rateReg-method
mcf_constFit <- mcf(constFit)
mcf_piecesFit <- mcf(piecesFit, newdata = rbind(NA, testDat[1, ]),
                     na.action = NULL)
mcf_splineFit <- mcf(splineFit,
                     newdata = rbind(NA, testDat[1:10, ]),
                     na.action = "na.exclude",
                     control = list(grid = seq.int(0, 168, by = 1)))
expect_equivalent(class(mcf_constFit), "mcf.rateReg")
expect_equivalent(class(mcf_piecesFit), "mcf.rateReg")
expect_equivalent(class(mcf_splineFit), "mcf.rateReg")
expect_error(mcf(splineFit, control = list(grid = factor(1:2))),
             "grid")
expect_error(mcf(splineFit, control = list(grid = seq.int(0, 200))),
             "boundary")
## test plot,mcf.rateReg-method
expect_equivalent(class(
    plot(mcf_constFit, conf.int = TRUE, lty = 2, col = "red")
), c("gg", "ggplot"))
expect_equivalent(class(
    plot(mcf_splineFit, conf.int = TRUE, lty = 1:4, col = 1:4)
), c("gg", "ggplot"))
expect_equivalent(class(plot(mcf_splineFit)), c("gg", "ggplot"))

## show methods
expect_equal(show(piecesFit), piecesFit)
expect_equal(show(splineFit), splineFit)
expect_equal(show(summary(piecesFit)), summary(piecesFit))
expect_equal(show(mcf(constFit)), mcf(constFit))
