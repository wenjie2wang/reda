context("Testing rateReg")

test_that("Testing exception handling of rateReg", {
    data(simuDat)

    ## error if formula is not specified
    expect_error(rateReg(data = simuDat), "formula", fixed = TRUE)

    ## error if subset is not logical
    expect_error(rateReg(Survr(ID, time, event) ~ group, simuDat, subset = 1),
                 "subset", fixed = TRUE)

    ## error if formula response is not of class 'Survr'
    expect_error(rateReg(ID ~ group, simuDat), "Survr", fixed = TRUE)

    ## warning if some spline basis does cover any event
    expect_error(
        rateReg(Survr(ID, time, event) ~ group, simuDat,
                knots = c(50, 100, 150, 170),
                control = list(Boundary.knots = c(0, 180))),
        "does not capture any event", fixed = TRUE
    )

    ## error if verbose is not logical vector of length one
    expect_error(
        rateReg(Survr(ID, time, event) ~ group, simuDat,
                control = list(verbose = 1)),
        "logical value", fixed = TRUE
    )

    ## error if something is wrong with the starting values
    expect_error(
        rateReg(Survr(ID, time, event) ~ group, simuDat,
                start = list(beta = c(0.1, 1))),
        "coefficients", fixed = TRUE
    )
    expect_error(
        rateReg(Survr(ID, time, event) ~ group, simuDat,
                start = list(theta = 0)),
        "frailty", fixed = TRUE
    )

})


test_that("Quick tests for normal usages", {
    ## try the case without any covariates
    expect_equal(coef(rateReg(Survr(ID, time, event) ~ 1, simuDat)),
                 numeric(0))
    expect_equal(coef(rateReg(Survr(ID, time, event) ~ 1, simuDat,
                              spline = "mSplines")),
                 numeric(0))

    ## test on subsetting
    expect_equal(length(coef(
        rateReg(Survr(ID, time, event) ~ group + gender,
                simuDat, subset = ID %in% seq_len(50))
    )), 2)

    ## test on na.action on missing values in covariates
    tmpDat <- subset(simuDat, ID %in% seq_len(50))
    tmpDat[6 : 8, "x1"] <- NA
    expect_equal(attr(
        rateReg(Survr(ID, time, event) ~ group + x1,
                tmpDat, na.action = na.exclude), "na.action"
    ), "na.exclude")
    expect_error(
        rateReg(Survr(ID, time, event) ~ group + x1,
                tmpDat, na.action = "na.fail"),
        "missing values", fixed = TRUE
    )

    ## test on contrasts
    expect_equal(names(coef(
        rateReg(Survr(ID, time, event) ~ x1 + group + gender,
                simuDat, ID %in% seq_len(50),
                contrasts = list(group = "contr.sum",
                                 gender = "contr.poly"))
    )), c("x1", "group1", "gender.L"))

    ## test related methods
    ## set up three fitted objects
    testDat <- base::subset(simuDat, ID %in% seq_len(50))
    constFit <- rateReg(Survr(ID, time, event) ~ x1 + group + gender,
                        testDat)
    piecesFit <- rateReg(Survr(ID, time, event) ~ x1 + group + gender,
                         testDat, knots = seq.int(28, 140, 28))
    splineFit <- rateReg(Survr(ID, time, event) ~ x1 + group + gender,
                         testDat, knots = c(60, 90, 120), degree = 3,
                         spline = "mSplines")
    ## test summary
    expect_equivalent(class(summary(constFit)), "summary.rateReg")

    ## skip testing coef since it has been covered in previous tests

    ## test confint
    expect_equal(isNumMatrix(confint(constFit), 3L, 2L), TRUE)
    expect_equal(isNumMatrix(confint(piecesFit, parm = 1:2), 2L, 2L), TRUE)
    expect_equal(isNumMatrix(confint(splineFit, parm = "x1"), 1L, 2L), TRUE)
    expect_error(confint(splineFit, factor(1)), "parm", fixed = TRUE)

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
    expect_error(baseRate(tmpFit), "Unknown splines type.", fixed = TRUE)
    ## test plot,baseRate.rateReg-method
    expect_equivalent(class(plot(br_constFit, conf.int = TRUE)),
                      c("gg", "ggplot"))
    ## trigger warnings
    ## set.seed(123)
    ## sinDat <- simEventData(100, rho = function(tVec) 1 - sin(tVec))
    ## sinFit <- rateReg(Survr(ID, time, event) ~ 1, sinDat,
    ##                   knots = c(2), degree = 3)
    ## expect_error(baseRate(sinFit), "variance-covariance", fixed = TRUE)

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
                 "grid", fixed = TRUE)
    expect_error(mcf(splineFit, control = list(grid = seq.int(0, 200))),
                 "boundary", fixed = TRUE)
    ## test plot,mcf.rateReg-method
    expect_equivalent(class(
        plot(mcf_constFit, conf.int = TRUE, lty = 2, col = "red")
    ), c("gg", "ggplot"))
    expect_equivalent(class(
        plot(mcf_splineFit, conf.int = TRUE, lty = 1:4, col = 1:4)
    ), c("gg", "ggplot"))
    expect_equivalent(class(plot(mcf_splineFit)), c("gg", "ggplot"))

    ## show methods
    expect_output(show(piecesFit), "Call:")
    expect_output(show(splineFit), "Call:")
    expect_output(show(summary(piecesFit)), "Call:")
    expect_output(show(mcf(constFit)), "MCF")

})
