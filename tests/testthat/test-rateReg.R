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
                         testDat, knots = c(60, 90, 120), degree = 3)
    ## test summary
    expect_equivalent(class(summary(constFit)), "summary.rateReg")

    ## skip testing coef since it has been covered in previous tests

    ## test confint
    expect_output(str(confint(constFit)),
                  "num [1:3, 1:2]", fixed = TRUE)
    expect_output(str(confint(constFit, parm = 1:2)),
                  "num [1:2, 1:2]", fixed = TRUE)
    expect_output(str(confint(constFit, parm = "x1")),
                  "num [1, 1:2]", fixed = TRUE)

    ## test AIC and BIC
    expect_output(str(AIC(constFit)), "num", fixed = TRUE)
    expect_output(str(BIC(constFit)), "num", fixed = TRUE)
    expect_output(str(AIC(constFit, piecesFit, splineFit)),
                  "'data.frame':\t3 obs. of  2 variables:",
                  fixed = TRUE)
    expect_output(str(BIC(constFit, piecesFit, splineFit)),
                  "'data.frame':\t3 obs. of  2 variables:",
                  fixed = TRUE)

    ## test baseRate
    br_constFit <- baseRate(constFit)
    expect_equivalent(class(br_constFit), "baseRate.rateReg")

    ## test plot,baseRate.rateReg-method
    expect_equivalent(class(plot(br_constFit, conf.int = TRUE)),
                      c("gg", "ggplot"))

    ## test mcf,rateReg-method
    mcf_constFit <- mcf(constFit)
    mcf_splineFit <- mcf(splineFit, newdata = testDat[1:10, ])
    expect_equivalent(class(mcf_constFit), "mcf.rateReg")
    expect_equivalent(class(mcf_splineFit), "mcf.rateReg")

    ## test plot,mcf.rateReg-method
    expect_equivalent(class(plot(mcf_constFit)), c("gg", "ggplot"))
    expect_equivalent(class(plot(mcf_splineFit)), c("gg", "ggplot"))

})
