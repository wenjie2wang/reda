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

    ## try the case without any covariates
    expect_equal(coef(rateReg(Survr(ID, time, event) ~ 1, simuDat)),
                 numeric(0))

    ## try a quick example
    expect_equal(length(coef(rateReg(Survr(ID, time, event) ~ group + gender,
                                     simuDat))), 2)

})
