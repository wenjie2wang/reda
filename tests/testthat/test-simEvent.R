context("Testing simEvent and simEventData")

test_that("call reda::simEvent", {
    ## check z
    expect_error(simEvent(z = factor(c(1, 2))),
                 "z", fixed = TRUE)
    expect_error(simEvent(z = NA_real_),
                 "z", fixed = TRUE)
    ## check zCoef
    expect_error(simEvent(zCoef = factor(c(1, 2))),
                 "zCoef", fixed = TRUE)
    expect_error(simEvent(zCoef = NA_real_),
                 "zCoef", fixed = TRUE)
    ## check rho
    expect_error(simEvent(rho = factor(c(1, 2))),
                 "rho", fixed = TRUE)
    expect_error(simEvent(rho = NA_real_),
                 "rho", fixed = TRUE)
    ## check rhoCoef
    expect_error(simEvent(rhoCoef = factor(c(1, 2))),
                 "rhoCoef", fixed = TRUE)
    expect_error(simEvent(rhoCoef = NA_real_),
                 "rhoCoef", fixed = TRUE)
    ## check origin
    expect_error(simEvent(origin = factor(c(1, 2))),
                 "origin", fixed = TRUE)
    expect_error(simEvent(origin = NA_real_),
                 "origin", fixed = TRUE)
    ## check endTime
    expect_error(simEvent(endTime = factor(c(1, 2))),
                 "endTime", fixed = TRUE)
    expect_error(simEvent(endTime = NA_real_),
                 "endTime", fixed = TRUE)
    ## check frailty
    expect_error(simEvent(frailty = factor(c(1, 2))),
                 "frailty", fixed = TRUE)
    expect_error(simEvent(frailty = NA_real_),
                 "frailty", fixed = TRUE)
})

test_that("call reda::simEventData", {
    ## check z
    expect_error(simEventData(z = factor(c(1, 2))),
                 "z", fixed = TRUE)
    expect_error(simEventData(z = data.frame(z = 1)),
                 "z", fixed = TRUE)
    expect_error(simEventData(z = NA_real_),
                 "z", fixed = TRUE)
    ## check origin
    expect_error(simEventData(origin = factor(c(1, 2))),
                 "origin", fixed = TRUE)
    expect_error(simEventData(origin = data.frame(origin = 1)),
                 "origin", fixed = TRUE)
    expect_error(simEventData(origin = NA_real_),
                 "origin", fixed = TRUE)
    ## check endTime
    expect_error(simEventData(endTime = factor(c(1, 2))),
                 "endTime", fixed = TRUE)
    expect_error(simEventData(endTime = data.frame(endTime = 1)),
                 "endTime", fixed = TRUE)
    expect_error(simEventData(endTime = NA_real_),
                 "endTime", fixed = TRUE)
    ## check frailty
    expect_error(simEventData(frailty = factor(c(1, 2))),
                 "frailty", fixed = TRUE)
    expect_error(simEventData(frailty = data.frame(frailty = 1)),
                 "frailty", fixed = TRUE)
    expect_error(simEventData(frailty = NA_real_),
                 "frailty", fixed = TRUE)
})
