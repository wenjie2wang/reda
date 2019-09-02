## check z
expect_error(simEvent(z = factor(c(1, 2))),
             "z")
## check zCoef
expect_error(simEvent(zCoef = factor(c(1, 2))),
             "zCoef")
## check rho
expect_error(simEvent(rho = factor(c(1, 2))),
             "rho")
## check rhoCoef
expect_error(simEvent(rhoCoef = factor(c(1, 2))),
             "rhoCoef")
## check origin
expect_error(simEvent(origin = factor(c(1, 2))),
             "origin")
## check endTime
expect_error(simEvent(endTime = factor(c(1, 2))),
             "endTime")
## check frailty
expect_error(simEvent(frailty = factor(c(1, 2))),
             "frailty")

## quick examples to increase test coverage
expect_equivalent(class(
    simEvent(interarrival = stats::rexp,
             relativeRisk = function(z, zCoef) exp(as.numeric(z %*% zCoef)))
), "simEvent")
expect_equivalent({
    my_rriskFun <- function(z, zCoef) exp(as.numeric(z %*% zCoef));
    class(simEvent(interarrival = stats::rexp,
                   relativeRisk = my_rriskFun))
}, "simEvent")
expect_error(simEvent(interarrival = stats::rexp,
                      relativeRisk = "rexp"),
             "relative risk function")
expect_error(simEvent(relativeRisk = 1), "relative risk")
expect_error(simEvent(relativeRisk = c("foo", "bar")),
             "relative risk")
expect_equivalent(class(
    simEvent(z = function(tVec) as.numeric(tVec > 1),
             zCoef = function(tVec) - sin(tVec) / 10,
             rho = 0.5,
             origin = rnorm,
             end = function(n, m = 5) rnorm(n, m, sd = 0.1),
             frailty = rlnorm)
), "simEvent")
expect_equivalent(class(
    simEvent(relativeRisk = "linear")
), "simEvent")
expect_equivalent(class(
    simEvent(relativeRisk = "excess")
), "simEvent")

## use the inversion method instead of thinning method
expect_equivalent(class(
    simEvent(method = "inversion")
), "simEvent")
expect_equivalent(class(
    simEvent(method = "inversion",
             interarrival = function(n, rate) runif(n, max = 2 / rate))
), "simEvent")
expect_equivalent(class(
    simEvent(method = "inversion", endTime = 0.1,
             interarrival = function(n, rate)
                 runif(n, min = 0.1, max = 2 / rate - 0.1))
), "simEvent")
expect_warning(simEvent(rho = function(tVec) 1 / sqrt(tVec),
                        origin = 0, endTime = 1),
               "infinite")

## show method
a <- simEvent(endTime = 1)
expect_equal(show(a), a)


## test simEventData
## check z
expect_error(simEventData(z = factor(c(1, 2))),
             "z")
expect_error(simEventData(z = data.frame(z = 1)),
             "z")
expect_error(simEventData(z = c(1, NA_real_)),
             "z")
## check origin
expect_error(simEventData(origin = factor(c(1, 2))),
             "origin")
expect_error(simEventData(origin = data.frame(origin = 1)),
             "origin")
expect_error(simEventData(origin = NA_real_),
             "origin")
## check endTime
expect_error(simEventData(endTime = factor(c(1, 2))),
             "endTime")
expect_error(simEventData(endTime = data.frame(endTime = 1)),
             "endTime")
expect_error(simEventData(endTime = NA_real_),
             "endTime")
## check frailty
expect_error(simEventData(frailty = factor(c(1, 2))),
             "frailty")
expect_error(simEventData(frailty = data.frame(frailty = 1)),
             "frailty")
expect_error(simEventData(frailty = NA_real_),
             "frailty")

## survival data instead of recurrent event data
expect_equal(dim(simEventData(recurrent = FALSE)), c(1L, 5L))
