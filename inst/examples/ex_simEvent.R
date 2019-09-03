library(reda)
set.seed(123)

### time-invariant covariates and coefficients
## one process
simEvent(z = c(0.5, 1), zCoef = c(1, 0))
simEvent(z = 1, zCoef = 0.5, recurrent = FALSE)

## simulated data
simEventData(z = c(0.5, 1), zCoef = c(1, 0), endTime = 2)
simEventData(z = cbind(rnorm(3), 1), zCoef = c(1, 0))
simEventData(z = matrix(rnorm(5)), zCoef = 0.5, recurrent = FALSE)


### time-varying covariates and time-varying coefficients
zFun <- function(time, intercept) {
    cbind(time / 10 + intercept, as.numeric(time > 1))
}
zCoefFun <- function(x, shift) {
    cbind(sqrt(x + shift), 1)
}
simEvent(z = zFun, zCoef = zCoefFun,
         arguments = list(z = list(intercept = 0.1),
                          zCoef = list(shift = 0.1)))

## same function of time for all processes
simEventData(3, z = zFun, zCoef = zCoefFun,
             arguments = list(z = list(intercept = 0.1),
                              zCoef = list(shift = 0.1)))

## same function within one process but different between processes
## use quote function in the arguments
simDat <- simEventData(3, z = zFun, zCoef = zCoefFun,
                       arguments = list(
                           z = list(intercept = quote(rnorm(1) / 10)),
                           zCoef = list(shift = 0.1)
                       ))
## check the intercept randomly generated,
## which should be the same within each ID but different between IDs.
unique(with(simDat, cbind(ID, intercept = round(X.1 - time / 10, 6))))


### non-negative time-varying baseline hazard rate function
simEvent(rho = function(timeVec) { sin(timeVec) + 1 })
simEventData(3, origin = rnorm(3), endTime = rnorm(3, 5),
             rho = function(timeVec) { sin(timeVec) + 1 })
## specify other arguments
simEvent(z = c(rnorm(1), rbinom(1, 1, 0.5)) / 10,
         rho = function(a, b) { sin(a + b) + 1 },
         arguments = list(rho = list(b = 0.5)))
simEventData(z = cbind(rnorm(3), rbinom(3, 1, 0.5)) / 10,
             rho = function(a, b) { sin(a + b) + 1 },
             arguments = list(rho = list(b = 0.5)))

## quadratic B-splines with one internal knot at "time = 1"
## (using function 'bSpline' from splines2 package)
simEvent(rho = splines2::bSpline, rhoCoef = c(0.8, 0.5, 1, 0.6),
         arguments = list(rho = list(degree = 2, knots = 1,
                                     intercept = TRUE,
                                     Boundary.knots = c(0, 3))))


### frailty effect
## Gamma distribution with mean one
simEvent(z = c(0.5, 1), zCoef = c(1, 0), frailty = rgamma,
         arguments = list(frailty = list(shape = 2, scale = 0.5)))

## lognormal with mean zero (on the log scale)
set.seed(123)
simEvent(z = c(0.5, 1), zCoef = c(1, 0), frailty = "rlnorm",
         arguments = list(frailty = list(sdlog = 1)))
## or equivalently
set.seed(123)
logNorm <- function(a) exp(rnorm(n = 1, mean = 0, sd = a))
simEvent(z = c(0.5, 1), zCoef = c(1, 0), frailty = logNorm,
         arguments = list(frailty = list(a = 1)))


### renewal process
## interarrival times following uniform distribution
rUnif <- function(n, rate, min) runif(n, min, max = 2 / rate)
simEvent(interarrival = rUnif,
         arguments = list(interarrival = list(min = 0)))

## interarrival times following Gamma distribution with scale one
set.seed(123)
simEvent(interarrival = function(n, rate) rgamma(n, shape = 1 / rate))
## or equivalently
set.seed(123)
simEvent(interarrival = function(rate) rgamma(n = 1, shape = 1 / rate))


### relative risk functioin
set.seed(123)
simEvent(relativeRisk = "linear")
## or equivalently
rriskFun <- function(z, zCoef, intercept) {
    as.numeric(z %*% zCoef) + intercept
}
set.seed(123)
simEvent(relativeRisk = rriskFun,
         arguments = list(relativeRisk = list(intercept = 1)))
