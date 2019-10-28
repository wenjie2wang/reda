library(reda)

### sample MCF
## Example 1. valve-seat data
## the default variance estimates by Lawless and Nadeau (1995) method
valveMcf0 <- mcf(Recur(Days, ID, No.) ~ 1, data = valveSeats)
plot(valveMcf0, conf.int = TRUE, mark.time = TRUE, addOrigin = TRUE) +
    ggplot2::xlab("Days") + ggplot2::theme_bw()

## variance estimates following Poisson process model
valveMcf1 <- mcf(Recur(Days, ID, No.) ~ 1,
                 data = valveSeats, variance = "Poisson")
## variance estimates by bootstrap method (with 1,000 bootstrap samples)
set.seed(123)
valveMcf2 <- mcf(Recur(Days, ID, No.) ~ 1,
                 data = valveSeats, variance = "bootstrap",
                 control = list(B = 1e3))

## comparing the variance estimates from different methods
library(ggplot2)
ciDat <- rbind(cbind(valveMcf0@MCF, Method = "Lawless & Nadeau"),
               cbind(valveMcf1@MCF, Method = "Poisson"),
               cbind(valveMcf2@MCF, Method = "Bootstrap"))
ggplot(ciDat, aes(x = time, y = se)) +
    geom_step(aes(color = Method, linetype = Method)) +
    xlab("Days") + ylab("SE estimates") + theme_bw()

## comparing the confidence interval estimates from different methods
ggplot(ciDat, aes(x = time)) +
    geom_step(aes(y = MCF)) +
    geom_step(aes(y = lower, color = Method, linetype = Method)) +
    geom_step(aes(y = upper, color = Method, linetype = Method)) +
    xlab("Days") + ylab("Confidence intervals") + theme_bw()


## Example 2. the simulated data
simuMcf <- mcf(Recur(time, ID, event) ~ group + gender,
               data = simuDat, ID %in% 1 : 50)
plot(simuMcf, conf.int = TRUE, lty = 1 : 4,
     legendName = "Treatment & Gender")

### estimate MCF difference between two groups
## one sample MCF object of two groups
mcf0 <- mcf(Recur(time, ID, event) ~ group, data = simuDat)
## two-sample pseudo-score tests
mcfDiff.test(mcf0)
## difference estimates over time
mcf0_diff <- mcfDiff(mcf0, testVariance = "none")
plot(mcf0_diff)

## or explicitly ask for the difference of two sample MCF
mcf1 <- mcf(Recur(time, ID, event) ~ 1, data = simuDat,
            subset = group %in% "Contr")
mcf2 <- mcf(Recur(time, ID, event) ~ 1, data = simuDat,
            subset = group %in% "Treat")
## perform two-sample tests and estimate difference at the same time
mcf12_diff1 <- mcfDiff(mcf1, mcf2)
mcf12_diff2 <- mcf1 - mcf2   # or equivalently using the `-` method
stopifnot(all.equal(mcf12_diff1, mcf12_diff2))
mcf12_diff1
plot(mcf12_diff1)

### For estimated MCF from a fitted model,
### see examples given in function rateReg.
