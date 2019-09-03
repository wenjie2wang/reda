library(reda)

## constant rate function
(constFit <- rateReg(Recur(time, ID, event) ~ group + x1, data = simuDat))

## six pieces' piecewise constant rate function
(piecesFit <- rateReg(Recur(time, ID, event) ~ group + x1,
                      data = simuDat, subset = ID %in% 1:50,
                      knots = seq.int(28, 140, by = 28)))

## fit rate function with cubic spline
(splineFit <- rateReg(Recur(time, ID, event) ~ group + x1, data = simuDat,
                      knots = c(56, 84, 112), degree = 3))

## more specific summary
summary(constFit)
summary(piecesFit)
summary(splineFit)

## model selection based on AIC or BIC
AIC(constFit, piecesFit, splineFit)
BIC(constFit, piecesFit, splineFit)

## estimated covariate coefficients
coef(piecesFit)
coef(splineFit)

## confidence intervals for covariate coefficients
confint(piecesFit)
confint(splineFit, "x1", 0.9)
confint(splineFit, 1, 0.975)

## estimated baseline rate function
splinesBase <- baseRate(splineFit)
plot(splinesBase, conf.int = TRUE)

## estimated baseline mean cumulative function (MCF) from a fitted model
piecesMcf <- mcf(piecesFit)
plot(piecesMcf, conf.int = TRUE, col = "blueviolet")

## estimated MCF for given new data
newDat <- data.frame(x1 = rep(0, 2), group = c("Treat", "Contr"))
splineMcf <- mcf(splineFit, newdata = newDat, groupName = "Group",
                 groupLevels = c("Treatment", "Control"))
plot(splineMcf, conf.int = TRUE, lty = c(1, 5))

## example of further customization by ggplot2
library(ggplot2)
plot(splineMcf) +
    geom_ribbon(aes(x = time, ymin = lower,
                    ymax = upper, fill = Group),
                data = splineMcf@MCF, alpha = 0.2) +
    xlab("Days")
