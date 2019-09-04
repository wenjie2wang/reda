library(reda)

data(valveSeats)
valveSeats$group <- cut(valveSeats$ID, c(250, 400, 450))
valveSeats$group3 <- cut(valveSeats$ID, c(250, 380, 410, 450))
mcf0 <- mcf(Recur(Days, ID, No.) ~ group, valveSeats)
mcf1 <- mcf(Recur(Days, ID, No.) ~ 1, valveSeats, ID < 400)
mcf2 <- mcf(Recur(Days, ID, No.) ~ 1, valveSeats, ID >= 400)

## error for objects not from mcf.formula class
expect_error(mcfDiff(with(valveSeats, Recur(Days, ID, No.) ~ group)),
             "mcf.formula")
expect_error(mcfDiff(mcf1 = mcf0, mcf2 = NA), "mcf.formula")
expect_error(mcfDiff(mcf1 = 1, mcf2 = mcf0), "mcf.formula")

## error if level is not a value between 0 and 1
expect_error(mcfDiff(mcf0, level = - 1), "level")
expect_error(mcfDiff(mcf0, level = c(0.1, 0.2)), "level")
expect_error(mcfDiff(mcf0, level = "NA"), "level")
expect_error(mcfDiff(mcf0, level = NA), "level")
expect_error(mcfDiff(mcf0, level = NA_real_), "level")

## error if mcf1 has more than two levels
expect_error(mcfDiff(mcf(Recur(Days, ID, No.) ~ group3, valveSeats)),
             "more than two groups")

## warning if mcf1 has two levels and mcf2 is not NULL
expect_warning(mcfDiff(mcf0, mcf1), "Only the 'mcf1'")

## error if mcf1 has only one level and mcf is NULL
expect_error(mcfDiff(mcf1), "'mcf2' cannot be missing")

## error if mcf2 contains more than one group
expect_error(mcfDiff(mcf1, mcf0), "only one group")

## warning if the methods for variance estimates are different
expect_warning(mcfDiff(mcf1,
                       mcf(Recur(Days, ID, No.) ~ 1, valveSeats,
                           ID < 400, variance = "Poisson"),
                       testVariance = "none"),
               "not consistent")

## warning if time origins are not the same
valveSeats$orig <- ifelse(valveSeats$ID < 400, 0, 1)
expect_warning(
    mcfDiff(mcf(Recur(Days, ID, No., origin = orig) ~ group,
                data = valveSeats), testVariance = "none"),
    "The earliest time origins")

## mcfDiff(mcf1, mcf2) should be equivalent to mcf1 - mcf2
diff12 <- mcfDiff(mcf1, mcf2)
expect_equal(mcf1 - mcf2, diff12)

## test plot,mcfDiff-method
expect_equal(class(
    plot(mcfDiff(mcf0, testVariance = "none"), addOrigin = TRUE)
), c("gg", "ggplot"))

## show method
expect_equal(show(diff12), diff12)


## test mcfDiff.test
mcf0_noData <- mcf(Recur(Days, ID, No.) ~ group, valveSeats,
                   control = list(keep.data = FALSE))
mcf2_noData <- mcf(Recur(Days, ID, No.) ~ 1, valveSeats, ID >= 400,
                   control = list(keep.data = FALSE))
emptyObj <- methods::new("mcfDiff.test")

## return empty object if testVariance == "none"
expect_equal(mcfDiff(mcf0, testVariance = "none")@test, emptyObj)
expect_equal(mcfDiff.test(mcf0, testVariance = "none"), emptyObj)

## warning if no process data in mcf1
expect_warning(mcfDiff.test(mcf0_noData), "No processed data")

## erro if input objects are not of class mcf.formula
expect_error(mcfDiff.test(1), "mcf.formula")
expect_error(mcfDiff.test(mcf1 = mcf0, mcf2 = NA),
             "mcf.formula")
expect_error(mcfDiff.test(mcf1 = 1, mcf2 = mcf0),
             "mcf.formula")

## error if mcf1 has more than two levels
expect_error(mcfDiff.test(mcf(Recur(Days, ID, No.) ~ group3, valveSeats)),
             "more than two groups")

## warning if mcf1 has two levels and mcf2 is not NULL
expect_warning(mcfDiff.test(mcf0, mcf1), "Only the 'mcf1'")

## error if mcf1 has only one level and mcf is NULL
expect_error(mcfDiff.test(mcf1), "'mcf2' cannot be missing")

## error if mcf2 contains more than one group
expect_error(mcfDiff.test(mcf1, mcf0), "only one group")

## warning if no process data in mcf2
expect_warning(mcfDiff.test(mcf1, mcf2_noData), "No processed data")

## warning if time origins are not the same
valveSeats$orig <- ifelse(valveSeats$ID < 400, 0, 1)
expect_warning(
    mcfDiff.test(mcf(Recur(Days, ID, No., origin = orig) ~ group,
                     data = valveSeats)),
    "The earliest time origins")

## try alternative option testVariance = "Poisson"
expect_equivalent(class(mcfDiff.test(mcf1, mcf2, testVariance = "Pois")),
                  "mcfDiff.test")

## try extreme cases
tmp <- mcf(Recur(Days, ID, No.) ~ 1, valveSeats, ID == 251)
tmpTest <- mcfDiff.test(tmp, tmp)
expect_equal(tmpTest@.Data[, 1L], tmpTest@.Data[, 2L])

## show method
expect_equal(show(tmpTest), tmpTest)
