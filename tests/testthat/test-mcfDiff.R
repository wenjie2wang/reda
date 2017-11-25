context("Testing mcfDiff and mcfDiff.test")

test_that("Testing reda::mcfDiff", {
    data(valveSeats)
    valveSeats$group <- cut(valveSeats$ID, c(250, 400, 450))
    valveSeats$group3 <- cut(valveSeats$ID, c(250, 380, 410, 450))
    mcf0 <- mcf(Survr(ID, Days, No.) ~ group, valveSeats)
    mcf1 <- mcf(Survr(ID, Days, No.) ~ 1, valveSeats, ID < 400)
    mcf2 <- mcf(Survr(ID, Days, No.) ~ 1, valveSeats, ID >= 400)

    ## error for objects not from mcf.formula class
    expect_error(mcfDiff(with(valveSeats, Survr(ID, Days, No.) ~ group)),
                 "mcf.formula", fixed = TRUE)
    expect_error(mcfDiff(mcf1 = mcf0, mcf2 = NA), "mcf.formula", fixed = TRUE)
    expect_error(mcfDiff(mcf1 = 1, mcf2 = mcf0), "mcf.formula", fixed = TRUE)

    ## error if level is not a value between 0 and 1
    expect_error(mcfDiff(mcf0, level = - 1), "level", fixed = TRUE)
    expect_error(mcfDiff(mcf0, level = c(0.1, 0.2)), "level", fixed = TRUE)
    expect_error(mcfDiff(mcf0, level = "NA"), "level", fixed = TRUE)
    expect_error(mcfDiff(mcf0, level = NA), "level", fixed = TRUE)
    expect_error(mcfDiff(mcf0, level = NA_real_), "level", fixed = TRUE)

    ## error if mcf1 has more than two levels
    expect_error(mcfDiff(mcf(Survr(ID, Days, No.) ~ group3, valveSeats)),
                 "more than two groups", fixed = TRUE)

    ## warning if mcf1 has two levels and mcf2 is not NULL
    expect_warning(mcfDiff(mcf0, mcf1), "Only the 'mcf1'", fixed = TRUE)

    ## error if mcf1 has only one level and mcf is NULL
    expect_error(mcfDiff(mcf1), "'mcf2' cannot be missing", fixed = TRUE)

    ## error if mcf2 contains more than one group
    expect_error(mcfDiff(mcf1, mcf0), "only one group", fixed = TRUE)

    ## warning if the methods for variance estimates are different
    expect_warning(mcfDiff(mcf1,
                           mcf(Survr(ID, Days, No.) ~ 1, valveSeats,
                               ID < 400, variance = "Poisson"),
                           testVariance = "none"),
                   "not consistent", fixed = TRUE)

    ## warning if time origins are not the same
    valveSeats$orig <- ifelse(valveSeats$ID < 400, 0, 1)
    expect_warning(
        mcfDiff(mcf(Survr(ID, Days, No., orig) ~ group, data = valveSeats),
                testVariance = "none"),
        "The earliest time origins", fixed = TRUE)

    ## mcfDiff(mcf1, mcf2) should be equivalent to mcf1 - mcf2
    diff12 <- mcfDiff(mcf1, mcf2)
    expect_equal(mcf1 - mcf2, diff12)

    ## test plot,mcfDiff-method
    expect_equal(class(
        plot(mcfDiff(mcf0, testVariance = "none"), addOrigin = TRUE)
    ), c("gg", "ggplot"))

    ## show method
    expect_output(show(diff12), "Pseudo-Score Tests")

})


test_that("Testing reda::mcfDiff.test", {
    data(valveSeats)
    valveSeats$group <- cut(valveSeats$ID, c(250, 400, 450))
    valveSeats$group3 <- cut(valveSeats$ID, c(250, 380, 410, 450))
    mcf0 <- mcf(Survr(ID, Days, No.) ~ group, valveSeats)
    mcf0_noData <- mcf(Survr(ID, Days, No.) ~ group, valveSeats,
                       control = list(keep.data = FALSE))
    mcf1 <- mcf(Survr(ID, Days, No.) ~ 1, valveSeats, ID < 400)
    mcf2 <- mcf(Survr(ID, Days, No.) ~ 1, valveSeats, ID >= 400)
    mcf2_noData <- mcf(Survr(ID, Days, No.) ~ 1, valveSeats, ID >= 400,
                       control = list(keep.data = FALSE))
    emptyObj <- methods::new("mcfDiff.test")

    ## return empty object if testVariance == "none"
    expect_equal(mcfDiff(mcf0, testVariance = "none")@test, emptyObj)
    expect_equal(mcfDiff.test(mcf0, testVariance = "none"), emptyObj)

    ## warning if no process data in mcf1
    expect_warning(mcfDiff.test(mcf0_noData), "No processed data")

    ## erro if input objects are not of class mcf.formula
    expect_error(mcfDiff.test(1), "mcf.formula", fixed = TRUE)
    expect_error(mcfDiff.test(mcf1 = mcf0, mcf2 = NA),
                 "mcf.formula", fixed = TRUE)
    expect_error(mcfDiff.test(mcf1 = 1, mcf2 = mcf0),
                 "mcf.formula", fixed = TRUE)

    ## error if mcf1 has more than two levels
    expect_error(mcfDiff.test(mcf(Survr(ID, Days, No.) ~ group3, valveSeats)),
                 "more than two groups", fixed = TRUE)

    ## warning if mcf1 has two levels and mcf2 is not NULL
    expect_warning(mcfDiff.test(mcf0, mcf1), "Only the 'mcf1'", fixed = TRUE)

    ## error if mcf1 has only one level and mcf is NULL
    expect_error(mcfDiff.test(mcf1), "'mcf2' cannot be missing", fixed = TRUE)

    ## error if mcf2 contains more than one group
    expect_error(mcfDiff.test(mcf1, mcf0), "only one group", fixed = TRUE)

    ## warning if no process data in mcf2
    expect_warning(mcfDiff.test(mcf1, mcf2_noData), "No processed data")

    ## warning if time origins are not the same
    valveSeats$orig <- ifelse(valveSeats$ID < 400, 0, 1)
    expect_warning(
        mcfDiff.test(mcf(Survr(ID, Days, No., orig) ~ group, valveSeats)),
        "The earliest time origins", fixed = TRUE)

    ## try alternative option testVariance = "Poisson"
    expect_equivalent(class(mcfDiff.test(mcf1, mcf2, testVariance = "Pois")),
                      "mcfDiff.test")

    ## try extreme cases
    tmp <- mcf(Survr(ID, Days, No.) ~ 1, valveSeats, ID == 251)
    tmpTest <- mcfDiff.test(tmp, tmp)
    expect_equal(tmpTest@.Data[, 1L], tmpTest@.Data[, 2L])

    ## show method
    expect_output(show(tmpTest), "Pseudo-Score Tests")

})
